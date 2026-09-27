//! Where Planet Nine could be tonight, and what has already been done about
//! it: Monte Carlo draws of an orbit solution, each scored by the ground
//! surveys that have searched for it and by Rubin, which will.
//!
//! **Priors.** [`PriorKind::Bb21`] is the Brown & Batygin (2021) posterior as
//! the workspace emulates it: mass, a, q and i drawn jointly by
//! `p9_core::data::posterior` and the per-object albedo U(0.2, 0.75) of the
//! reference population, with the orbit *oriented* by the (ω, Ω) of the
//! catalogued 2021 solution in `p9_survey::studies`. (The reference
//! population used for the exclusion bookkeeping draws ω and Ω uniformly,
//! which is fine for a searched *fraction* but carries no information about
//! *where* to point.) The other priors are the catalogued solutions with
//! their own spreads.
//!
//! **Ground surveys.** Each draw is pushed through the reproduced ZTF, DES
//! and Pan-STARRS1 detection models (`p9-2021-ztf`, `p9-2022-des`,
//! `p9-2024-panstarrs`). Those models know each survey's footprint, depth
//! and cadence but not the Galaxy; every survey's probability is therefore
//! multiplied by its crowding completeness (`crate::crowding`) with a mask
//! radius set by its image quality.
//!
//! **Rubin.** The `p9-2023-lsst-strategy` linking model over the wide-fast-deep
//! footprint (Dec ≤ +12°), optionally extended along the ecliptic by the
//! North Ecliptic Spur (|β| ≤ 10°, Dec ≤ +30°), with Rubin's own crowding
//! completeness in place of the hard |b| > 10° cut.

use p9_2021_ztf::detection_efficiency::detection_probability_for_orbit as ztf_probability;
use p9_2021_ztf::survey_model::ZtfSurvey;
use p9_2022_des::color_models::fiducial as des_fiducial_colors;
use p9_2022_des::survey_model::DesSurvey;
use p9_2023_lsst_strategy::{discovery_probability, LsstStrategy, P9Sample};
use p9_2024_panstarrs::detection_pipeline::detection_probability_for_orbit as ps1_probability;
use p9_2024_panstarrs::survey_model::Ps1Survey;
use p9_core::analysis::photometry::{bb21_apparent_magnitude, BB21_ALBEDO_MAX, BB21_ALBEDO_MIN};
use p9_core::constants::{GM_SUN, RAD2DEG};
use p9_core::coords::sky::{ecliptic_vec_to_equatorial_deg, equatorial_to_galactic};
use p9_core::data::posterior::{mcmc_2021_posterior, sample_from_posterior};
use p9_core::types::{elements_to_cartesian, OrbitalElements};
use p9_survey::ephemeris::nu_weight;
use p9_survey::sampling::{draw_elements, draw_orientation};
use p9_survey::schema::OrbitSolution;
use p9_survey::studies::catalog;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::Serialize;

use crate::crowding::completeness;

/// Orbit solutions a plan can be optimised for or tested against.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
pub enum PriorKind {
    /// Brown & Batygin (2021) posterior, oriented.
    Bb21,
    /// Siraj, Chyba & Tremaine (2024): closer, lighter, nearly coplanar.
    Siraj2024,
    /// Batygin et al. (2019) review best fit.
    Review2019,
    /// Batygin & Brown (2016) nominal.
    Bb2016,
}

impl PriorKind {
    pub const ALL: [PriorKind; 4] = [
        PriorKind::Bb21,
        PriorKind::Siraj2024,
        PriorKind::Review2019,
        PriorKind::Bb2016,
    ];

    pub fn label(self) -> &'static str {
        match self {
            PriorKind::Bb21 => "Brown & Batygin 2021",
            PriorKind::Siraj2024 => "Siraj, Chyba & Tremaine 2024",
            PriorKind::Review2019 => "Batygin et al. 2019",
            PriorKind::Bb2016 => "Batygin & Brown 2016",
        }
    }

    fn catalog_key(self) -> &'static str {
        match self {
            PriorKind::Bb21 => "2021 Brown & Batygin",
            PriorKind::Siraj2024 => "2024 Siraj",
            PriorKind::Review2019 => "2019 Batygin",
            PriorKind::Bb2016 => "2016 Batygin & Brown (nominal)",
        }
    }

    /// The catalogued solution behind this prior.
    pub fn solution(self) -> OrbitSolution {
        catalog()
            .into_iter()
            .find(|s| s.name.starts_with(self.catalog_key()))
            .unwrap_or_else(|| panic!("no catalogued solution for {self:?}"))
    }
}

/// Image quality of each survey, as the radius lost around a star (arcsec),
/// and the depth to which stars matter.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Seeing {
    pub mask_radius_arcsec: f64,
    pub star_depth: f64,
}

pub const ZTF_SEEING: Seeing = Seeing {
    mask_radius_arcsec: 3.0,
    star_depth: 20.5,
};
pub const PS1_SEEING: Seeing = Seeing {
    mask_radius_arcsec: 2.0,
    star_depth: 21.5,
};
pub const DES_SEEING: Seeing = Seeing {
    mask_radius_arcsec: 1.5,
    star_depth: 23.8,
};
pub const RUBIN_SEEING: Seeing = Seeing {
    mask_radius_arcsec: 1.5,
    star_depth: 24.5,
};

/// Rubin's northern reach: wide-fast-deep limit, and the North Ecliptic Spur.
pub const RUBIN_WFD_DEC_MAX_DEG: f64 = 12.0;
pub const RUBIN_NES_DEC_MAX_DEG: f64 = 30.0;
pub const RUBIN_NES_BETA_MAX_DEG: f64 = 10.0;

/// One scored draw.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Draw {
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub gal_l_deg: f64,
    pub gal_b_deg: f64,
    pub ecl_lon_deg: f64,
    pub ecl_lat_deg: f64,
    pub dist_au: f64,
    pub v_mag: f64,
    /// Weight under the Cassini-ranging orbital-phase prior (1 inside the
    /// favoured true-anomaly interval).
    pub nu_weight: f64,
    /// Probability ZTF, DES or Pan-STARRS1 would already have found it.
    pub p_ground: f64,
    /// Probability Rubin finds it: wide-fast-deep only, and with the North
    /// Ecliptic Spur.
    pub p_rubin_wfd: f64,
    pub p_rubin_nes: f64,
}

struct Surveys {
    ztf: ZtfSurvey,
    des: DesSurvey,
    ps1: Ps1Survey,
    rubin_wfd: LsstStrategy,
    rubin_nes: LsstStrategy,
}

impl Surveys {
    fn new() -> Self {
        let mut rubin_wfd = LsstStrategy::baseline();
        rubin_wfd.galactic_lat_min_deg = 0.0;
        rubin_wfd.dec_max_deg = RUBIN_WFD_DEC_MAX_DEG;
        let mut rubin_nes = rubin_wfd.clone();
        rubin_nes.dec_max_deg = RUBIN_NES_DEC_MAX_DEG;
        Self {
            ztf: ZtfSurvey::default(),
            des: DesSurvey::default(),
            ps1: Ps1Survey::default(),
            rubin_wfd,
            rubin_nes,
        }
    }
}

fn draw_orbit<R: Rng>(kind: PriorKind, sol: &OrbitSolution, rng: &mut R) -> (OrbitalElements, f64) {
    match kind {
        PriorKind::Bb21 => {
            let p = sample_from_posterior(&mcmc_2021_posterior(), rng);
            let (omega, omega_big) = draw_orientation(sol, rng);
            (
                OrbitalElements {
                    a: p.a,
                    e: p.e,
                    i: p.i,
                    omega,
                    omega_big,
                    mean_anomaly: p.mean_anomaly,
                },
                p.mass_earth,
            )
        }
        _ => (draw_elements(sol, rng), sol.mass_earth),
    }
}

fn score<R: Rng>(kind: PriorKind, sol: &OrbitSolution, s: &Surveys, rng: &mut R) -> Draw {
    let (elements, mass) = draw_orbit(kind, sol, rng);
    let albedo = rng.gen_range(BB21_ALBEDO_MIN..BB21_ALBEDO_MAX);
    let pos = elements_to_cartesian(&elements, GM_SUN).pos;
    let dist = pos.norm();
    let v = bb21_apparent_magnitude(mass, albedo, dist);
    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
    let (l, b) = equatorial_to_galactic(ra.to_radians(), dec.to_radians());
    let (l, b) = (l * RAD2DEG, b * RAD2DEG);
    let ecl_lon = pos.y.atan2(pos.x).to_degrees().rem_euclid(360.0);
    let ecl_lat = (pos.z / dist).asin().to_degrees();

    let crowd = |see: Seeing| completeness(see.star_depth, see.mask_radius_arcsec, l, b);
    let p_ztf = ztf_probability(&s.ztf, &elements, v) * crowd(ZTF_SEEING);
    let mags = des_fiducial_colors().band_magnitudes_with_albedo(mass, dist, albedo);
    let p_des = s.des.detection_probability_for_orbit(&elements, &mags) * crowd(DES_SEEING);
    let p_ps1 = ps1_probability(&s.ps1, &elements, v) * crowd(PS1_SEEING);
    let p_ground = 1.0 - (1.0 - p_ztf) * (1.0 - p_des) * (1.0 - p_ps1);

    let sample = P9Sample {
        elements,
        v_magnitude: v,
    };
    let rubin_crowd = crowd(RUBIN_SEEING);
    let p_rubin_wfd = discovery_probability(&s.rubin_wfd, &sample) * rubin_crowd;
    let p_rubin_nes = if dec > RUBIN_WFD_DEC_MAX_DEG {
        if ecl_lat.abs() <= RUBIN_NES_BETA_MAX_DEG {
            discovery_probability(&s.rubin_nes, &sample) * rubin_crowd
        } else {
            0.0
        }
    } else {
        p_rubin_wfd
    };

    // True anomaly from the conic; sin M fixes the half-plane.
    let cos_nu = (((elements.a * (1.0 - elements.e * elements.e) / dist) - 1.0) / elements.e)
        .clamp(-1.0, 1.0);
    let nu = if elements.mean_anomaly.sin() >= 0.0 {
        cos_nu.acos()
    } else {
        std::f64::consts::TAU - cos_nu.acos()
    };

    Draw {
        ra_deg: ra,
        dec_deg: dec,
        gal_l_deg: l,
        gal_b_deg: b,
        ecl_lon_deg: ecl_lon,
        ecl_lat_deg: ecl_lat,
        dist_au: dist,
        v_mag: v,
        nu_weight: nu_weight(nu.to_degrees()),
        p_ground,
        p_rubin_wfd,
        p_rubin_nes,
    }
}

/// `n` scored draws of a prior (parallel, reproducible for a given seed).
pub fn sample(kind: PriorKind, n: usize, seed: u64) -> Vec<Draw> {
    const CHUNK: usize = 2_000;
    let sol = kind.solution();
    let surveys = Surveys::new();
    (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .flat_map_iter(|k| {
            let mut rng =
                rand::rngs::StdRng::seed_from_u64(seed ^ (k as u64).wrapping_mul(0x9E37_79B9));
            let m = CHUNK.min(n - k * CHUNK);
            (0..m)
                .map(|_| score(kind, &sol, &surveys, &mut rng))
                .collect::<Vec<_>>()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_prior_has_a_catalogued_solution() {
        for k in PriorKind::ALL {
            assert!(k.solution().a_au > 100.0, "{k:?}");
        }
    }

    #[test]
    fn draws_are_reproducible_and_physical() {
        let a = sample(PriorKind::Bb21, 4_000, 7);
        let b = sample(PriorKind::Bb21, 4_000, 7);
        assert_eq!(a.len(), 4_000);
        for (x, y) in a.iter().zip(&b) {
            assert_eq!(x.ra_deg, y.ra_deg);
            assert!((0.0..=1.0).contains(&x.p_ground));
            assert!((0.0..=1.0).contains(&x.p_rubin_nes));
            assert!(x.p_rubin_nes >= x.p_rubin_wfd - 1e-12);
            assert!(x.dist_au > 50.0 && x.v_mag > 14.0 && x.v_mag < 30.0);
        }
    }

    #[test]
    fn the_oriented_prior_concentrates_toward_aphelion() {
        // ϖ ≈ 250° puts aphelion near ecliptic longitude 70°; the planet
        // dwells there, so that half of the sky holds most of the draws.
        let d = sample(PriorKind::Bb21, 20_000, 11);
        let near = d
            .iter()
            .filter(|x| ((x.ecl_lon_deg - 70.0 + 540.0).rem_euclid(360.0) - 180.0).abs() < 90.0)
            .count();
        assert!(near as f64 > 0.6 * d.len() as f64, "{near} of {}", d.len());
    }

    #[test]
    fn ground_surveys_found_the_bright_ones() {
        let d = sample(PriorKind::Bb21, 20_000, 13);
        let mean = |f: &dyn Fn(&&Draw) -> bool| {
            let sel: Vec<&Draw> = d.iter().filter(f).collect();
            sel.iter().map(|x| x.p_ground).sum::<f64>() / sel.len().max(1) as f64
        };
        let bright = mean(&|x| x.v_mag < 20.0 && x.gal_b_deg.abs() > 30.0 && x.dec_deg > -25.0);
        let faint = mean(&|x| x.v_mag > 23.0);
        let plane = mean(&|x| x.v_mag < 20.0 && x.gal_b_deg.abs() < 3.0 && x.dec_deg > -25.0);
        assert!(bright > 0.9, "bright high-latitude draws found: {bright}");
        assert!(faint < 0.2, "faint draws found: {faint}");
        assert!(plane < bright, "plane {plane} vs clean {bright}");
    }
}

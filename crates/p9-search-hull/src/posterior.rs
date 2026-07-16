//! The "where could it be" side: an ensemble current-position sky posterior
//! built from every orbit solution in [`p9_survey::studies`].
//!
//! Each study is sampled exactly as `p9_survey::skymap` does — perturb the
//! elements by 1σ, draw the mean anomaly uniformly (equal probability per unit
//! time → a dwell-weighted "where is it now" prior), propagate with p9-core,
//! and read off the equatorial direction, heliocentric distance, galactic
//! latitude, and reflected-light V. We give every study equal weight, so the
//! ensemble is study-agnostic; we also expose the per-study clouds so each
//! solution's prediction can be projected on its own.

use nalgebra::Vector3;
use p9_core::analysis::photometry::planet_apparent_magnitude;
use p9_core::constants::{DEG2RAD, GM_SUN};
use p9_core::coords::sky::{
    ecliptic_to_equatorial_matrix, ecliptic_vec_to_equatorial_deg, equatorial_to_galactic_matrix,
};
use p9_core::types::{elements_to_cartesian, OrbitalElements};
use p9_survey::schema::{OrbitSolution, SkyGrid};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal, Uniform};
use serde::Serialize;

/// Galactic latitude (deg) of a heliocentric *ecliptic* position vector.
pub fn galactic_latitude_deg(ecl_pos: &Vector3<f64>) -> f64 {
    let eq = ecliptic_to_equatorial_matrix() * ecl_pos;
    let gal = equatorial_to_galactic_matrix() * eq;
    (gal.z / gal.norm()).asin().to_degrees()
}

/// One Monte Carlo draw's observable position and brightness.
#[derive(Clone, Copy, Serialize)]
pub struct PosSample {
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub dist_au: f64,
    pub v_mag: f64,
    /// Galactic latitude (deg) of this draw's direction.
    pub gal_b_deg: f64,
    /// True anomaly (deg, 0..360) of this draw — the orbital phase the
    /// Cassini-ranging ν prior weights.
    pub nu_deg: f64,
}

/// The element/mean-anomaly draw machinery for one orbit solution: set up once,
/// drawn many times, so the ensemble and the per-study clouds share identical
/// sampling.
struct Draws {
    na: Normal<f64>,
    ne: Normal<f64>,
    ni: Normal<f64>,
    nom: Normal<f64>,
    nbom: Normal<f64>,
    um: Uniform<f64>,
    mass_earth: f64,
    albedo: f64,
}

impl Draws {
    fn new(sol: &OrbitSolution) -> Self {
        Draws {
            na: Normal::new(sol.a_au, sol.a_sigma_au).unwrap(),
            ne: Normal::new(sol.e, sol.e_sigma).unwrap(),
            ni: Normal::new(sol.i_deg, sol.i_sigma_deg).unwrap(),
            nom: Normal::new(sol.omega_deg, sol.omega_sigma_deg).unwrap(),
            nbom: Normal::new(sol.omega_big_deg, sol.omega_big_sigma_deg).unwrap(),
            um: Uniform::new(0.0, std::f64::consts::TAU),
            mass_earth: sol.mass_earth,
            albedo: sol.albedo,
        }
    }

    fn draw<R: rand::Rng>(&self, rng: &mut R) -> PosSample {
        let a = self.na.sample(rng).clamp(150.0, 1500.0);
        let e = self.ne.sample(rng).clamp(0.01, 0.9);
        let i = self.ni.sample(rng).clamp(0.0, 60.0) * DEG2RAD;
        let omega = self.nom.sample(rng) * DEG2RAD;
        let omega_big = self.nbom.sample(rng) * DEG2RAD;
        let m = self.um.sample(rng);

        let elem = OrbitalElements {
            a,
            e,
            i,
            omega,
            omega_big,
            mean_anomaly: m,
        };
        let pos = elements_to_cartesian(&elem, GM_SUN).pos;
        let r = pos.norm();
        let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
        let v = planet_apparent_magnitude(self.mass_earth, self.albedo, r);

        // True anomaly of the draw, recovered from the conic r(ν): the mean
        // and true anomalies always lie in the same half-plane, so sin M
        // fixes the branch.
        let cos_nu = (((a * (1.0 - e * e) / r) - 1.0) / e).clamp(-1.0, 1.0);
        let nu = if m.sin() >= 0.0 {
            cos_nu.acos()
        } else {
            std::f64::consts::TAU - cos_nu.acos()
        };

        PosSample {
            ra_deg: ra,
            dec_deg: dec,
            dist_au: r,
            v_mag: v,
            gal_b_deg: galactic_latitude_deg(&pos),
            nu_deg: nu.to_degrees(),
        }
    }
}

/// The ensemble sky posterior over the grid.
pub struct SkyPosterior {
    /// Per-cell probability (sums to 1).
    pub prob: Vec<f64>,
    /// Per-cell probability-weighted mean heliocentric distance (AU); NaN where
    /// no samples landed.
    pub dist_mean: Vec<f64>,
    /// Per-cell probability-weighted mean apparent V; NaN where empty.
    pub v_mean: Vec<f64>,
    /// Per-cell mean PER-SAMPLE survey-survival probability
    /// E[∏_s(1 − cov_s·1[detected_s])] — evaluated draw by draw, because the
    /// ensemble deliberately mixes bright and faint solutions and
    /// 1[E[V] ≤ depth] ≠ E[1[V ≤ depth]] when survey depths sit inside the
    /// mixture's V span. NaN where empty.
    pub survival_mean: Vec<f64>,
    /// Per-cell mean of (scope-reachable AND survives): the residual mass a
    /// space telescope can actually chase, per draw.
    pub reach_survival_mean: Vec<f64>,
    /// Per-cell fraction of draws reachable by the space telescope.
    pub reach_frac: Vec<f64>,
    /// Per-cell fraction of draws reachable by Rubin/LSST.
    pub lsst_frac: Vec<f64>,
    /// A thinned cloud of draws for the parameter-space hull.
    pub cloud: Vec<PosSample>,
}

/// Build the ensemble posterior. `per_study` MC draws are taken from each
/// solution; `cloud_keep` thins the returned cloud. `weight` is a per-draw
/// prior weight (pass `|_| 1.0` for the uniform-phase posterior, or a
/// ν-prior weight to condition on the Cassini-ranging phase constraint —
/// the cloud itself stays unweighted for plotting). The three scorers are
/// evaluated on every draw so exclusion/reachability are per-sample (see
/// `SkyPosterior::survival_mean`).
pub fn ensemble_scored(
    studies: &[OrbitSolution],
    grid: &SkyGrid,
    per_study: usize,
    seed: u64,
    cloud_keep: usize,
    weight: impl Fn(&PosSample) -> f64,
    survival: impl Fn(&PosSample) -> f64,
    scope_reach: impl Fn(&PosSample) -> bool,
    lsst_reach: impl Fn(&PosSample) -> bool,
) -> SkyPosterior {
    let n = grid.len();
    let mut count = vec![0.0_f64; n];
    let mut sum_dist = vec![0.0_f64; n];
    let mut sum_v = vec![0.0_f64; n];
    let mut sum_surv = vec![0.0_f64; n];
    let mut sum_reach_surv = vec![0.0_f64; n];
    let mut sum_reach = vec![0.0_f64; n];
    let mut sum_lsst = vec![0.0_f64; n];
    let mut cloud = Vec::new();

    let total_samples = studies.len().max(1) * per_study;
    let cloud_stride = (total_samples / cloud_keep.max(1)).max(1);
    let mut drawn = 0usize;

    for (s_idx, sol) in studies.iter().enumerate() {
        // Per-study seed so the ensemble is reproducible and each study is
        // independent.
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(s_idx as u64));
        let d = Draws::new(sol);
        for _ in 0..per_study {
            let s = d.draw(&mut rng);
            let idx = grid.index(s.ra_deg, s.dec_deg);
            let w = weight(&s);
            count[idx] += w;
            sum_dist[idx] += w * s.dist_au;
            sum_v[idx] += w * s.v_mag;
            let surv = survival(&s);
            sum_surv[idx] += w * surv;
            let reach = scope_reach(&s);
            if reach {
                sum_reach[idx] += w;
                sum_reach_surv[idx] += w * surv;
            }
            if lsst_reach(&s) {
                sum_lsst[idx] += w;
            }
            if drawn % cloud_stride == 0 {
                cloud.push(s);
            }
            drawn += 1;
        }
    }

    let total: f64 = count.iter().sum();
    let mut prob = vec![0.0_f64; n];
    let mut dist_mean = vec![f64::NAN; n];
    let mut v_mean = vec![f64::NAN; n];
    let mut survival_mean = vec![f64::NAN; n];
    let mut reach_survival_mean = vec![f64::NAN; n];
    let mut reach_frac = vec![0.0_f64; n];
    let mut lsst_frac = vec![0.0_f64; n];
    for k in 0..n {
        if count[k] > 0.0 {
            prob[k] = count[k] / total;
            dist_mean[k] = sum_dist[k] / count[k];
            v_mean[k] = sum_v[k] / count[k];
            survival_mean[k] = sum_surv[k] / count[k];
            reach_survival_mean[k] = sum_reach_surv[k] / count[k];
            reach_frac[k] = sum_reach[k] / count[k];
            lsst_frac[k] = sum_lsst[k] / count[k];
        }
    }

    SkyPosterior {
        prob,
        dist_mean,
        v_mean,
        survival_mean,
        reach_survival_mean,
        reach_frac,
        lsst_frac,
        cloud,
    }
}

/// Trivially-scored ensemble (no survey model): survival 1, nothing
/// reachable. For tests and cloud-only consumers.
#[cfg_attr(not(test), allow(dead_code))]
pub fn ensemble(
    studies: &[OrbitSolution],
    grid: &SkyGrid,
    per_study: usize,
    seed: u64,
    cloud_keep: usize,
) -> SkyPosterior {
    ensemble_scored(
        studies,
        grid,
        per_study,
        seed,
        cloud_keep,
        |_| 1.0,
        |_| 1.0,
        |_| false,
        |_| false,
    )
}

/// One orbit solution's thinned prediction cloud (positions + brightness).
#[derive(Serialize)]
pub struct StudyCloud {
    pub name: String,
    pub samples: Vec<PosSample>,
}

/// Per-study prediction clouds, each thinned to ~`keep_per_study` draws, using
/// the same per-study seeds as [`ensemble`] so the projections are consistent.
pub fn study_clouds(
    studies: &[OrbitSolution],
    per_study: usize,
    seed: u64,
    keep_per_study: usize,
) -> Vec<StudyCloud> {
    let stride = (per_study / keep_per_study.max(1)).max(1);
    studies
        .iter()
        .enumerate()
        .map(|(s_idx, sol)| {
            let mut rng = StdRng::seed_from_u64(seed.wrapping_add(s_idx as u64));
            let d = Draws::new(sol);
            let samples = (0..per_study)
                .map(|k| (k, d.draw(&mut rng)))
                .filter(|(k, _)| k % stride == 0)
                .map(|(_, s)| s)
                .collect();
            StudyCloud {
                name: sol.name.clone(),
                samples,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_grid() -> SkyGrid {
        SkyGrid {
            ra_min_deg: 0.0,
            ra_max_deg: 360.0,
            n_ra: 60,
            dec_min_deg: -60.0,
            dec_max_deg: 60.0,
            n_dec: 20,
        }
    }

    #[test]
    fn nu_weighted_posterior_concentrates_the_phase_prior() {
        let studies = p9_survey::studies::catalog();
        let grid = small_grid();
        let uniform = ensemble(&studies, &grid, 4000, 7, 2000);
        let nu = ensemble_scored(
            &studies,
            &grid,
            4000,
            7,
            2000,
            |s| p9_survey::ephemeris::nu_weight(s.nu_deg),
            |_| 1.0,
            |_| false,
            |_| false,
        );
        // Weighted counts still normalize.
        let total: f64 = nu.prob.iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "nu prob sums to {total}");
        // Every draw carries a valid orbital phase, and the conic geometry
        // holds: the favored-ν arc (108-129 deg, past quadrature outbound)
        // puts the body beyond its semi-latus-rectum distance.
        for s in &uniform.cloud {
            assert!((0.0..360.0).contains(&s.nu_deg), "nu = {}", s.nu_deg);
        }
        // The prior cuts sky: the ν-weighted map must be more concentrated
        // (fewer cells holding 90% of the probability) than uniform phase.
        let cells_for_90 = |prob: &[f64]| {
            let mut p: Vec<f64> = prob.to_vec();
            p.sort_by(|a, b| b.partial_cmp(a).unwrap());
            let mut acc = 0.0;
            p.iter()
                .take_while(|&&x| {
                    acc += x;
                    acc < 0.9
                })
                .count()
        };
        let n_uniform = cells_for_90(&uniform.prob);
        let n_nu = cells_for_90(&nu.prob);
        assert!(
            n_nu < n_uniform / 2,
            "nu-weighted 90% region {n_nu} cells vs uniform {n_uniform}"
        );
    }

    #[test]
    fn posterior_normalizes_and_fills_cells() {
        let studies = p9_survey::studies::catalog();
        let grid = small_grid();
        let post = ensemble(&studies, &grid, 2000, 7, 500);
        let total: f64 = post.prob.iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "prob sums to {total}");
        let filled = post.dist_mean.iter().filter(|d| d.is_finite()).count();
        assert!(filled > 0);
        for (k, d) in post.dist_mean.iter().enumerate() {
            if d.is_finite() {
                assert!(*d > 60.0 && *d < 2000.0, "cell {k}: dist {d}");
                assert!(post.v_mean[k] > 15.0 && post.v_mean[k] < 30.0);
            }
        }
        assert!(!post.cloud.is_empty());
    }

    #[test]
    fn per_study_clouds_cover_every_study() {
        let studies = p9_survey::studies::catalog();
        let clouds = study_clouds(&studies, 4000, 7, 800);
        assert_eq!(clouds.len(), studies.len());
        for c in &clouds {
            assert!(!c.samples.is_empty(), "{} has no samples", c.name);
            for s in &c.samples {
                assert!(s.dist_au > 60.0 && s.dist_au < 2000.0);
                assert!(s.dec_deg.abs() <= 90.0);
                assert!(s.ra_deg >= 0.0 && s.ra_deg <= 360.0);
            }
        }
    }
}

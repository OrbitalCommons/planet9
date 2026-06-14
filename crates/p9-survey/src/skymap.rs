//! Monte Carlo engine: turn an orbit solution into a current-position sky
//! probability map, a brightness distribution, and per-telescope detection
//! odds — all from p9-core primitives.
//!
//! Method. The planet's current orbital phase is unknown, so we draw the mean
//! anomaly uniformly (equal probability per unit *time* → the body spends most
//! of it near aphelion, the natural "where is it now" prior). Each sample also
//! perturbs the orbital elements by the solution's assumed 1σ. For every draw
//! we build `OrbitalElements`, propagate with `p9_core::elements_to_cartesian`
//! (heliocentric ecliptic position; at hundreds of AU the geocentric direction
//! differs by < 0.2°, the Earth-baseline parallax), convert to equatorial
//! RA/Dec and to galactic latitude with `p9_core::coords::sky`, and evaluate
//! the reflected-light apparent magnitude with `p9_core::analysis::photometry`.

use crate::schema::{OrbitSolution, Overlays, SkyGrid, StudyDetection, StudyResult, Telescope};
use nalgebra::Vector3;
use p9_core::analysis::photometry::planet_apparent_magnitude;
use p9_core::constants::{DEG2RAD, GM_SUN};
use p9_core::coords::sky::{
    ecliptic_to_equatorial_matrix, ecliptic_vec_to_equatorial_deg, equatorial_to_galactic_matrix,
    galactic_to_ecliptic_matrix,
};
use p9_core::types::{elements_to_cartesian, OrbitalElements};
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Normal, Uniform};

/// Galactic latitude (deg) of a heliocentric *ecliptic* position vector.
fn galactic_latitude_deg(ecl_pos: &Vector3<f64>) -> f64 {
    let eq = ecliptic_to_equatorial_matrix() * ecl_pos;
    let gal = equatorial_to_galactic_matrix() * eq;
    (gal.z / gal.norm()).asin().to_degrees()
}

fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let idx = ((p * (sorted.len() as f64 - 1.0)).round() as usize).min(sorted.len() - 1);
    sorted[idx]
}

/// Run the Monte Carlo for one study against the telescope catalog.
///
/// Returns the study's sky/brightness summary plus one [`StudyDetection`] per
/// telescope, in the telescopes' order.
pub fn run_study(
    sol: &OrbitSolution,
    grid: &SkyGrid,
    telescopes: &[Telescope],
    n_samples: usize,
    seed: u64,
) -> (StudyResult, Vec<StudyDetection>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let na = Normal::new(sol.a_au, sol.a_sigma_au).unwrap();
    let ne = Normal::new(sol.e, sol.e_sigma).unwrap();
    let ni = Normal::new(sol.i_deg, sol.i_sigma_deg).unwrap();
    let nom = Normal::new(sol.omega_deg, sol.omega_sigma_deg).unwrap();
    let nbom = Normal::new(sol.omega_big_deg, sol.omega_big_sigma_deg).unwrap();
    let um = Uniform::new(0.0, std::f64::consts::TAU);

    let mut prob = vec![0.0_f64; grid.len()];
    let mut dists = Vec::with_capacity(n_samples);
    let mut vmags = Vec::with_capacity(n_samples);

    // Per-telescope detection counters: (in_footprint, bright, both).
    let mut counts = vec![(0u64, 0u64, 0u64); telescopes.len()];

    for _ in 0..n_samples {
        let a = na.sample(&mut rng).clamp(150.0, 1500.0);
        let e = ne.sample(&mut rng).clamp(0.01, 0.9);
        let i = ni.sample(&mut rng).clamp(0.0, 60.0) * DEG2RAD;
        let omega = nom.sample(&mut rng) * DEG2RAD;
        let omega_big = nbom.sample(&mut rng) * DEG2RAD;
        let m = um.sample(&mut rng);

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
        let gal_b = galactic_latitude_deg(&pos);
        let v = planet_apparent_magnitude(sol.mass_earth, sol.albedo, r);

        prob[grid.index(ra, dec)] += 1.0;
        dists.push(r);
        vmags.push(v);

        for (k, t) in telescopes.iter().enumerate() {
            let in_fp = t.footprint.accepts(dec, gal_b);
            let bright = v <= t.limiting_mag;
            if in_fp {
                counts[k].0 += 1;
            }
            if bright {
                counts[k].1 += 1;
            }
            if in_fp && bright {
                counts[k].2 += 1;
            }
        }
    }

    // Normalize the probability grid.
    let total: f64 = prob.iter().sum();
    if total > 0.0 {
        for p in prob.iter_mut() {
            *p /= total;
        }
    }

    // Peak cell center.
    let (peak_idx, _) = prob
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .unwrap();
    let peak_iy = peak_idx / grid.n_ra;
    let peak_ix = peak_idx % grid.n_ra;
    let peak_ra_deg = grid.ra_min_deg + (peak_ix as f64 + 0.5) * grid.dra();
    let peak_dec_deg = grid.dec_center(peak_iy);

    // Smallest cos(Dec)-weighted area enclosing 68% / 95%.
    let mut cells: Vec<(f64, f64)> = prob
        .iter()
        .enumerate()
        .map(|(idx, &p)| (p, grid.dec_center(idx / grid.n_ra)))
        .collect();
    cells.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
    let cell_base = grid.dra() * grid.ddec();
    let (mut area68, mut area95, mut acc) = (0.0, 0.0, 0.0);
    for (p, dec_c) in &cells {
        let solid = cell_base * (dec_c * DEG2RAD).cos();
        if acc < 0.68 {
            area68 += solid;
        }
        if acc < 0.95 {
            area95 += solid;
        }
        acc += p;
        if acc >= 0.95 {
            break;
        }
    }

    dists.sort_by(|a, b| a.partial_cmp(b).unwrap());
    vmags.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let n = n_samples as f64;
    let detections = telescopes
        .iter()
        .zip(counts.iter())
        .map(|(t, &(in_fp, bright, both))| {
            let cov = t.footprint.coverage_fraction;
            StudyDetection {
                study_name: sol.name.clone(),
                detection_prob: cov * both as f64 / n,
                prob_in_footprint: cov * in_fp as f64 / n,
                prob_bright_enough: bright as f64 / n,
            }
        })
        .collect();

    let result = StudyResult {
        solution: sol.clone(),
        prob,
        peak_ra_deg,
        peak_dec_deg,
        area68_deg2: area68,
        area95_deg2: area95,
        dist_p16_au: percentile(&dists, 0.16),
        dist_median_au: percentile(&dists, 0.50),
        dist_p84_au: percentile(&dists, 0.84),
        v_p16: percentile(&vmags, 0.16),
        v_median: percentile(&vmags, 0.50),
        v_p84: percentile(&vmags, 0.84),
    };
    (result, detections)
}

/// Precompute sky overlays (galactic plane and ecliptic) as RA/Dec polylines,
/// using p9-core's SPICE-derived frame matrices so the plotter does no
/// astronomy of its own.
pub fn overlays(n_points: usize) -> Overlays {
    let g2e = galactic_to_ecliptic_matrix();
    let mut galactic_plane = Vec::with_capacity(n_points);
    let mut ecliptic = Vec::with_capacity(n_points);
    for k in 0..n_points {
        let lon = std::f64::consts::TAU * k as f64 / n_points as f64;
        // Galactic equator (b = 0).
        let gal = Vector3::new(lon.cos(), lon.sin(), 0.0);
        let ecl = g2e * gal;
        let (ra, dec) = ecliptic_vec_to_equatorial_deg(&ecl);
        galactic_plane.push([ra, dec]);
        // Ecliptic (β = 0).
        let ecl0 = Vector3::new(lon.cos(), lon.sin(), 0.0);
        let (ra_e, dec_e) = ecliptic_vec_to_equatorial_deg(&ecl0);
        ecliptic.push([ra_e, dec_e]);
    }
    Overlays {
        galactic_plane,
        ecliptic,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{studies, telescope};

    fn test_grid() -> SkyGrid {
        SkyGrid {
            ra_min_deg: 0.0,
            ra_max_deg: 360.0,
            n_ra: 120,
            dec_min_deg: -60.0,
            dec_max_deg: 60.0,
            n_dec: 40,
        }
    }

    #[test]
    fn probability_grid_normalized() {
        let grid = test_grid();
        let sol = &studies::catalog()[3]; // 2021 MCMC
        let (res, _) = run_study(sol, &grid, &telescope::catalog(), 40_000, 7);
        let s: f64 = res.prob.iter().sum();
        assert!((s - 1.0).abs() < 1e-9, "grid sum = {s}");
        assert_eq!(res.prob.len(), grid.len());
    }

    #[test]
    fn mcmc_2021_peak_near_aphelion_region() {
        let grid = test_grid();
        let sol = &studies::catalog()[3]; // 2021 MCMC
        let (res, _) = run_study(sol, &grid, &telescope::catalog(), 80_000, 2021);
        // Anti-aligned aphelion lobe: RA ~ 4.5–7h (65–105°), Dec ~ +5..+30.
        assert!(
            (55.0..110.0).contains(&res.peak_ra_deg),
            "peak RA = {}",
            res.peak_ra_deg
        );
        assert!(
            (0.0..35.0).contains(&res.peak_dec_deg),
            "peak Dec = {}",
            res.peak_dec_deg
        );
        // Brightness consistent with the hand derivation (V ~ 19–22).
        assert!(
            res.v_median > 18.5 && res.v_median < 22.5,
            "V_median = {}",
            res.v_median
        );
        // Most-probable distance near aphelion (a(1+e) ≈ 494 AU).
        assert!(
            res.dist_median_au > 380.0 && res.dist_median_au < 560.0,
            "d_median = {}",
            res.dist_median_au
        );
    }

    #[test]
    fn newer_solutions_are_brighter_and_tighter() {
        let grid = test_grid();
        let cat = studies::catalog();
        let tel = telescope::catalog();
        let (b2016, _) = run_study(&cat[0], &grid, &tel, 60_000, 11);
        let (b2021, _) = run_study(&cat[3], &grid, &tel, 60_000, 11);
        // 2021 (380 AU) is closer/brighter than 2016 (700 AU).
        assert!(
            b2021.v_median < b2016.v_median,
            "2021 V {} should be < 2016 V {}",
            b2021.v_median,
            b2016.v_median
        );
        assert!(b2021.dist_median_au < b2016.dist_median_au);
        // And its 95% region is smaller (tighter constraints).
        assert!(b2021.area95_deg2 < b2016.area95_deg2);
    }

    #[test]
    fn space_survey_beats_ground_for_detection() {
        let grid = test_grid();
        let tel = telescope::catalog();
        let sol = &studies::catalog()[3];
        let (_, det) = run_study(sol, &grid, &tel, 80_000, 99);
        let rubin = det
            .iter()
            .zip(&tel)
            .find(|(_, t)| t.name.contains("Rubin"))
            .unwrap()
            .0;
        let space = det
            .iter()
            .zip(&tel)
            .find(|(_, t)| t.name.contains("Proposed"))
            .unwrap()
            .0;
        // All detection probabilities are valid.
        for d in &det {
            assert!((0.0..=1.0).contains(&d.detection_prob));
            assert!(d.detection_prob <= d.prob_in_footprint + 1e-12);
        }
        // The deep all-sky space survey detects more of the posterior than the
        // plane-avoiding southern ground survey.
        assert!(
            space.detection_prob > rubin.detection_prob,
            "space {} vs rubin {}",
            space.detection_prob,
            rubin.detection_prob
        );
    }

    #[test]
    fn overlays_are_closed_curves() {
        let o = overlays(360);
        assert_eq!(o.galactic_plane.len(), 360);
        assert_eq!(o.ecliptic.len(), 360);
        for p in o.galactic_plane.iter().chain(o.ecliptic.iter()) {
            assert!((0.0..=360.0).contains(&p[0]));
            assert!((-90.0..=90.0).contains(&p[1]));
        }
    }
}

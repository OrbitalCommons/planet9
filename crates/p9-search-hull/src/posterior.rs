//! The "where could it be" side: an ensemble current-position sky posterior
//! built from every orbit solution in [`p9_survey::studies`].
//!
//! Each study is sampled exactly as `p9_survey::skymap` does — perturb the
//! elements by 1σ, draw the mean anomaly uniformly (equal probability per unit
//! time → a dwell-weighted "where is it now" prior), propagate with p9-core,
//! and read off the equatorial direction, heliocentric distance, galactic
//! latitude, and reflected-light V. We give every study equal weight, so the
//! result is a study-agnostic ensemble. Beyond the probability grid we keep, per
//! cell, the probability-weighted mean distance and mean V — the inputs the
//! coverage hull needs to decide whether that direction has been searched
//! deeply enough.

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

/// Galactic latitude (deg) of a heliocentric *ecliptic* position vector.
pub fn galactic_latitude_deg(ecl_pos: &Vector3<f64>) -> f64 {
    let eq = ecliptic_to_equatorial_matrix() * ecl_pos;
    let gal = equatorial_to_galactic_matrix() * eq;
    (gal.z / gal.norm()).asin().to_degrees()
}

/// One Monte Carlo draw's observable position and brightness.
#[derive(Clone, Copy)]
pub struct Sample {
    pub dist_au: f64,
    pub v_mag: f64,
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
    /// A thinned cloud of (distance, V) draws for the parameter-space hull.
    pub cloud: Vec<Sample>,
}

/// Build the equal-weight ensemble posterior. `per_study` MC draws are taken
/// from each solution; `cloud_keep` thins the returned (distance, V) cloud.
pub fn ensemble(
    studies: &[OrbitSolution],
    grid: &SkyGrid,
    per_study: usize,
    seed: u64,
    cloud_keep: usize,
) -> SkyPosterior {
    let n = grid.len();
    let mut count = vec![0.0_f64; n];
    let mut sum_dist = vec![0.0_f64; n];
    let mut sum_v = vec![0.0_f64; n];
    let mut cloud = Vec::new();

    let total_samples = studies.len().max(1) * per_study;
    let cloud_stride = (total_samples / cloud_keep.max(1)).max(1);
    let mut drawn = 0usize;

    for (s_idx, sol) in studies.iter().enumerate() {
        // Per-study seed so the ensemble is reproducible and each study is
        // independent.
        let mut rng = StdRng::seed_from_u64(seed.wrapping_add(s_idx as u64));
        let na = Normal::new(sol.a_au, sol.a_sigma_au).unwrap();
        let ne = Normal::new(sol.e, sol.e_sigma).unwrap();
        let ni = Normal::new(sol.i_deg, sol.i_sigma_deg).unwrap();
        let nom = Normal::new(sol.omega_deg, sol.omega_sigma_deg).unwrap();
        let nbom = Normal::new(sol.omega_big_deg, sol.omega_big_sigma_deg).unwrap();
        let um = Uniform::new(0.0, std::f64::consts::TAU);

        for _ in 0..per_study {
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
            let v = planet_apparent_magnitude(sol.mass_earth, sol.albedo, r);

            let idx = grid.index(ra, dec);
            count[idx] += 1.0;
            sum_dist[idx] += r;
            sum_v[idx] += v;

            if drawn % cloud_stride == 0 {
                cloud.push(Sample {
                    dist_au: r,
                    v_mag: v,
                });
            }
            drawn += 1;
        }
    }

    let total: f64 = count.iter().sum();
    let mut prob = vec![0.0_f64; n];
    let mut dist_mean = vec![f64::NAN; n];
    let mut v_mean = vec![f64::NAN; n];
    for k in 0..n {
        if count[k] > 0.0 {
            prob[k] = count[k] / total;
            dist_mean[k] = sum_dist[k] / count[k];
            v_mean[k] = sum_v[k] / count[k];
        }
    }

    SkyPosterior {
        prob,
        dist_mean,
        v_mean,
        cloud,
    }
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
    fn posterior_normalizes_and_fills_cells() {
        let studies = p9_survey::studies::catalog();
        let grid = small_grid();
        let post = ensemble(&studies, &grid, 2000, 7, 500);
        let total: f64 = post.prob.iter().sum();
        assert!((total - 1.0).abs() < 1e-9, "prob sums to {total}");
        // At least some cells are populated with finite predicted distance / V.
        let filled = post.dist_mean.iter().filter(|d| d.is_finite()).count();
        assert!(filled > 0);
        // Predicted distances are in the outer-solar-system range.
        for (k, d) in post.dist_mean.iter().enumerate() {
            if d.is_finite() {
                assert!(*d > 60.0 && *d < 2000.0, "cell {k}: dist {d}");
                assert!(post.v_mean[k] > 15.0 && post.v_mean[k] < 30.0);
            }
        }
        assert!(!post.cloud.is_empty());
    }
}

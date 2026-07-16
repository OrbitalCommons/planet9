//! Model 2 of Brown & Batygin (2021), Section 3: kernel density estimation
//! of the conditional distributions P_{j,k}[(i, delta_varpi, delta_Omega) |
//! (a, e) bin].
//!
//! For each simulation grid point j the pooled post-burn particle samples
//! (Model 1) are binned by (a, e); within the bin matching an observed
//! ETNO's (a, e), a 3-D kernel density estimate over (i, delta_varpi,
//! delta_Omega) gives the likelihood of that object's observed angles.
//!
//! Bandwidth law (Brown & Batygin 2021, Section 3): the distributions are
//! smoothed with a fractional bandwidth of 5% below a = 230 AU, rising
//! linearly to 30% at a = 730 AU. The paper does not spell out how the
//! fraction maps onto each coordinate; here it multiplies the natural domain
//! width of each dimension (documented interpretation):
//!
//! - delta_varpi, delta_Omega: von Mises kernels with sigma = f(a) * 2 pi,
//!   concentration kappa = 1/sigma^2 (the small-angle Gaussian limit);
//! - inclination: Gaussian kernel with sigma = f(a) * pi/2, reflected at
//!   i = 0 so density does not leak to negative inclinations.
//!
//! Empty/sparse bins fall back to the full sample set, and a log-density
//! floor keeps the likelihood finite (both documented below) — at reduced
//! simulation scale sparse bins are common.

use std::f64::consts::PI;

use serde::{Deserialize, Serialize};

use p9_core::constants::{DEG2RAD, TWO_PI};
use p9_core::data::etno::Etno;

use crate::sim_grid::ElementSample;
use p9_core::analysis::circular::bessel_i0;

/// Below this semi-major axis the fractional bandwidth is 5%
/// (Brown & Batygin 2021, Section 3).
pub const BANDWIDTH_A_LOW: f64 = 230.0;
/// At this semi-major axis the fractional bandwidth has risen to 30%.
pub const BANDWIDTH_A_HIGH: f64 = 730.0;
/// Fractional bandwidth below `BANDWIDTH_A_LOW`.
pub const BANDWIDTH_FRAC_LOW: f64 = 0.05;
/// Fractional bandwidth at and above `BANDWIDTH_A_HIGH`.
pub const BANDWIDTH_FRAC_HIGH: f64 = 0.30;

/// Floor on the per-object log-density. Keeps the total log-likelihood
/// finite when an observation falls in a region with no simulation support
/// (exp(-30) ~ 1e-13, far below any resolved density).
pub const LOG_DENSITY_FLOOR: f64 = -30.0;

/// Fewer samples than this in the matched (a, e) bin triggers the
/// full-sample fallback (documented reduced-scale necessity).
pub const MIN_BIN_SAMPLES: usize = 10;

/// The paper's bandwidth law: 5% fractional bandwidth for a < 230 AU,
/// linear in a up to 30% at a = 730 AU, constant beyond.
pub fn fractional_bandwidth(a: f64) -> f64 {
    let t = ((a - BANDWIDTH_A_LOW) / (BANDWIDTH_A_HIGH - BANDWIDTH_A_LOW)).clamp(0.0, 1.0);
    BANDWIDTH_FRAC_LOW + t * (BANDWIDTH_FRAC_HIGH - BANDWIDTH_FRAC_LOW)
}

/// (a, e) bin edges for the conditional distributions.
///
/// The defaults bracket both the initial disk (a in 150-500 AU, e in
/// ~0.67-0.94) and the observed ETNO sample (a in 228-506 AU, e in
/// 0.69-0.93). Out-of-range values are clamped into the nearest bin
/// (documented: at reduced scale discarding diffused particles starves the
/// estimate).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinSpec {
    pub a_edges: Vec<f64>,
    pub e_edges: Vec<f64>,
}

impl Default for BinSpec {
    fn default() -> Self {
        Self {
            a_edges: vec![150.0, 250.0, 350.0, 450.0, 600.0],
            e_edges: vec![0.5, 0.7, 0.8, 0.9, 1.0],
        }
    }
}

impl BinSpec {
    /// Bin index for (a, e), clamping out-of-range values into the nearest
    /// edge bin.
    pub fn bin_index(&self, a: f64, e: f64) -> (usize, usize) {
        (clamped_bin(&self.a_edges, a), clamped_bin(&self.e_edges, e))
    }

    fn n_bins(&self) -> (usize, usize) {
        (self.a_edges.len() - 1, self.e_edges.len() - 1)
    }
}

fn clamped_bin(edges: &[f64], x: f64) -> usize {
    let n = edges.len() - 1;
    if x < edges[0] {
        return 0;
    }
    for k in 0..n {
        if x < edges[k + 1] {
            return k;
        }
    }
    n - 1
}

/// Von Mises density on [0, 2 pi).
fn von_mises_pdf(theta: f64, mu: f64, kappa: f64) -> f64 {
    (kappa * (theta - mu).cos()).exp() / (TWO_PI * bessel_i0(kappa))
}

/// Gaussian density reflected at i = 0 (inclination support is i >= 0).
fn reflected_gaussian_pdf(x: f64, mu: f64, sigma: f64) -> f64 {
    let norm = 1.0 / (sigma * (TWO_PI).sqrt());
    let g = |m: f64| norm * (-0.5 * ((x - m) / sigma).powi(2)).exp();
    g(mu) + g(-mu)
}

/// 3-D conditional KDE over (i, delta_varpi, delta_Omega) for one (a, e)
/// bin, with kernel widths set by the paper's bandwidth law evaluated at
/// the *observed* object's semi-major axis.
#[derive(Debug, Clone)]
pub struct ConditionalKde {
    samples: Vec<(f64, f64, f64)>,
    /// Gaussian sigma for inclination (radians)
    sigma_i: f64,
    /// Von Mises concentration for the two angle dimensions
    kappa_angle: f64,
}

impl ConditionalKde {
    /// Build a KDE from bin samples, with bandwidths from
    /// `fractional_bandwidth(a_obs)`. Returns None for an empty sample set.
    pub fn new(samples: Vec<(f64, f64, f64)>, a_obs: f64) -> Option<Self> {
        if samples.is_empty() {
            return None;
        }
        let f = fractional_bandwidth(a_obs);
        let sigma_angle = f * TWO_PI;
        let sigma_i = f * (PI / 2.0);
        Some(Self {
            samples,
            sigma_i,
            kappa_angle: 1.0 / (sigma_angle * sigma_angle),
        })
    }

    /// Density at (i, delta_varpi, delta_Omega), angles in radians.
    pub fn density(&self, i: f64, delta_varpi: f64, delta_omega_big: f64) -> f64 {
        let mut total = 0.0;
        for &(si, svarpi, somega) in &self.samples {
            total += reflected_gaussian_pdf(i, si, self.sigma_i)
                * von_mises_pdf(delta_varpi, svarpi, self.kappa_angle)
                * von_mises_pdf(delta_omega_big, somega, self.kappa_angle);
        }
        total / self.samples.len() as f64
    }
}

/// All conditional KDE inputs for one simulation grid point: samples binned
/// by (a, e), with the documented sparse-bin fallback.
#[derive(Debug, Clone)]
pub struct GridPointKde {
    spec: BinSpec,
    /// bins\[ia\]\[ie\] -> (i, delta_varpi, delta_Omega) tuples
    bins: Vec<Vec<Vec<(f64, f64, f64)>>>,
    all: Vec<(f64, f64, f64)>,
}

impl GridPointKde {
    pub fn from_samples(samples: &[ElementSample], spec: BinSpec) -> Self {
        let (na, ne) = spec.n_bins();
        let mut bins = vec![vec![Vec::new(); ne]; na];
        let mut all = Vec::with_capacity(samples.len());
        for s in samples {
            let tuple = (s.i, s.delta_varpi, s.delta_omega_big);
            let (ia, ie) = spec.bin_index(s.a, s.e);
            bins[ia][ie].push(tuple);
            all.push(tuple);
        }
        Self { spec, bins, all }
    }

    /// Number of samples in the bin matched by (a, e).
    pub fn bin_count(&self, a: f64, e: f64) -> usize {
        let (ia, ie) = self.spec.bin_index(a, e);
        self.bins[ia][ie].len()
    }

    /// Log-density of observed (i, delta_varpi, delta_Omega) given the
    /// observation's (a, e) bin. Falls back to the full sample set when the
    /// bin holds fewer than `MIN_BIN_SAMPLES` samples; floored at
    /// `LOG_DENSITY_FLOOR`.
    pub fn log_density(
        &self,
        a_obs: f64,
        e_obs: f64,
        i: f64,
        delta_varpi: f64,
        delta_omega_big: f64,
    ) -> f64 {
        let (ia, ie) = self.spec.bin_index(a_obs, e_obs);
        let bin = &self.bins[ia][ie];
        let samples = if bin.len() >= MIN_BIN_SAMPLES {
            bin
        } else {
            &self.all
        };
        match ConditionalKde::new(samples.clone(), a_obs) {
            Some(kde) => kde
                .density(i, delta_varpi, delta_omega_big)
                .ln()
                .max(LOG_DENSITY_FLOOR),
            None => LOG_DENSITY_FLOOR,
        }
    }
}

/// Total log-likelihood of the observed ETNO sample given this grid point's
/// conditional distributions and a hypothesized Planet Nine orientation
/// (varpi9, Omega9) (radians).
///
/// The simulations record particle angles *relative* to Planet Nine, so the
/// observed relative angles are delta_varpi = varpi_etno - varpi9 and
/// delta_Omega = Omega_etno - Omega9.
pub fn etno_log_likelihood(
    kde: &GridPointKde,
    etnos: &[Etno],
    varpi9: f64,
    omega_big9: f64,
) -> f64 {
    etnos
        .iter()
        .map(|etno| {
            let delta_varpi = (etno.longitude_of_perihelion() - varpi9).rem_euclid(TWO_PI);
            let delta_omega_big = (etno.omega_big_deg * DEG2RAD - omega_big9).rem_euclid(TWO_PI);
            kde.log_density(
                etno.a,
                etno.e,
                etno.i_deg * DEG2RAD,
                delta_varpi,
                delta_omega_big,
            )
        })
        .sum()
}

/// Profile the ETNO log-likelihood over the Planet Nine orientation angles.
///
/// The paper's MCMC includes (varpi9, Omega9) as free parameters; the
/// reduced 4-parameter pipeline here (m9, a9, e9, i9) instead *profiles*
/// them out: a coarse `n_scan` x `n_scan` scan over [0, 2 pi)^2 keeps the
/// maximizing orientation (documented reduction; 24 x 24 = 15 deg
/// resolution is finer than any kernel bandwidth in use).
///
/// Returns (max log-likelihood, best varpi9, best Omega9).
pub fn profiled_log_likelihood(
    kde: &GridPointKde,
    etnos: &[Etno],
    n_scan: usize,
) -> (f64, f64, f64) {
    let mut best = (f64::NEG_INFINITY, 0.0, 0.0);
    for jv in 0..n_scan {
        let varpi9 = TWO_PI * jv as f64 / n_scan as f64;
        for jo in 0..n_scan {
            let omega_big9 = TWO_PI * jo as f64 / n_scan as f64;
            let ll = etno_log_likelihood(kde, etnos, varpi9, omega_big9);
            if ll > best.0 {
                best = (ll, varpi9, omega_big9);
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::data::etno::BROWN_2017_SAMPLE;

    #[test]
    fn test_bandwidth_law() {
        // Paper values: 5% below 230 AU, 30% at 730 AU, linear between.
        assert_relative_eq!(fractional_bandwidth(150.0), 0.05);
        assert_relative_eq!(fractional_bandwidth(230.0), 0.05);
        assert_relative_eq!(fractional_bandwidth(480.0), 0.175); // midpoint
        assert_relative_eq!(fractional_bandwidth(730.0), 0.30);
        assert_relative_eq!(fractional_bandwidth(1000.0), 0.30);
    }

    #[test]
    fn test_bessel_i0_reference_values() {
        // Abramowitz & Stegun Table 9.8: I0(0)=1, I0(1)=1.2660658,
        // I0(2)=2.2795853, I0(5)=27.239872.
        assert_relative_eq!(bessel_i0(0.0), 1.0, epsilon = 1e-6);
        assert_relative_eq!(bessel_i0(1.0), 1.2660658, epsilon = 1e-5);
        assert_relative_eq!(bessel_i0(2.0), 2.2795853, epsilon = 1e-5);
        assert_relative_eq!(bessel_i0(5.0), 27.239872, max_relative = 1e-5);
    }

    #[test]
    fn test_von_mises_normalization() {
        // Numerically integrate the von Mises pdf over the circle.
        for &kappa in &[0.3, 2.0, 10.0] {
            let n = 5000;
            let sum: f64 = (0..n)
                .map(|k| von_mises_pdf(TWO_PI * k as f64 / n as f64, 1.0, kappa))
                .sum::<f64>()
                * TWO_PI
                / n as f64;
            assert_relative_eq!(sum, 1.0, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_reflected_gaussian_integrates_to_one_on_half_line() {
        let (mu, sigma) = (0.05, 0.1);
        let n = 20000;
        let dx = 2.0 / n as f64;
        let sum: f64 = (0..n)
            .map(|k| reflected_gaussian_pdf((k as f64 + 0.5) * dx, mu, sigma))
            .sum::<f64>()
            * dx;
        assert_relative_eq!(sum, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_bin_clamping() {
        let spec = BinSpec::default();
        // In-range values land in the right bins.
        assert_eq!(spec.bin_index(160.0, 0.55), (0, 0));
        assert_eq!(spec.bin_index(500.0, 0.95), (3, 3));
        // Out-of-range values clamp to the nearest edge bin.
        assert_eq!(spec.bin_index(100.0, 0.2), (0, 0));
        assert_eq!(spec.bin_index(2000.0, 0.999), (3, 3));
    }

    #[test]
    fn test_kde_peaks_near_samples() {
        // Tight cluster of samples at i=0.3, dvarpi=pi, dOmega=0.
        let samples: Vec<(f64, f64, f64)> = (0..50)
            .map(|k| {
                let eps = 0.01 * (k as f64 - 25.0) / 25.0;
                (0.3 + eps, PI + eps, eps.rem_euclid(TWO_PI))
            })
            .collect();
        let kde = ConditionalKde::new(samples, 300.0).unwrap();
        let at_cluster = kde.density(0.3, PI, 0.0);
        let far_away = kde.density(1.2, 0.0, PI);
        assert!(
            at_cluster > far_away,
            "KDE not peaked at samples: {at_cluster} <= {far_away}"
        );
    }

    #[test]
    fn test_grid_point_kde_fallback_and_floor() {
        use crate::sim_grid::ElementSample;
        // 5 samples in one bin: below MIN_BIN_SAMPLES, so every lookup uses
        // the full-set fallback and still returns a finite value.
        let samples: Vec<ElementSample> = (0..5)
            .map(|k| ElementSample {
                a: 300.0,
                e: 0.85,
                i: 0.3,
                delta_varpi: PI + 0.01 * k as f64,
                delta_omega_big: 0.1,
            })
            .collect();
        let kde = GridPointKde::from_samples(&samples, BinSpec::default());
        assert_eq!(kde.bin_count(300.0, 0.85), 5);
        let ll = kde.log_density(450.0, 0.75, 0.3, PI, 0.1);
        assert!(ll.is_finite() && ll > LOG_DENSITY_FLOOR);

        // Empty sample set -> floor.
        let empty = GridPointKde::from_samples(&[], BinSpec::default());
        assert_eq!(
            empty.log_density(300.0, 0.85, 0.3, PI, 0.1),
            LOG_DENSITY_FLOOR
        );
    }

    #[test]
    fn test_profiled_likelihood_finds_clustered_orientation() {
        use crate::sim_grid::ElementSample;
        // Synthetic anti-aligned population: delta_varpi concentrated at pi,
        // delta_Omega at 0, inclinations near the ETNO range.
        let mut samples = Vec::new();
        for k in 0..400 {
            let jitter = 0.25 * ((k * 37 % 100) as f64 / 100.0 - 0.5);
            samples.push(ElementSample {
                a: 250.0 + (k % 30) as f64 * 8.0,
                e: 0.7 + (k % 25) as f64 * 0.01,
                i: (15.0 + 10.0 * ((k * 13 % 50) as f64 / 50.0 - 0.5)) * DEG2RAD,
                delta_varpi: (PI + jitter).rem_euclid(TWO_PI),
                delta_omega_big: jitter.rem_euclid(TWO_PI),
            });
        }
        let kde = GridPointKde::from_samples(&samples, BinSpec::default());
        let (ll_best, varpi9, _) = profiled_log_likelihood(&kde, &BROWN_2017_SAMPLE, 24);
        assert!(ll_best.is_finite());

        // The clustered ETNO varpi values sit near ~50 deg (circular mean of
        // the 8-object cluster); anti-alignment puts the preferred varpi9
        // roughly opposite. The profile maximum must beat a deliberately
        // orthogonal orientation.
        let ll_off = etno_log_likelihood(&kde, &BROWN_2017_SAMPLE, varpi9 + PI / 2.0, 0.0);
        assert!(ll_best >= ll_off);
    }
}

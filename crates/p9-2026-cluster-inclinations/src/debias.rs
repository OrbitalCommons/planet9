//! Brown (2001) debiasing of the observed inclination distribution
//! (paper Section 3.1, Figure 1).
//!
//! Objects are found preferentially near the ecliptic, so their raw
//! inclinations are biased low. For an object discovered at latitude `β`
//! from the reference plane, the probability that its inclination is at most
//! `i_obs` under an intrinsic distribution `f_t(i) = sin i · exp(−i²/2w²)`
//! is
//!
//! ```text
//!   P = ∫_β^{i_obs} f_t(i) (sin²i − sin²β)^{-1/2} di
//!     / ∫_β^{π/2}   f_t(i) (sin²i − sin²β)^{-1/2} di ,
//! ```
//!
//! the (sin²i − sin²β)^{-1/2} factor being the fraction of time an orbit of
//! inclination `i` spends at latitude `β`. If `w` is right, the `P_j` of the
//! sample are uniform on [0, 1]; the Kuiper statistic `√N·V` of the `P_j`
//! against uniformity, scanned over `w`, is minimised at the intrinsic width.
//! Confidence intervals come from the Monte Carlo calibration of Brown
//! (2001): the statistic's null distribution at each candidate `w` is built
//! from synthetic samples of `N` inclinations drawn from `f_t(·; w)` at the
//! observed latitudes.

use p9_core::analysis::circular::kuiper_statistic;
use p9_core::constants::{DEG2RAD, TWO_PI};
use rand::Rng;
use rayon::prelude::*;

use crate::width::intrinsic_density;

/// Tabulated conditional CDF of inclination given discovery latitude.
#[derive(Debug, Clone)]
pub struct ConditionalCdf {
    /// Inclination grid on [β, π/2] (radians).
    pub i_grid: Vec<f64>,
    /// Cumulative probability at each grid point, 0 at β and 1 at π/2.
    pub cdf: Vec<f64>,
}

impl ConditionalCdf {
    /// Build the table for width `w` and latitude `beta` (radians, β ≥ 0),
    /// with `n` grid intervals. The integrable singularity at i = β is
    /// removed with the substitution i = β + x².
    pub fn new(w: f64, beta: f64, n: usize) -> Self {
        let beta = beta.abs().min(std::f64::consts::FRAC_PI_2 - 1e-9);
        let x_max = (std::f64::consts::FRAC_PI_2 - beta).sqrt();
        let h = x_max / n as f64;
        let sin2b = beta.sin().powi(2);
        let integrand = |x: f64| -> f64 {
            if x <= 0.0 {
                // Limit x → 0: (sin²(β+x²) − sin²β)^{1/2} → x·sqrt(sin 2β).
                return if beta > 1e-12 {
                    2.0 * intrinsic_density(beta, w) / (2.0 * beta).sin().sqrt()
                } else {
                    0.0
                };
            }
            let i = beta + x * x;
            let d = (i.sin().powi(2) - sin2b).max(0.0).sqrt();
            if d <= 0.0 {
                0.0
            } else {
                2.0 * x * intrinsic_density(i, w) / d
            }
        };
        let mut i_grid = Vec::with_capacity(n + 1);
        let mut cdf = Vec::with_capacity(n + 1);
        let mut acc = 0.0;
        let mut prev = integrand(0.0);
        i_grid.push(beta);
        cdf.push(0.0);
        for k in 1..=n {
            let x = k as f64 * h;
            let cur = integrand(x);
            acc += 0.5 * (prev + cur) * h;
            prev = cur;
            i_grid.push(beta + x * x);
            cdf.push(acc);
        }
        let total = acc.max(1e-300);
        for c in cdf.iter_mut() {
            *c /= total;
        }
        Self { i_grid, cdf }
    }

    /// P(i ≤ i_obs | β), linearly interpolated.
    pub fn eval(&self, i_obs: f64) -> f64 {
        interp(&self.i_grid, &self.cdf, i_obs).clamp(0.0, 1.0)
    }

    /// Inverse CDF: the inclination at cumulative probability `u`.
    pub fn quantile(&self, u: f64) -> f64 {
        interp(&self.cdf, &self.i_grid, u.clamp(0.0, 1.0))
    }
}

fn interp(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    if x <= xs[0] {
        return ys[0];
    }
    let last = xs.len() - 1;
    if x >= xs[last] {
        return ys[last];
    }
    let k = xs.partition_point(|&v| v <= x).max(1);
    let (x0, x1) = (xs[k - 1], xs[k]);
    let t = if x1 > x0 { (x - x0) / (x1 - x0) } else { 0.0 };
    ys[k - 1] + t * (ys[k] - ys[k - 1])
}

/// Number of grid intervals in each conditional CDF table.
pub const CDF_GRID: usize = 400;

/// Kuiper `√N·V` of the `P_j` for width `w`, given discovery latitudes and
/// observed inclinations (radians).
pub fn kuiper_sqrt_n(w: f64, betas: &[f64], incs: &[f64]) -> f64 {
    let tables: Vec<ConditionalCdf> = betas
        .iter()
        .map(|&b| ConditionalCdf::new(w, b, CDF_GRID))
        .collect();
    kuiper_sqrt_n_with_tables(&tables, incs)
}

fn kuiper_sqrt_n_with_tables(tables: &[ConditionalCdf], incs: &[f64]) -> f64 {
    let p: Vec<f64> = tables
        .iter()
        .zip(incs)
        .map(|(t, &i)| t.eval(i) * TWO_PI)
        .collect();
    (p.len() as f64).sqrt() * kuiper_statistic(&p)
}

/// Result of the width scan.
#[derive(Debug, Clone)]
pub struct WidthFit {
    /// Candidate widths scanned (degrees).
    pub w_grid_deg: Vec<f64>,
    /// `√N·V` at each candidate width.
    pub statistic: Vec<f64>,
    /// Monte-Carlo p-value at each candidate width: the fraction of synthetic
    /// samples drawn from `f_t(·; w)` whose statistic is at least the observed
    /// one.
    pub p_value: Vec<f64>,
    /// Width minimising the statistic (degrees).
    pub w_best_deg: f64,
}

impl WidthFit {
    /// Smallest and largest scanned widths not rejected at the given
    /// two-sided-equivalent threshold (`0.159` for 1σ, `0.0228` for 2σ,
    /// `0.00135` for 3σ), following Brown (2001).
    pub fn interval(&self, p_threshold: f64) -> (f64, f64) {
        let accepted: Vec<f64> = self
            .w_grid_deg
            .iter()
            .zip(&self.p_value)
            .filter(|(_, &p)| p >= p_threshold)
            .map(|(&w, _)| w)
            .collect();
        match (accepted.first(), accepted.last()) {
            (Some(&lo), Some(&hi)) => (lo, hi),
            _ => (f64::NAN, f64::NAN),
        }
    }

    /// p-value at the scanned width nearest `w_deg`.
    pub fn p_at(&self, w_deg: f64) -> f64 {
        let k = self
            .w_grid_deg
            .iter()
            .enumerate()
            .min_by(|a, b| {
                (a.1 - w_deg)
                    .abs()
                    .partial_cmp(&(b.1 - w_deg).abs())
                    .unwrap()
            })
            .map(|(k, _)| k)
            .unwrap();
        self.p_value[k]
    }
}

/// One-sigma acceptance threshold on the calibrated p-value (84.1% level).
pub const P_ONE_SIGMA: f64 = 1.0 - 0.841;
/// Three-sigma acceptance threshold (99.87% level).
pub const P_THREE_SIGMA: f64 = 1.0 - 0.99865;

/// Scan candidate widths, computing the statistic and its Monte-Carlo
/// calibrated p-value. `betas` and `incs` in radians, `w_grid_deg` in
/// degrees, `n_mc` synthetic samples per width.
pub fn fit_intrinsic_width(
    betas: &[f64],
    incs: &[f64],
    w_grid_deg: &[f64],
    n_mc: usize,
    seed: u64,
) -> WidthFit {
    assert_eq!(betas.len(), incs.len());
    let rows: Vec<(f64, f64)> = w_grid_deg
        .par_iter()
        .enumerate()
        .map(|(k, &w_deg)| {
            let w = w_deg * DEG2RAD;
            let tables: Vec<ConditionalCdf> = betas
                .iter()
                .map(|&b| ConditionalCdf::new(w, b, CDF_GRID))
                .collect();
            let observed = kuiper_sqrt_n_with_tables(&tables, incs);
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed.wrapping_add(k as u64));
            let mut exceed = 0usize;
            let mut synth = vec![0.0; incs.len()];
            for _ in 0..n_mc {
                for (s, t) in synth.iter_mut().zip(&tables) {
                    *s = t.quantile(rng.gen::<f64>());
                }
                if kuiper_sqrt_n_with_tables(&tables, &synth) >= observed {
                    exceed += 1;
                }
            }
            (observed, exceed as f64 / n_mc as f64)
        })
        .collect();
    let statistic: Vec<f64> = rows.iter().map(|r| r.0).collect();
    let p_value: Vec<f64> = rows.iter().map(|r| r.1).collect();
    let best = statistic
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(k, _)| k)
        .unwrap();
    WidthFit {
        w_grid_deg: w_grid_deg.to_vec(),
        statistic,
        p_value,
        w_best_deg: w_grid_deg[best],
    }
}

use rand::SeedableRng;

/// Default scan: 3° to 45° in 1° steps.
pub fn default_w_grid() -> Vec<f64> {
    (3..=45).map(|k| k as f64).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conditional_cdf_is_monotone_and_normalised() {
        let t = ConditionalCdf::new(12.0 * DEG2RAD, 5.0 * DEG2RAD, CDF_GRID);
        assert!((t.cdf[0]).abs() < 1e-12);
        assert!((t.cdf[t.cdf.len() - 1] - 1.0).abs() < 1e-12);
        for k in 1..t.cdf.len() {
            assert!(t.cdf[k] >= t.cdf[k - 1]);
        }
        // Round trip.
        let i = 20.0 * DEG2RAD;
        assert!((t.quantile(t.eval(i)) - i).abs() < 1e-3);
    }

    #[test]
    fn ecliptic_latitude_reduces_to_plain_cdf() {
        // At β = 0 the weight is 1/sin i, so the conditional density is
        // exp(−i²/2w²) — a half-Gaussian truncated at π/2.
        let w = 10.0 * DEG2RAD;
        let t = ConditionalCdf::new(w, 0.0, CDF_GRID);
        let i = 10.0 * DEG2RAD;
        let erf = |x: f64| {
            // Abramowitz–Stegun 7.1.26
            let s = x.signum();
            let x = x.abs();
            let t = 1.0 / (1.0 + 0.3275911 * x);
            let y = 1.0
                - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
                    + 0.254829592)
                    * t
                    * (-x * x).exp();
            s * y
        };
        let expect =
            erf(i / (w * 2f64.sqrt())) / erf(std::f64::consts::FRAC_PI_2 / (w * 2f64.sqrt()));
        assert!(
            (t.eval(i) - expect).abs() < 2e-3,
            "{} vs {}",
            t.eval(i),
            expect
        );
    }

    #[test]
    fn scan_recovers_synthetic_width() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(3);
        let w_true = 15.0 * DEG2RAD;
        let betas: Vec<f64> = (0..60).map(|_| rng.gen_range(0.0..8.0 * DEG2RAD)).collect();
        let incs: Vec<f64> = betas
            .iter()
            .map(|&b| ConditionalCdf::new(w_true, b, CDF_GRID).quantile(rng.gen()))
            .collect();
        let fit = fit_intrinsic_width(&betas, &incs, &default_w_grid(), 200, 1);
        assert!(
            (fit.w_best_deg - 15.0).abs() <= 4.0,
            "best w = {}",
            fit.w_best_deg
        );
        let (lo, hi) = fit.interval(P_ONE_SIGMA);
        assert!(
            lo <= 15.0 && hi >= 15.0,
            "1σ interval [{lo}, {hi}] misses 15°"
        );
    }
}

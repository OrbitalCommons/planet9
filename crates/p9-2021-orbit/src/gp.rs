//! Model 3 of Brown & Batygin (2021), Section 4: Gaussian-process emulation
//! of the simulation-grid likelihood.
//!
//! The paper interpolates the per-grid-point likelihood over continuous
//! (m9, a9, e9, i9) with a GP using a Matern nu = 3/2 kernel, so the MCMC
//! (Model 4) can evaluate the likelihood off-grid. The grid is small
//! (<= 121 points), so exact GP regression via Cholesky factorization
//! (nalgebra) is used — no sparse/inducing-point machinery.
//!
//! Hyperparameters (shared length scale over the normalized inputs, noise
//! variance) are fit by grid-search maximization of the log marginal
//! likelihood — documented choice: with <= 121 training points and a
//! 2-parameter search space, an optimizer buys nothing over an exhaustive
//! scan. Targets are standardized internally (zero mean, unit variance), so
//! the signal variance is fixed at 1.
//!
//! Leave-one-out cross-validation uses the closed-form identities
//! mu_i = y_i - [K^-1 y]_i / [K^-1]_ii (Rasmussen & Williams 2006, eq. 5.12).

use nalgebra::{Cholesky, DMatrix, DVector, Dyn};
use serde::{Deserialize, Serialize};

/// Matern nu = 3/2 correlation as a function of the scaled distance r:
/// k(r) = (1 + sqrt(3) r) exp(-sqrt(3) r).
pub fn matern32(r: f64) -> f64 {
    let s = 3.0_f64.sqrt() * r;
    (1.0 + s) * (-s).exp()
}

/// Hyperparameter search grids (normalized-input units / standardized-target
/// units). Length scales span "wiggly within the box" to "nearly linear";
/// noise spans interpolation to heavy regularization.
const LENGTH_SCALE_GRID: [f64; 10] = [0.1, 0.15, 0.25, 0.4, 0.6, 0.9, 1.3, 2.0, 3.0, 5.0];
const NOISE_VAR_GRID: [f64; 6] = [1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 0.5];

/// Exact GP regressor with a Matern 3/2 kernel over normalized inputs.
#[derive(Debug, Clone)]
pub struct Matern32Gp {
    x: Vec<[f64; 4]>,
    /// Standardized targets
    y: DVector<f64>,
    y_mean: f64,
    y_std: f64,
    /// Fitted hyperparameters
    pub length_scale: f64,
    pub noise_var: f64,
    /// Log marginal likelihood at the fitted hyperparameters
    pub log_marginal: f64,
    chol: Cholesky<f64, Dyn>,
    alpha: DVector<f64>,
}

/// Summary of the fitted GP, serializable for output records.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GpFitSummary {
    pub n_points: usize,
    pub length_scale: f64,
    pub noise_var: f64,
    pub log_marginal: f64,
    pub loo_relative_error: f64,
}

fn kernel_matrix(x: &[[f64; 4]], length_scale: f64, noise_var: f64) -> DMatrix<f64> {
    let n = x.len();
    DMatrix::from_fn(n, n, |i, j| {
        let r = scaled_distance(&x[i], &x[j], length_scale);
        matern32(r) + if i == j { noise_var } else { 0.0 }
    })
}

fn scaled_distance(a: &[f64; 4], b: &[f64; 4], length_scale: f64) -> f64 {
    let mut s = 0.0;
    for d in 0..4 {
        let dx = (a[d] - b[d]) / length_scale;
        s += dx * dx;
    }
    s.sqrt()
}

impl Matern32Gp {
    /// Fit the GP to (normalized input, target) pairs, grid-searching the
    /// length scale and noise variance by maximum log marginal likelihood.
    ///
    /// Panics if fewer than 3 points are supplied or the targets are
    /// degenerate (zero variance).
    pub fn fit(inputs: &[[f64; 4]], targets: &[f64]) -> Self {
        assert!(inputs.len() == targets.len());
        assert!(inputs.len() >= 3, "GP needs at least 3 training points");
        let n = inputs.len();

        let y_mean = targets.iter().sum::<f64>() / n as f64;
        let var = targets.iter().map(|y| (y - y_mean).powi(2)).sum::<f64>() / n as f64;
        assert!(var > 0.0, "degenerate (constant) GP targets");
        let y_std = var.sqrt();
        let y = DVector::from_iterator(n, targets.iter().map(|t| (t - y_mean) / y_std));

        struct Candidate {
            lml: f64,
            ell: f64,
            noise: f64,
            chol: Cholesky<f64, Dyn>,
            alpha: DVector<f64>,
        }
        let mut best: Option<Candidate> = None;
        for &ell in LENGTH_SCALE_GRID.iter() {
            for &noise in NOISE_VAR_GRID.iter() {
                let k = kernel_matrix(inputs, ell, noise);
                let Some(chol) = Cholesky::new(k) else {
                    continue; // numerically non-PD; skip this combination
                };
                let alpha = chol.solve(&y);
                let log_det: f64 = chol.l().diagonal().iter().map(|d| 2.0 * d.ln()).sum();
                let lml = -0.5 * y.dot(&alpha)
                    - 0.5 * log_det
                    - 0.5 * n as f64 * (2.0 * std::f64::consts::PI).ln();
                if best.as_ref().is_none_or(|b| lml > b.lml) {
                    best = Some(Candidate {
                        lml,
                        ell,
                        noise,
                        chol,
                        alpha,
                    });
                }
            }
        }
        let best = best.expect("no positive-definite kernel matrix in the search grid");

        Self {
            x: inputs.to_vec(),
            y,
            y_mean,
            y_std,
            length_scale: best.ell,
            noise_var: best.noise,
            log_marginal: best.lml,
            chol: best.chol,
            alpha: best.alpha,
        }
    }

    /// Predictive mean and variance at a normalized input point
    /// (de-standardized to target units).
    pub fn predict(&self, x: &[f64; 4]) -> (f64, f64) {
        let n = self.x.len();
        let k_star = DVector::from_iterator(
            n,
            self.x
                .iter()
                .map(|xi| matern32(scaled_distance(xi, x, self.length_scale))),
        );
        let mean_std = k_star.dot(&self.alpha);
        let v = self.chol.solve(&k_star);
        let var_std = (1.0 + self.noise_var - k_star.dot(&v)).max(0.0);
        (
            self.y_mean + self.y_std * mean_std,
            self.y_std * self.y_std * var_std,
        )
    }

    /// Leave-one-out residuals (target units), via the closed-form
    /// K^-1 identities — no n refits.
    pub fn loo_residuals(&self) -> Vec<f64> {
        let n = self.x.len();
        let k_inv = self.chol.inverse();
        let k_inv_y = &k_inv * &self.y;
        (0..n)
            .map(|i| self.y_std * k_inv_y[i] / k_inv[(i, i)])
            .collect()
    }

    /// LOO relative error: mean |residual| divided by the target range.
    ///
    /// Normalizing by the range (not |y_i|) keeps the metric meaningful for
    /// log-likelihood targets, whose absolute values are offset-arbitrary.
    pub fn loo_relative_error(&self) -> f64 {
        let res = self.loo_residuals();
        let n = res.len() as f64;
        let mean_abs = res.iter().map(|r| r.abs()).sum::<f64>() / n;
        let y_min = self.y.iter().cloned().fold(f64::INFINITY, f64::min);
        let y_max = self.y.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let range = (y_max - y_min) * self.y_std;
        mean_abs / range
    }

    pub fn summary(&self) -> GpFitSummary {
        GpFitSummary {
            n_points: self.x.len(),
            length_scale: self.length_scale,
            noise_var: self.noise_var,
            log_marginal: self.log_marginal,
            loo_relative_error: self.loo_relative_error(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Smooth test surface over the unit box.
    fn smooth_fn(x: &[f64; 4]) -> f64 {
        3.0 * x[0] - 2.0 * (x[1] * std::f64::consts::PI).sin() + x[2] * x[3] - 0.5 * x[2]
    }

    /// Deterministic low-discrepancy-ish points in the unit box.
    fn test_inputs(n: usize) -> Vec<[f64; 4]> {
        (0..n)
            .map(|k| {
                let t = k as f64;
                [
                    (t * 0.618_033_988_75).fract(),
                    (t * 0.414_213_562_37).fract(),
                    (t * 0.732_050_807_57).fract(),
                    (t * 0.236_067_977_50).fract(),
                ]
            })
            .collect()
    }

    #[test]
    fn test_matern32_limits() {
        assert_relative_eq!(matern32(0.0), 1.0);
        assert!(matern32(0.5) > matern32(1.0));
        assert!(matern32(10.0) < 1e-6);
    }

    #[test]
    fn test_gp_interpolates_training_points() {
        let x = test_inputs(40);
        let y: Vec<f64> = x.iter().map(smooth_fn).collect();
        let gp = Matern32Gp::fit(&x, &y);
        for (xi, yi) in x.iter().zip(&y) {
            let (mean, _) = gp.predict(xi);
            assert!(
                (mean - yi).abs() < 0.15,
                "training point miss: {mean} vs {yi}"
            );
        }
    }

    #[test]
    fn test_gp_predicts_held_out_points() {
        let x = test_inputs(60);
        let y: Vec<f64> = x.iter().map(smooth_fn).collect();
        let gp = Matern32Gp::fit(&x, &y);
        let y_range = 6.0; // approximate range of smooth_fn on the box
        for probe in test_inputs(97).iter().skip(60) {
            let (mean, var) = gp.predict(probe);
            let truth = smooth_fn(probe);
            assert!(
                (mean - truth).abs() < 0.2 * y_range,
                "off-grid miss: {mean} vs {truth}"
            );
            assert!(var >= 0.0);
        }
    }

    #[test]
    fn test_loo_residuals_match_explicit_refits() {
        // The closed-form LOO must agree with literally refitting without
        // each point (same fixed hyperparameters).
        let x = test_inputs(20);
        let y: Vec<f64> = x.iter().map(smooth_fn).collect();
        let gp = Matern32Gp::fit(&x, &y);
        let res = gp.loo_residuals();

        for hold in [0usize, 7, 19] {
            let xs: Vec<[f64; 4]> = x
                .iter()
                .enumerate()
                .filter(|(k, _)| *k != hold)
                .map(|(_, v)| *v)
                .collect();
            // Standardize with the *full* stats and reuse the fitted
            // hyperparameters, mirroring the closed-form identity.
            let y_mean = y.iter().sum::<f64>() / y.len() as f64;
            let y_std =
                (y.iter().map(|v| (v - y_mean).powi(2)).sum::<f64>() / y.len() as f64).sqrt();
            let ys = DVector::from_iterator(
                xs.len(),
                y.iter()
                    .enumerate()
                    .filter(|(k, _)| *k != hold)
                    .map(|(_, v)| (v - y_mean) / y_std),
            );
            let k = kernel_matrix(&xs, gp.length_scale, gp.noise_var);
            let chol = Cholesky::new(k).unwrap();
            let alpha = chol.solve(&ys);
            let k_star = DVector::from_iterator(
                xs.len(),
                xs.iter()
                    .map(|xi| matern32(scaled_distance(xi, &x[hold], gp.length_scale))),
            );
            let mu = y_mean + y_std * k_star.dot(&alpha);
            assert_relative_eq!(y[hold] - mu, res[hold], epsilon = 1e-8);
        }
    }

    #[test]
    fn test_loo_relative_error_small_on_smooth_function() {
        // 12 points (the Reduced grid size) on the smooth surface: measured
        // LOO relative error 0.203 — 12 points in a 4-D box is genuinely
        // sparse. (The paper quotes < 5% for its 121-point grid; see
        // REPRODUCTION_NOTES.md for the scale difference.)
        let x = test_inputs(12);
        let y: Vec<f64> = x.iter().map(smooth_fn).collect();
        let gp = Matern32Gp::fit(&x, &y);
        let err = gp.loo_relative_error();
        assert!(err < 0.25, "12-point LOO relative error = {err:.4}");

        // At 60 points the same surface must emulate below the paper's 5%.
        let x = test_inputs(60);
        let y: Vec<f64> = x.iter().map(smooth_fn).collect();
        let gp = Matern32Gp::fit(&x, &y);
        let err = gp.loo_relative_error();
        assert!(err < 0.05, "60-point LOO relative error = {err:.4}");
    }
}

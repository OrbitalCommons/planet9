//! Model 4 of Brown & Batygin (2021), Section 5: affine-invariant ensemble
//! MCMC (Goodman & Weare 2010; the emcee "stretch move" with the red-black
//! parallel update of Foreman-Mackey et al. 2013).
//!
//! The sampler is generic over a log-posterior closure. The paper ran 100
//! chains of 20,890 samples, burn-in 260, thinning 42, yielding 49,100
//! uncorrelated posterior samples; those numbers are pinned in
//! [`PaperMcmcSettings`] for the `#[ignore]`d paper-scale entry point, while
//! reduced runs use shorter chains (documented at the call sites).
//!
//! Convergence diagnostics:
//! - Gelman-Rubin R-hat treating each walker's post-burn chain as a chain.
//!   Walkers in an ensemble are not independent chains, so this is a health
//!   metric (R-hat near 1 is necessary, not sufficient) — documented caveat.
//! - Effective sample size from the integrated autocorrelation time
//!   (Sokal windowing, c = 5, autocorrelations averaged across walkers).

use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

/// The paper's MCMC run settings (Brown & Batygin 2021, Section 5).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PaperMcmcSettings {
    pub n_walkers: usize,
    pub n_steps: usize,
    pub burn_in: usize,
    pub thin: usize,
}

impl PaperMcmcSettings {
    pub fn paper() -> Self {
        Self {
            n_walkers: 100,
            n_steps: 20_890,
            burn_in: 260,
            thin: 42,
        }
    }
}

/// Ensemble chain: `samples[step][walker]` is a parameter vector.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnsembleChain {
    pub samples: Vec<Vec<Vec<f64>>>,
    pub log_probs: Vec<Vec<f64>>,
    pub acceptance_rate: f64,
    pub seed: u64,
}

impl EnsembleChain {
    pub fn n_steps(&self) -> usize {
        self.samples.len()
    }

    pub fn n_walkers(&self) -> usize {
        self.samples.first().map_or(0, |s| s.len())
    }

    pub fn dim(&self) -> usize {
        self.samples
            .first()
            .and_then(|s| s.first())
            .map_or(0, |w| w.len())
    }

    /// Flattened post-burn, thinned samples (all walkers interleaved).
    pub fn flat_samples(&self, burn_in: usize, thin: usize) -> Vec<Vec<f64>> {
        let mut out = Vec::new();
        for step in (burn_in..self.n_steps()).step_by(thin.max(1)) {
            for walker in &self.samples[step] {
                out.push(walker.clone());
            }
        }
        out
    }

    /// Gelman-Rubin R-hat per dimension over the post-burn chain, treating
    /// each walker as a chain (health metric; see module docs for the
    /// ensemble caveat).
    pub fn gelman_rubin(&self, burn_in: usize) -> Vec<f64> {
        let m = self.n_walkers();
        let n = self.n_steps().saturating_sub(burn_in);
        let d = self.dim();
        assert!(
            m >= 2 && n >= 4,
            "need >= 2 walkers and >= 4 post-burn steps"
        );

        (0..d)
            .map(|dim| {
                // Per-walker means and variances.
                let mut means = vec![0.0; m];
                for step in burn_in..self.n_steps() {
                    for (w, mean) in means.iter_mut().enumerate() {
                        *mean += self.samples[step][w][dim];
                    }
                }
                for mean in means.iter_mut() {
                    *mean /= n as f64;
                }
                let grand = means.iter().sum::<f64>() / m as f64;

                let mut within = 0.0;
                for (w, &mean) in means.iter().enumerate() {
                    let mut s2 = 0.0;
                    for step in burn_in..self.n_steps() {
                        s2 += (self.samples[step][w][dim] - mean).powi(2);
                    }
                    within += s2 / (n - 1) as f64;
                }
                within /= m as f64;

                let between = n as f64 / (m - 1) as f64
                    * means.iter().map(|&mu| (mu - grand).powi(2)).sum::<f64>();

                if within <= 0.0 {
                    return 1.0; // all walkers frozen at the same value
                }
                let var_plus = (n - 1) as f64 / n as f64 * within + between / n as f64;
                (var_plus / within).sqrt()
            })
            .collect()
    }

    /// Effective sample size per dimension: ESS = (m * n) / tau, with the
    /// integrated autocorrelation time tau estimated per walker, averaged
    /// across walkers, using Sokal's adaptive window (c = 5).
    pub fn effective_sample_size(&self, burn_in: usize) -> Vec<f64> {
        let m = self.n_walkers();
        let n = self.n_steps().saturating_sub(burn_in);
        let d = self.dim();
        assert!(n >= 8, "need >= 8 post-burn steps for an ESS estimate");

        (0..d)
            .map(|dim| {
                // Mean autocorrelation function across walkers.
                let max_lag = n / 2;
                let mut acf_sum = vec![0.0; max_lag];
                for w in 0..m {
                    let series: Vec<f64> = (burn_in..self.n_steps())
                        .map(|s| self.samples[s][w][dim])
                        .collect();
                    let mean = series.iter().sum::<f64>() / n as f64;
                    let var = series.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64;
                    if var <= 0.0 {
                        continue;
                    }
                    for (lag, acf) in acf_sum.iter_mut().enumerate() {
                        let mut c = 0.0;
                        for t in 0..(n - lag) {
                            c += (series[t] - mean) * (series[t + lag] - mean);
                        }
                        *acf += c / (n as f64 * var);
                    }
                }
                for acf in acf_sum.iter_mut() {
                    *acf /= m as f64;
                }

                // Sokal window: smallest W with W >= c * tau(W), c = 5.
                let mut tau = 1.0;
                for window in 1..max_lag {
                    tau = 1.0 + 2.0 * acf_sum[1..=window].iter().sum::<f64>();
                    if (window as f64) >= 5.0 * tau {
                        break;
                    }
                }
                let tau = tau.max(1.0);
                (m * n) as f64 / tau
            })
            .collect()
    }
}

/// Affine-invariant ensemble sampler (Goodman & Weare 2010 stretch move).
#[derive(Debug, Clone)]
pub struct EnsembleSampler {
    pub n_walkers: usize,
    /// Stretch parameter a (paper and emcee default: 2).
    pub stretch_a: f64,
}

impl EnsembleSampler {
    pub fn new(n_walkers: usize) -> Self {
        assert!(
            n_walkers >= 4 && n_walkers.is_multiple_of(2),
            "need an even number of walkers >= 4"
        );
        Self {
            n_walkers,
            stretch_a: 2.0,
        }
    }

    /// Run the sampler for `n_steps` ensemble updates from the given initial
    /// walker positions. `log_prob` may return -inf (rejected region).
    ///
    /// Red-black update: each half of the ensemble is moved using partners
    /// drawn only from the other half (Foreman-Mackey et al. 2013), which
    /// preserves detailed balance for the whole ensemble update.
    pub fn run<F: Fn(&[f64]) -> f64>(
        &self,
        log_prob: F,
        initial: &[Vec<f64>],
        n_steps: usize,
        seed: u64,
    ) -> EnsembleChain {
        assert_eq!(initial.len(), self.n_walkers);
        let dim = initial[0].len();
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

        let mut positions: Vec<Vec<f64>> = initial.to_vec();
        let mut log_probs: Vec<f64> = positions.iter().map(|p| log_prob(p)).collect();
        assert!(
            log_probs.iter().any(|lp| lp.is_finite()),
            "no walker starts at finite log-probability"
        );

        let half = self.n_walkers / 2;
        let mut samples = Vec::with_capacity(n_steps);
        let mut chain_log_probs = Vec::with_capacity(n_steps);
        let mut n_accept = 0usize;
        let mut n_propose = 0usize;

        for _ in 0..n_steps {
            for side in 0..2 {
                let (lo, hi) = if side == 0 {
                    (0, half)
                } else {
                    (half, self.n_walkers)
                };
                let (other_lo, other_len) = if side == 0 {
                    (half, self.n_walkers - half)
                } else {
                    (0, half)
                };
                for k in lo..hi {
                    let partner = other_lo + rng.gen_range(0..other_len);
                    // z ~ g(z) prop. 1/sqrt(z) on [1/a, a]:
                    // z = ((a-1) u + 1)^2 / a.
                    let u: f64 = rng.gen();
                    let z = ((self.stretch_a - 1.0) * u + 1.0).powi(2) / self.stretch_a;
                    let proposal: Vec<f64> = (0..dim)
                        .map(|d| {
                            positions[partner][d] + z * (positions[k][d] - positions[partner][d])
                        })
                        .collect();
                    let lp_new = log_prob(&proposal);
                    let ln_accept = (dim as f64 - 1.0) * z.ln() + lp_new - log_probs[k];
                    n_propose += 1;
                    if ln_accept >= 0.0 || rng.gen::<f64>().ln() < ln_accept {
                        positions[k] = proposal;
                        log_probs[k] = lp_new;
                        n_accept += 1;
                    }
                }
            }
            samples.push(positions.clone());
            chain_log_probs.push(log_probs.clone());
        }

        EnsembleChain {
            samples,
            log_probs: chain_log_probs,
            acceptance_rate: n_accept as f64 / n_propose as f64,
            seed,
        }
    }
}

/// Initialize walkers uniformly inside a box (seeded).
pub fn init_walkers_in_box(n_walkers: usize, bounds: &[(f64, f64)], seed: u64) -> Vec<Vec<f64>> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    (0..n_walkers)
        .map(|_| {
            bounds
                .iter()
                .map(|&(lo, hi)| rng.gen_range(lo..hi))
                .collect()
        })
        .collect()
}

/// Quantile of a (will be sorted) sample vector, linear interpolation.
pub fn quantile(values: &mut [f64], q: f64) -> f64 {
    assert!(!values.is_empty());
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pos = q.clamp(0.0, 1.0) * (values.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    let frac = pos - lo as f64;
    values[lo] * (1.0 - frac) + values[hi] * frac
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2D correlated Gaussian: mean (1, -2), var (2, 0.5), corr 0.6.
    fn gaussian_2d_logp(x: &[f64]) -> f64 {
        let (mx, my) = (1.0, -2.0);
        let (sx2, sy2, rho): (f64, f64, f64) = (2.0, 0.5, 0.6);
        let (sx, sy) = (sx2.sqrt(), sy2.sqrt());
        let dx = (x[0] - mx) / sx;
        let dy = (x[1] - my) / sy;
        -0.5 * (dx * dx - 2.0 * rho * dx * dy + dy * dy) / (1.0 - rho * rho)
    }

    #[test]
    fn test_recovers_gaussian_mean_and_covariance() {
        let sampler = EnsembleSampler::new(32);
        let init = init_walkers_in_box(32, &[(0.0, 2.0), (-3.0, -1.0)], 7);
        let chain = sampler.run(gaussian_2d_logp, &init, 3000, 11);

        let flat = chain.flat_samples(500, 1);
        let n = flat.len() as f64;
        let mean_x = flat.iter().map(|s| s[0]).sum::<f64>() / n;
        let mean_y = flat.iter().map(|s| s[1]).sum::<f64>() / n;
        assert!((mean_x - 1.0).abs() < 0.1, "mean x = {mean_x}");
        assert!((mean_y + 2.0).abs() < 0.05, "mean y = {mean_y}");

        let var_x = flat.iter().map(|s| (s[0] - mean_x).powi(2)).sum::<f64>() / n;
        let var_y = flat.iter().map(|s| (s[1] - mean_y).powi(2)).sum::<f64>() / n;
        let cov_xy = flat
            .iter()
            .map(|s| (s[0] - mean_x) * (s[1] - mean_y))
            .sum::<f64>()
            / n;
        assert!((var_x - 2.0).abs() < 0.25, "var x = {var_x}");
        assert!((var_y - 0.5).abs() < 0.07, "var y = {var_y}");
        let corr = cov_xy / (var_x * var_y).sqrt();
        assert!((corr - 0.6).abs() < 0.08, "corr = {corr}");

        // Healthy diagnostics on a well-sampled unimodal target.
        let r_hat = chain.gelman_rubin(500);
        assert!(r_hat.iter().all(|&r| r < 1.05), "R-hat = {r_hat:?}");
        let ess = chain.effective_sample_size(500);
        assert!(ess.iter().all(|&e| e > 500.0), "ESS = {ess:?}");
        assert!(chain.acceptance_rate > 0.1 && chain.acceptance_rate < 0.9);
    }

    #[test]
    fn test_rosenbrock_smoke() {
        // Rosenbrock density (log pi = -(100 (y - x^2)^2 + (1 - x)^2) / 20):
        // banana-shaped, the standard hard smoke target.
        let logp =
            |x: &[f64]| -(100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2)) / 20.0;
        let sampler = EnsembleSampler::new(32);
        let init = init_walkers_in_box(32, &[(-2.0, 2.0), (-1.0, 3.0)], 3);
        let chain = sampler.run(logp, &init, 6000, 5);

        let flat = chain.flat_samples(2000, 5);
        assert!(!flat.is_empty());
        // The mode is at (1, 1); the sampler must populate its vicinity.
        let near_mode = flat
            .iter()
            .filter(|s| (s[0] - 1.0).abs() < 1.0 && (s[1] - 1.0).abs() < 2.0)
            .count();
        assert!(
            near_mode as f64 > 0.2 * flat.len() as f64,
            "only {near_mode}/{} samples near the Rosenbrock mode",
            flat.len()
        );
        // Loose smoke bound: the banana's long tails give the stretch move
        // a large autocorrelation time (measured R-hat = [1.37, 1.29] at
        // 6000 steps / burn 2000, seed 5; 1.59 at 2000 steps — this is a
        // machinery check, not a convergence claim).
        let r_hat = chain.gelman_rubin(2000);
        assert!(
            r_hat.iter().all(|&r| r.is_finite() && r < 1.5),
            "R-hat = {r_hat:?}"
        );
    }

    #[test]
    fn test_sampler_is_deterministic_for_fixed_seed() {
        let sampler = EnsembleSampler::new(8);
        let init = init_walkers_in_box(8, &[(0.0, 1.0), (0.0, 1.0)], 1);
        let c1 = sampler.run(gaussian_2d_logp, &init, 50, 2);
        let c2 = sampler.run(gaussian_2d_logp, &init, 50, 2);
        assert_eq!(c1.samples, c2.samples);
        let c3 = sampler.run(gaussian_2d_logp, &init, 50, 3);
        assert_ne!(c1.samples, c3.samples);
    }

    #[test]
    fn test_walkers_never_enter_forbidden_region() {
        // Uniform box prior: -inf outside. No accepted sample may be outside.
        let logp = |x: &[f64]| {
            if x.iter().all(|&v| (0.0..1.0).contains(&v)) {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        };
        let sampler = EnsembleSampler::new(16);
        let init = init_walkers_in_box(16, &[(0.2, 0.8), (0.2, 0.8)], 9);
        let chain = sampler.run(logp, &init, 500, 13);
        for step in &chain.samples {
            for walker in step {
                assert!(walker.iter().all(|&v| (0.0..1.0).contains(&v)));
            }
        }
    }

    #[test]
    fn test_quantile() {
        let mut v: Vec<f64> = (0..101).map(|k| k as f64).collect();
        assert_eq!(quantile(&mut v, 0.5), 50.0);
        assert_eq!(quantile(&mut v, 0.0), 0.0);
        assert_eq!(quantile(&mut v, 1.0), 100.0);
        assert!((quantile(&mut v, 0.975) - 97.5).abs() < 1e-12);
    }

    #[test]
    fn test_paper_settings_pinned() {
        // Brown & Batygin (2021) Section 5: 100 chains, 20,890 samples,
        // burn-in 260, thin 42 -> 49,100 uncorrelated samples.
        let s = PaperMcmcSettings::paper();
        assert_eq!(s.n_walkers, 100);
        let kept_per_walker = (s.n_steps - s.burn_in) / s.thin;
        let total = kept_per_walker * s.n_walkers;
        assert!(
            (total as i64 - 49_100).abs() < 100,
            "paper sample count: {total}"
        );
    }
}

//! Film export for `p9-2026-apsidal-clustering`: the numbers its scene and ledger entry draw.
//!
//! The crate's 21- and 25-object samples are synthetic stand-ins tuned to the
//! published significances, so they are not drawn. The estimator is run on
//! the real vetted sample, before and after the two 2025 discoveries.

use p9_2025_new_discoveries::discoveries::{ammonite, of201};
use p9_2026_apsidal_clustering::estimator::{fit, log_likelihood_ratio};
use p9_2026_apsidal_clustering::samples::{stable_21, stable_25, vetted_etno_varpi};
use p9_2026_apsidal_clustering::significance::{lambda_to_p_value, sigma};
use p9_2026_apsidal_clustering::{PUBLISHED_SIGMA_21, PUBLISHED_SIGMA_25};
use p9_core::analysis::circular::bessel_i0;
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::{DEG2RAD, RAD2DEG, TWO_PI};
use serde_json::{Value, json};

/// One sample through the estimator: the fit, each object's contribution to
/// the log-likelihood ratio, both significance conventions, and the fitted
/// von Mises density for drawing.
fn analyse(varpi: &[f64]) -> Value {
    let f = fit(varpi);
    let p = lambda_to_p_value(f.lambda);
    let lon_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
    let norm = TWO_PI * bessel_i0(f.kappa);
    let density: Vec<f64> = lon_deg
        .iter()
        .map(|&l| (f.kappa * (l * DEG2RAD - f.mu).cos()).exp() / norm)
        .collect();
    let votes: Vec<f64> = varpi
        .iter()
        .map(|&w| log_likelihood_ratio(&[w], f.mu, f.kappa))
        .collect();
    json!({
        "n": f.n,
        "varpi_deg": varpi.iter().map(|v| v.rem_euclid(TWO_PI) * RAD2DEG).collect::<Vec<_>>(),
        "votes": votes,
        "mu_deg": f.mu.rem_euclid(TWO_PI) * RAD2DEG,
        "kappa": f.kappa,
        "r_bar": f.r_bar,
        "lambda": f.lambda,
        "p": p,
        "sigma_two_sided": sigma(f.lambda),
        "lon_deg": lon_deg,
        "density": density,
    })
}

pub fn export() -> Value {
    let real = vetted_etno_varpi();
    let mut real_plus = real.clone();
    real_plus.push(of201().etno.longitude_of_perihelion());
    real_plus.push(ammonite().etno.longitude_of_perihelion());

    let before = fit(&real);
    let after = fit(&real_plus);
    let p_after = lambda_to_p_value(after.lambda);

    json!({
        "real": analyse(&real),
        "real_plus": analyse(&real_plus),
        "sigma_real_two_sided": sigma(before.lambda),
        "sigma_real_plus_two_sided": sigma(after.lambda),
        "sigma": p_value_to_sigma(p_after),
        "paper_sigma21": PUBLISHED_SIGMA_21,
        "paper_sigma25": PUBLISHED_SIGMA_25,
        "paper_n_before": stable_21().len(),
        "paper_n_sample": stable_25().len(),
        "uniform_density": 1.0 / TWO_PI,
    })
}

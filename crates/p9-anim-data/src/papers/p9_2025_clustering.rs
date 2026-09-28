//! Film export for `p9-2025-clustering`: the numbers its scene and ledger entry draw.

use p9_2025_clustering::clone_generation::distant_tno_sample;
use p9_2025_clustering::clustering::{
    fit_von_mises, stable_clustering_paper, unstable_clustering_paper, von_mises_pdf,
};
use p9_2025_clustering::orbital_poles::{
    OrbitalPole, laplace_pole, mean_pole, misalignment, pole_scatter,
};
use p9_2025_clustering::pipeline::{PipelineConfig, run_pipeline};
use p9_2025_clustering::stability::{D_CRIT, StabilityClass};
use p9_core::analysis::circular::{circular_std, mean_resultant_length, rayleigh_p_value};
use p9_core::analysis::resonance::neptune_diffusion_coefficient;
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::{DEG2RAD, RAD2DEG, TWO_PI};
use serde_json::{Value, json};

/// Monte Carlo draws for the bias-resampling null (the smoke default is 200).
const N_MC: usize = 20_000;

fn class_name(c: StabilityClass) -> &'static str {
    match c {
        StabilityClass::Stable => "stable",
        StabilityClass::Metastable => "metastable",
        StabilityClass::Unstable => "unstable",
    }
}

pub fn export() -> Value {
    let config = PipelineConfig {
        n_mc: N_MC,
        ..PipelineConfig::smoke()
    };
    let result = run_pipeline(&config);
    let sample = distant_tno_sample();
    let laplace = laplace_pole();

    let objects: Vec<Value> = result
        .per_tno
        .iter()
        .zip(sample.iter())
        .map(|(r, tno)| {
            let pole = OrbitalPole::from_elements(tno.name, tno.i, tno.omega_big);
            json!({
                "name": r.name,
                "a": tno.a,
                "e": tno.e,
                "q": tno.a * (1.0 - tno.e),
                "i_deg": tno.i * RAD2DEG,
                "node_deg": (tno.omega_big * RAD2DEG).rem_euclid(360.0),
                "varpi_deg": (r.varpi * RAD2DEG).rem_euclid(360.0),
                "d_mean": r.d_mean,
                "d_std": r.d_std,
                "d_analytical": r.d_analytical,
                "class": class_name(r.class),
                "n_escaped": r.n_escaped,
                "pole_x_deg": pole.x * RAD2DEG,
                "pole_y_deg": pole.y * RAD2DEG,
                "pole_offset_deg": misalignment(&pole, &laplace) * RAD2DEG,
            })
        })
        .collect();

    // Mean orbital pole of the sample against the giant planets' plane.
    let poles: Vec<OrbitalPole> = sample
        .iter()
        .map(|t| OrbitalPole::from_elements(t.name, t.i, t.omega_big))
        .collect();
    let (mean_x, mean_y) = mean_pole(&poles);
    let mean = OrbitalPole::from_elements("mean pole", mean_x.hypot(mean_y), mean_y.atan2(mean_x));

    // Longitudes of perihelion of the objects the pipeline finds (meta)stable.
    let kept: Vec<f64> = result
        .per_tno
        .iter()
        .filter(|t| t.class != StabilityClass::Unstable)
        .map(|t| t.varpi)
        .collect();
    let fit = fit_von_mises(&kept);
    let paper_fit = stable_clustering_paper();
    let paper_unstable = unstable_clustering_paper();

    let lon_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
    let fit_density: Vec<f64> = lon_deg
        .iter()
        .map(|&l| von_mises_pdf(l * DEG2RAD, fit.mu, fit.kappa))
        .collect();
    let paper_density: Vec<f64> = lon_deg
        .iter()
        .map(|&l| von_mises_pdf(l * DEG2RAD, paper_fit.mu, paper_fit.kappa))
        .collect();

    // The analytical Neptune-scattering diffusion law against perihelion.
    let q_au: Vec<f64> = (0..=100).map(|k| 30.0 + 0.55 * k as f64).collect();
    let d_analytical: Vec<f64> = q_au
        .iter()
        .map(|&q| neptune_diffusion_coefficient(q))
        .collect();

    json!({
        "config": {
            "n_clones": config.n_clones,
            "t_total_yr": config.t_total_yr,
            "paper_n_clones": PipelineConfig::paper().n_clones,
            "paper_t_total_yr": PipelineConfig::paper().t_total_yr,
        },
        "objects": objects,
        "n_sample": result.per_tno.len(),
        "n_stable": result.summary.n_stable,
        "n_metastable": result.summary.n_metastable,
        "n_unstable": result.summary.n_unstable,
        "d_crit": D_CRIT,
        "diffusion_law": {"q_au": q_au, "d": d_analytical},
        "n_kept": kept.len(),
        "r_bar": mean_resultant_length(&kept),
        "rayleigh_p": rayleigh_p_value(&kept),
        "bias_mc_p": result.bias_mc_p,
        "sigma": p_value_to_sigma(result.bias_mc_p),
        "mean_varpi_deg": fit.mu.rem_euclid(TWO_PI) * RAD2DEG,
        "kappa": fit.kappa,
        "spread_rad": circular_std(&kept),
        "paper_mean_varpi_deg": paper_fit.mu * RAD2DEG,
        "paper_spread_rad": 1.0 / paper_fit.kappa.sqrt(),
        "paper_unstable_modes_deg": [
            paper_unstable.comp1.mu * RAD2DEG,
            paper_unstable.comp2.mu * RAD2DEG,
        ],
        "fit": {"lon_deg": lon_deg, "density": fit_density, "paper_density": paper_density},
        "laplace_pole": {
            "x_deg": laplace.x * RAD2DEG,
            "y_deg": laplace.y * RAD2DEG,
            "i_deg": laplace.x.hypot(laplace.y) * RAD2DEG,
        },
        "mean_pole": {
            "x_deg": mean_x * RAD2DEG,
            "y_deg": mean_y * RAD2DEG,
            "i_deg": mean_x.hypot(mean_y) * RAD2DEG,
            "offset_deg": misalignment(&mean, &laplace) * RAD2DEG,
            "scatter_deg": pole_scatter(&poles) * RAD2DEG,
        },
    })
}

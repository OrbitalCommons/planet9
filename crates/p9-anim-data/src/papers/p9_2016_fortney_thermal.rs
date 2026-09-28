//! Film export for `p9-2016-fortney-thermal`: the numbers its scene and ledger entry draw.

use p9_2016_fortney_thermal::published::{MASS_MAX_EARTH, MASS_MIN_EARTH, TEFF_MAX_K, TEFF_MIN_K};
use p9_2016_fortney_thermal::{P9Thermal, W1_WAVELENGTH_M, W2_WAVELENGTH_M};
use p9_core::analysis::thermal::R_EARTH_M;
use serde_json::{Value, json};

/// Heliocentric distance at which the spectrum and the mass scan are evaluated.
const DISTANCE_AU: f64 = 500.0;

/// Points on the log-wavelength grid of the spectrum (1 µm to 3 mm).
const N_SED: usize = 141;

/// Published excess of the model-atmosphere 3-5 µm flux over a blackbody
/// (orders of magnitude; Fortney et al. 2016 abstract).
const PUBLISHED_NEAR_IR_EXCESS_DEX: f64 = 20.0;

pub fn export() -> Value {
    let nominal = P9Thermal::new(10.0, DISTANCE_AU);

    let masses: Vec<f64> = (0..=36).map(|k| 3.0 + 0.5 * k as f64).collect();
    let t_eff: Vec<f64> = masses
        .iter()
        .map(|&m| P9Thermal::new(m, DISTANCE_AU).effective_temp())
        .collect();

    let distances: Vec<f64> = (0..=48).map(|k| 40.0 + 20.0 * k as f64).collect();
    let t_solar: Vec<f64> = distances
        .iter()
        .map(|&d| P9Thermal::new(10.0, d).solar_equilibrium_temp())
        .collect();
    let t_total: Vec<f64> = distances
        .iter()
        .map(|&d| P9Thermal::new(10.0, d).effective_temp())
        .collect();

    let log_um: Vec<f64> = (0..N_SED)
        .map(|k| 3.5 * k as f64 / (N_SED - 1) as f64)
        .collect();
    let log_flux: Vec<f64> = log_um
        .iter()
        .map(|&lg| {
            nominal
                .flux_density_jy(10f64.powf(lg) * 1.0e-6)
                .max(1.0e-300)
                .log10()
        })
        .collect();

    let log_w1 = nominal.w1_flux_jy().max(1.0e-300).log10();
    let log_w2 = nominal.w2_flux_jy().max(1.0e-300).log10();

    json!({
        "distance_au": DISTANCE_AU,
        "teff_10me": nominal.effective_temp(),
        "teff_5me": P9Thermal::new(MASS_MIN_EARTH, DISTANCE_AU).effective_temp(),
        "teff_20me": P9Thermal::new(MASS_MAX_EARTH, DISTANCE_AU).effective_temp(),
        "t_solar_10me": nominal.solar_equilibrium_temp(),
        "radius_10me_earth": nominal.radius_m() / R_EARTH_M,
        "internal_to_solar": nominal.internal_luminosity_w() / nominal.absorbed_sunlight_w(),
        "published": {
            "teff_min_k": TEFF_MIN_K,
            "teff_max_k": TEFF_MAX_K,
            "mass_min_earth": MASS_MIN_EARTH,
            "mass_max_earth": MASS_MAX_EARTH,
            "near_ir_excess_dex": PUBLISHED_NEAR_IR_EXCESS_DEX,
        },
        "teff_vs_mass": {"mass_earth": masses, "t_eff_k": t_eff},
        "temp_vs_distance": {
            "distance_au": distances,
            "solar_only_k": t_solar,
            "with_internal_heat_k": t_total,
        },
        "sed": {"log_wavelength_um": log_um, "log_flux_jy": log_flux},
        "sed_peak_um": nominal.sed_peak_wavelength_m() * 1.0e6,
        "far_ir_flux_jy": nominal.far_ir_flux_jy(),
        "w1": {"wavelength_um": W1_WAVELENGTH_M * 1.0e6, "log_flux_jy": log_w1},
        "w2": {"wavelength_um": W2_WAVELENGTH_M * 1.0e6, "log_flux_jy": log_w2},
    })
}

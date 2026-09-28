//! Film export for `p9-2016-cowan-thermal`: the numbers its scene and ledger entry draw.

use p9_2016_cowan_thermal::sed::{
    COWAN_FAINT_FLUX_1MM_MJY, COWAN_FLUX_1MM_MJY, ONE_MM_M, P9Sed, sensitivity::CMB_1MM_MJY,
};
use p9_core::analysis::thermal::{C_LIGHT, H_PLANCK, K_BOLTZ, R_EARTH_M, WIEN_X_NU};
use serde_json::{Value, json};

/// Points on the log-wavelength grid of the spectrum (0.3 µm to 1 cm).
const N_SED: usize = 181;

/// CMB-experiment band centres (GHz) marked on the spectrum.
const CMB_BANDS_GHZ: [f64; 3] = [98.0, 150.0, 229.0];

/// A smaller, colder body: the faint end of the paper's range.
const FAINT_MASS_EARTH: f64 = 5.0;
const FAINT_TEMP_K: f64 = 30.0;

pub fn export() -> Value {
    let fid = P9Sed::cowan_fiducial();

    let log_um: Vec<f64> = (0..N_SED)
        .map(|k| -0.5 + 4.5 * k as f64 / (N_SED - 1) as f64)
        .collect();
    let lambda_m = |lg: f64| 10f64.powf(lg) * 1.0e-6;
    let thermal: Vec<f64> = log_um
        .iter()
        .map(|&lg| fid.thermal_flux_mjy(lambda_m(lg)).max(1.0e-300).log10())
        .collect();
    let reflected: Vec<f64> = log_um
        .iter()
        .map(|&lg| fid.reflected_flux_mjy(lambda_m(lg)).max(1.0e-300).log10())
        .collect();

    let cmb_bands: Vec<Value> = CMB_BANDS_GHZ
        .iter()
        .map(|&ghz| {
            let lam = C_LIGHT / (ghz * 1.0e9);
            json!({
                "ghz": ghz,
                "wavelength_um": lam * 1.0e6,
                "flux_mjy": fid.total_flux_mjy(lam),
            })
        })
        .collect();

    // Flux at 1 mm against distance: the fiducial body and the faint case.
    let distances: Vec<f64> = (0..=60).map(|k| 200.0 + 20.0 * k as f64).collect();
    let flux_at = |mass: f64, temp: f64| -> Vec<f64> {
        distances
            .iter()
            .map(|&d| P9Sed::new(mass, d, temp).total_flux_mjy(ONE_MM_M))
            .collect()
    };

    json!({
        "mass_earth": fid.mass_earth,
        "distance_au": fid.distance_au,
        "temp_k": fid.temp_k,
        "radius_earth": fid.radius_m() / R_EARTH_M,
        "flux_1mm_mjy": fid.total_flux_mjy(ONE_MM_M),
        "published_flux_1mm_mjy": COWAN_FLUX_1MM_MJY,
        "published_faint_flux_1mm_mjy": COWAN_FAINT_FLUX_1MM_MJY,
        "cmb_threshold_mjy": CMB_1MM_MJY,
        // The plotted spectrum is per unit frequency, so mark its B_nu peak.
        "wien_peak_um": 1.0e6 * C_LIGHT * H_PLANCK / (WIEN_X_NU * K_BOLTZ * fid.temp_k),
        "crossover_um": fid.crossover_wavelength_m().map(|m| m * 1.0e6),
        "v_band_flux_mjy": fid.total_flux_mjy(0.55e-6),
        "parallax_arcmin": (1.0 / fid.distance_au).atan().to_degrees() * 60.0,
        "sed": {
            "log_wavelength_um": log_um,
            "log_thermal_mjy": thermal,
            "log_reflected_mjy": reflected,
        },
        "cmb_bands": cmb_bands,
        "flux_vs_distance": {
            "distance_au": distances,
            "fiducial_mjy": flux_at(fid.mass_earth, fid.temp_k),
            "faint_mjy": flux_at(FAINT_MASS_EARTH, FAINT_TEMP_K),
            "faint_mass_earth": FAINT_MASS_EARTH,
            "faint_temp_k": FAINT_TEMP_K,
        },
    })
}

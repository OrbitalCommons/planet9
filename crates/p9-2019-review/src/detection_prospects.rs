//! Detection prospects for revised Planet Nine parameters.
//!
//! Computes V-band magnitude estimates for different P9 configurations
//! and compares with survey depth limits.

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::types::P9Params;

/// Physical properties of Planet Nine for brightness estimation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct P9PhysicalProperties {
    /// Physical radius in Earth radii
    pub radius_earth: f64,
    /// Geometric albedo (0 to 1)
    pub albedo: f64,
}

impl P9PhysicalProperties {
    /// Conservative estimate for 5 M_Earth.
    pub fn five_me_conservative() -> Self {
        Self {
            radius_earth: 2.4,
            albedo: 0.40,
        }
    }

    /// Optimistic estimate for 5 M_Earth.
    pub fn five_me_optimistic() -> Self {
        Self {
            radius_earth: 3.5,
            albedo: 0.75,
        }
    }

    /// Conservative estimate for 10 M_Earth.
    pub fn ten_me_conservative() -> Self {
        Self {
            radius_earth: 1.9,
            albedo: 0.40,
        }
    }

    /// Optimistic estimate for 10 M_Earth.
    pub fn ten_me_optimistic() -> Self {
        Self {
            radius_earth: 3.7,
            albedo: 0.75,
        }
    }
}

/// Compute apparent V-band magnitude at a given heliocentric distance.
///
/// Uses the standard asteroid photometry formula:
///   V = H + 5 * log10(r * delta) - 2.5 * log10(phase_function)
///
/// where H is the absolute magnitude (r=1 AU, delta=1 AU, zero phase angle).
///
/// For distant objects (r >> 1 AU), delta ≈ r, so:
///   V ≈ H + 10 * log10(r)
pub fn apparent_magnitude(r_au: f64, physical: &P9PhysicalProperties) -> f64 {
    // Physical radius in AU
    let r_earth_au = 4.26e-5_f64; // Earth radius in AU
    let radius_au = physical.radius_earth * r_earth_au;

    // Absolute magnitude from radius and albedo
    // H = -5 * log10(diameter_km) - 2.5 * log10(albedo) + 15.618
    let diameter_km = radius_au * 2.0 * AU_KM;
    let h = -5.0 * diameter_km.log10() - 2.5 * physical.albedo.log10() + 15.618;

    // At large distances, geocentric distance ≈ heliocentric distance
    let delta = r_au;

    h + 5.0 * (r_au * delta).log10()
}

/// Compute V magnitude at perihelion and aphelion.
pub fn brightness_at_extremes(params: &P9Params, physical: &P9PhysicalProperties) -> (f64, f64) {
    let q = params.a * (1.0 - params.e);
    let aphelion = params.a * (1.0 + params.e);

    let v_perihelion = apparent_magnitude(q, physical);
    let v_aphelion = apparent_magnitude(aphelion, physical);

    (v_perihelion, v_aphelion)
}

/// Survey depth limits.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SurveyLimit {
    pub name: &'static str,
    pub v_limit: f64,
}

/// Known survey depth limits.
pub fn survey_limits() -> Vec<SurveyLimit> {
    vec![
        SurveyLimit {
            name: "Pan-STARRS",
            v_limit: 21.5,
        },
        SurveyLimit {
            name: "DECam/Blanco",
            v_limit: 23.5,
        },
        SurveyLimit {
            name: "Subaru HSC",
            v_limit: 25.0,
        },
        SurveyLimit {
            name: "LSST (Rubin)",
            v_limit: 26.5,
        },
    ]
}

/// Brightness estimate for the revised P9 parameters.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BrightnessEstimate {
    pub label: String,
    pub mass_earth: f64,
    pub v_perihelion_bright: f64,
    pub v_perihelion_faint: f64,
    pub v_aphelion_bright: f64,
    pub v_aphelion_faint: f64,
}

/// Compute brightness estimates for the paper's best-fit parameters.
pub fn brightness_table() -> Vec<BrightnessEstimate> {
    use crate::revised_parameters::{best_fit_10me, best_fit_5me};

    let mut estimates = Vec::new();

    for params in best_fit_5me() {
        let (v_peri_con, v_aph_con) =
            brightness_at_extremes(&params, &P9PhysicalProperties::five_me_conservative());
        let (v_peri_opt, v_aph_opt) =
            brightness_at_extremes(&params, &P9PhysicalProperties::five_me_optimistic());

        estimates.push(BrightnessEstimate {
            label: format!("5 ME, a={}, e={}", params.a, params.e),
            mass_earth: 5.0,
            v_perihelion_bright: v_peri_opt,
            v_perihelion_faint: v_peri_con,
            v_aphelion_bright: v_aph_opt,
            v_aphelion_faint: v_aph_con,
        });
    }

    for params in best_fit_10me() {
        let (v_peri_con, v_aph_con) =
            brightness_at_extremes(&params, &P9PhysicalProperties::ten_me_conservative());
        let (v_peri_opt, v_aph_opt) =
            brightness_at_extremes(&params, &P9PhysicalProperties::ten_me_optimistic());

        estimates.push(BrightnessEstimate {
            label: format!("10 ME, a={}, e={}", params.a, params.e),
            mass_earth: 10.0,
            v_perihelion_bright: v_peri_opt,
            v_perihelion_faint: v_peri_con,
            v_aphelion_bright: v_aph_opt,
            v_aphelion_faint: v_aph_con,
        });
    }

    estimates
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_apparent_magnitude_nearby() {
        let physical = P9PhysicalProperties::five_me_optimistic();
        let v_100 = apparent_magnitude(100.0, &physical);
        let v_500 = apparent_magnitude(500.0, &physical);

        // Should be fainter at larger distances
        assert!(
            v_500 > v_100,
            "V at 500 AU ({:.1}) should be fainter than at 100 AU ({:.1})",
            v_500,
            v_100
        );
    }

    #[test]
    fn test_brightness_at_extremes() {
        let params = P9Params {
            mass_earth: 5.0,
            a: 500.0,
            e: 0.25,
            i: 20.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        let physical = P9PhysicalProperties::five_me_optimistic();
        let (v_peri, v_aph) = brightness_at_extremes(&params, &physical);

        assert!(v_peri < v_aph, "Should be brighter at perihelion");
        // Paper reports V ~ 19-22 for 5 ME
        assert!(
            v_peri > 15.0 && v_peri < 25.0,
            "V perihelion = {:.1}",
            v_peri
        );
    }

    #[test]
    fn test_brightness_table() {
        let table = brightness_table();
        assert_eq!(table.len(), 4); // 2 for 5 ME + 2 for 10 ME

        for entry in &table {
            assert!(entry.v_perihelion_bright < entry.v_aphelion_bright);
            assert!(entry.v_perihelion_bright < entry.v_perihelion_faint);
        }
    }

    #[test]
    fn test_survey_limits() {
        let limits = survey_limits();
        assert_eq!(limits.len(), 4);
        // LSST should be deepest
        assert!(limits[3].v_limit > limits[0].v_limit);
    }
}

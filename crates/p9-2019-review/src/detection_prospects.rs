//! Detection prospects for revised Planet Nine parameters.
//!
//! Computes V-band magnitude estimates for different P9 configurations
//! and compares with survey depth limits. Photometry goes through
//! `p9_core::analysis::photometry` (single workspace implementation of the
//! reflected-sunlight distance law and absolute magnitude); survey depths
//! shared with `p9_core::analysis::surveys` where available.

use p9_core::analysis::photometry;
use p9_core::analysis::surveys::limiting_magnitude;
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

/// Compute apparent V-band magnitude at a given heliocentric distance,
/// observed at opposition, via the shared p9-core photometry:
///   H = 5 log10(1329 / (√p · D_km)),  m = H + 5 log10(r·Δ),  Δ = r − 1 AU.
pub fn apparent_magnitude(r_au: f64, physical: &P9PhysicalProperties) -> f64 {
    let radius_km = physical.radius_earth * EARTH_RADIUS_KM;
    let h = photometry::absolute_magnitude(radius_km, physical.albedo);
    photometry::apparent_magnitude(h, r_au, photometry::opposition_delta(r_au))
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

/// Known survey depth limits. Pan-STARRS and DES come from the shared
/// p9-core survey table; Subaru HSC and LSST are not in the shared table
/// (deep targeted surveys) and are kept here.
pub fn survey_limits() -> Vec<SurveyLimit> {
    vec![
        SurveyLimit {
            name: "Pan-STARRS",
            v_limit: limiting_magnitude("PS1 3pi").expect("PS1 in shared table"),
        },
        SurveyLimit {
            name: "DECam/DES",
            v_limit: limiting_magnitude("DES").expect("DES in shared table"),
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

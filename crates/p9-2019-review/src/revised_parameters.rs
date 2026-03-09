//! Revised Planet Nine parameter estimates from the 2019 review.
//!
//! Compares original Batygin & Brown (2016) parameters with the refined
//! estimates from the comprehensive simulation ensemble.

use p9_core::constants::*;
use p9_core::types::P9Params;

/// Parameter range for a single orbital element.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct ParameterRange {
    pub min: f64,
    pub max: f64,
    pub best: f64,
}

/// Complete set of parameter ranges for Planet Nine.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct P9ParameterSet {
    pub label: &'static str,
    pub mass_earth: ParameterRange,
    pub a: ParameterRange,
    pub e: ParameterRange,
    pub i_deg: ParameterRange,
}

/// Original Batygin & Brown (2016) parameters.
pub fn original_2016() -> P9ParameterSet {
    P9ParameterSet {
        label: "Batygin & Brown (2016)",
        mass_earth: ParameterRange {
            min: 10.0,
            max: 10.0,
            best: 10.0,
        },
        a: ParameterRange {
            min: 700.0,
            max: 700.0,
            best: 700.0,
        },
        e: ParameterRange {
            min: 0.6,
            max: 0.6,
            best: 0.6,
        },
        i_deg: ParameterRange {
            min: 30.0,
            max: 30.0,
            best: 30.0,
        },
    }
}

/// Revised parameters from the 2019 review.
pub fn revised_2019() -> P9ParameterSet {
    P9ParameterSet {
        label: "Batygin+ (2019)",
        mass_earth: ParameterRange {
            min: 5.0,
            max: 10.0,
            best: 5.0,
        },
        a: ParameterRange {
            min: 400.0,
            max: 800.0,
            best: 500.0,
        },
        e: ParameterRange {
            min: 0.2,
            max: 0.5,
            best: 0.25,
        },
        i_deg: ParameterRange {
            min: 15.0,
            max: 25.0,
            best: 20.0,
        },
    }
}

/// Best-fit parameters for m₉ = 5 M_Earth (Table in Section 5).
pub fn best_fit_5me() -> Vec<P9Params> {
    vec![
        P9Params {
            mass_earth: 5.0,
            a: 400.0,
            e: 0.15,
            i: 20.0 * DEG2RAD,
            omega: 140.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        },
        P9Params {
            mass_earth: 5.0,
            a: 500.0,
            e: 0.25,
            i: 20.0 * DEG2RAD,
            omega: 140.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        },
    ]
}

/// Best-fit parameters for m₉ = 10 M_Earth.
pub fn best_fit_10me() -> Vec<P9Params> {
    vec![
        P9Params {
            mass_earth: 10.0,
            a: 700.0,
            e: 0.35,
            i: 20.0 * DEG2RAD,
            omega: 140.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        },
        P9Params {
            mass_earth: 10.0,
            a: 800.0,
            e: 0.45,
            i: 15.0 * DEG2RAD,
            omega: 140.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        },
    ]
}

/// Compute perihelion distance from parameters.
pub fn perihelion(a: f64, e: f64) -> f64 {
    a * (1.0 - e)
}

/// Compute aphelion distance from parameters.
pub fn aphelion(a: f64, e: f64) -> f64 {
    a * (1.0 + e)
}

/// Check if a set of P9 parameters falls within the allowed region.
///
/// Constraints: q₉ > 100 AU (too close → detected) and q₉ < 500 AU (too far → no effect).
pub fn is_viable(params: &P9Params) -> bool {
    let q = perihelion(params.a, params.e);
    q > 100.0 && q < 500.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_original_vs_revised() {
        let orig = original_2016();
        let rev = revised_2019();

        // Revised mass is lower
        assert!(rev.mass_earth.best <= orig.mass_earth.best);

        // Revised a range is broader and includes smaller values
        assert!(rev.a.min < orig.a.best);
    }

    #[test]
    fn test_best_fit_viable() {
        for params in best_fit_5me() {
            assert!(is_viable(&params), "5 ME fit should be viable");
        }
        for params in best_fit_10me() {
            assert!(is_viable(&params), "10 ME fit should be viable");
        }
    }

    #[test]
    fn test_perihelion_aphelion() {
        let p = perihelion(500.0, 0.25);
        let a = aphelion(500.0, 0.25);
        assert!((p - 375.0).abs() < 1e-10);
        assert!((a - 625.0).abs() < 1e-10);
    }
}

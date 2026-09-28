//! Perihelion clustering strength (paper Section 3.2): the von Mises
//! concentration κ of the longitudes of perihelion ϖ of the clustering
//! sample (a ≥ 250 AU, q ∈ (40, 100) AU, i ≤ 40°), κ = A⁻¹(R̄) with R̄ the
//! mean resultant length of the ϖ.

use p9_core::analysis::circular::{kappa_from_r_bar, mean_resultant_length};
use p9_core::types::OrbitalElements;

use crate::width::{select, SelectionCuts};

/// von Mises κ of the longitudes of perihelion of the orbits passing
/// [`SelectionCuts::clustering_sample`]. Returns NaN with fewer than three
/// orbits.
pub fn perihelion_concentration(elements: &[OrbitalElements]) -> f64 {
    let sel = select(elements, &SelectionCuts::clustering_sample());
    if sel.len() < 3 {
        return f64::NAN;
    }
    let varpi: Vec<f64> = sel.iter().map(|e| e.longitude_of_perihelion()).collect();
    kappa_from_r_bar(mean_resultant_length(&varpi))
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    fn orbit(varpi_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: 400.0,
            e: 0.85,
            i: 10.0 * DEG2RAD,
            omega_big: 0.0,
            omega: varpi_deg * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    #[test]
    fn aligned_apsides_are_concentrated_and_spread_ones_are_not() {
        let tight: Vec<_> = (0..20).map(|k| orbit(60.0 + k as f64)).collect();
        let spread: Vec<_> = (0..20).map(|k| orbit(18.0 * k as f64)).collect();
        assert!(perihelion_concentration(&tight) > 5.0);
        assert!(perihelion_concentration(&spread) < 0.5);
    }
}

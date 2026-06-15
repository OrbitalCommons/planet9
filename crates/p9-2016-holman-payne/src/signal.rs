//! Earth–Saturn range-residual amplitude induced by Planet Nine, expressed in
//! metres for direct comparison with the Cassini range precision.
//!
//! The underlying physics — the differential (tidal) acceleration P9 exerts on
//! Saturn relative to the Sun, propagated over Saturn's real Cassini-epoch arc
//! and reduced to a post-fit range residual — is computed entirely by the
//! sibling [`p9_2016_cassini_ranging`] crate. Here we only re-express it in the
//! units Holman & Payne use (metres) and expose the favored-zone selection.

use p9_core::types::P9Params;

use p9_2016_cassini_ranging::perturbation::{
    favored_true_anomaly as sibling_favored, range_perturbation_amplitude,
};

/// Metres per kilometre.
const M_PER_KM: f64 = 1_000.0;

/// Post-fit Earth–Saturn range-residual amplitude (METRES) for P9 at true
/// anomaly `nu` (radians) on the orbit `params`.
///
/// Thin unit wrapper over the sibling crate's
/// [`range_perturbation_amplitude`] (which returns kilometres). The amplitude
/// scales LINEARLY with P9 mass (through `GM_p9`) and falls STEEPLY (≈`1/d³`)
/// with P9 heliocentric distance — the two scalings Holman & Payne emphasize.
pub fn range_residual_m(params: &P9Params, nu: f64) -> f64 {
    range_perturbation_amplitude(params, nu) * M_PER_KM
}

/// The favored true anomaly (radians) on the near (post-perihelion) arc — the
/// detectable post-fit minimum, identical to the sibling crate's selection
/// logic. `n` is the number of uniform ν samples scanned.
pub fn favored_true_anomaly(params: &P9Params, n: usize) -> f64 {
    sibling_favored(params, n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use p9_core::data::ephemeris_constraint::brown_batygin_orbit;
    use p9_core::types::position_at_true_anomaly;

    #[test]
    fn residual_scales_linearly_with_mass() {
        // amplitude ∝ GM_p9 ∝ mass, at fixed geometry.
        let mut p10 = brown_batygin_orbit();
        p10.mass_earth = 10.0;
        let mut p20 = brown_batygin_orbit();
        p20.mass_earth = 20.0;
        let nu = 20.0_f64.to_radians(); // perihelion-facing, high signal
        let r10 = range_residual_m(&p10, nu);
        let r20 = range_residual_m(&p20, nu);
        assert!((r20 / r10 - 2.0).abs() < 1e-9, "ratio = {}", r20 / r10);
    }

    #[test]
    fn residual_falls_steeply_with_distance() {
        // Push P9 farther by enlarging the orbit; the tidal range signal should
        // collapse much faster than linearly (≈1/d³). Compare at a fixed ν on
        // the high-signal perihelion-facing arc, where the distance ratio at
        // that phase tracks the semi-major-axis ratio.
        let nu = 20.0_f64.to_radians();
        let mut p_close = brown_batygin_orbit();
        p_close.a = 400.0;
        let mut p_far = brown_batygin_orbit();
        p_far.a = 800.0;
        let close = range_residual_m(&p_close, nu);
        let far = range_residual_m(&p_far, nu);
        let r_close = position_at_true_anomaly(&p_close, nu).distance;
        let r_far = position_at_true_anomaly(&p_far, nu).distance;
        // Effective falloff exponent from the two samples.
        let n_eff = (close / far).ln() / (r_far / r_close).ln();
        assert!(
            n_eff > 2.5,
            "distance falloff too shallow: n_eff = {n_eff:.2} (r {r_close:.0}→{r_far:.0} AU)"
        );
    }

    #[test]
    fn residual_is_metres_scale_of_sibling_km() {
        // Exactly 1000× the sibling km amplitude (unit conversion only).
        let p = brown_batygin_orbit();
        let nu = 30.0_f64.to_radians();
        let km = range_perturbation_amplitude(&p, nu);
        assert!((range_residual_m(&p, nu) - km * 1000.0).abs() < 1e-12);
    }

    #[test]
    fn favored_anomaly_near_published() {
        let p = brown_batygin_orbit();
        let nu_deg = favored_true_anomaly(&p, 1440).to_degrees();
        let diff = (nu_deg - published::PREFERRED_TRUE_ANOMALY_DEG).abs();
        assert!(
            diff <= published::TRUE_ANOMALY_TOLERANCE_DEG,
            "favored ν = {nu_deg:.1}° not within ±{}° of {}°",
            published::TRUE_ANOMALY_TOLERANCE_DEG,
            published::PREFERRED_TRUE_ANOMALY_DEG
        );
    }
}

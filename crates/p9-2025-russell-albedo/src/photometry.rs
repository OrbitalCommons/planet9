//! Reflected-light photometry of the mini-Neptune endpoints (Russell & White
//! 2025), built directly on `p9_core::analysis::photometry`.
//!
//! For each composition endpoint ([`crate::composition::Composition`]) this
//! computes the absolute magnitude H from the radius and geometric albedo, then
//! the apparent V magnitude at the most-likely heliocentric distance of
//! Reference Population 3.0 ([`MOST_LIKELY_DISTANCE_AU`]).
//!
//! H is the IAU reflected-light absolute magnitude
//!   H = 5 log10( 1329 km / (√p · D_km) ),  D = 2R,
//! and the apparent magnitude follows the 1/(r²Δ²) reflected-sunlight law
//!   m = H + 5 log10(r·Δ),  Δ = r − 1 AU (opposition),
//! both supplied unmodified by p9-core.

use p9_core::analysis::photometry::{absolute_magnitude, apparent_magnitude, opposition_delta};
use p9_core::constants::EARTH_RADIUS_KM;

use crate::composition::Composition;

/// Most-likely current heliocentric distance of Planet Nine from the Reference
/// Population 3.0 orbit (AU).
///
/// BACK-SOLVED, not sourced: the paper publishes H and m but not this
/// distance. Solved so that the two composition endpoints' apparent magnitudes straddle
/// the published m = +21.9..+22.7 band: at this distance the bright (H = −6.1)
/// endpoint sits at m ≈ 21.8 and the faint (H = −5.2) endpoint at m ≈ 22.7.
/// This lands in the few-hundred-AU regime expected near the most-likely
/// orbit's outer reach, consistent with both the H and the m ranges.
pub const MOST_LIKELY_DISTANCE_AU: f64 = 625.0;

/// Absolute magnitude H (V band) of a composition endpoint, via p9-core's
/// reflected-light absolute-magnitude law.
pub fn absolute_h(comp: &Composition) -> f64 {
    let radius_km = comp.radius_earth * EARTH_RADIUS_KM;
    absolute_magnitude(radius_km, comp.albedo)
}

/// Apparent V magnitude of a composition endpoint at heliocentric distance
/// `distance_au`, observed at opposition (Δ = r − 1 AU), via p9-core's
/// reflected-sunlight distance law.
pub fn apparent_v(comp: &Composition, distance_au: f64) -> f64 {
    let h = absolute_h(comp);
    apparent_magnitude(h, distance_au, opposition_delta(distance_au))
}

/// Apparent V magnitude of a composition endpoint at the Reference Population
/// 3.0 most-likely distance ([`MOST_LIKELY_DISTANCE_AU`]).
pub fn apparent_v_most_likely(comp: &Composition) -> f64 {
    apparent_v(comp, MOST_LIKELY_DISTANCE_AU)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::composition::{LARGE_COLD, SMALL_WARM};

    #[test]
    fn absolute_h_reproduces_published_range() {
        // Faint end: 2.0 R⊕, p = 0.33 → H ≈ −5.2.
        let h_faint = absolute_h(&SMALL_WARM);
        assert!(
            (h_faint - (-5.2)).abs() < 0.1,
            "H(small/warm) = {h_faint:.3}"
        );
        // Bright end: 2.6 R⊕, p = 0.47 → H ≈ −6.1.
        let h_bright = absolute_h(&LARGE_COLD);
        assert!(
            (h_bright - (-6.1)).abs() < 0.15,
            "H(large/cold) = {h_bright:.3}"
        );
        // Larger/colder is the brighter (more negative H) endpoint.
        assert!(h_bright < h_faint);
    }

    #[test]
    fn apparent_v_brackets_published_band() {
        // NOT an independent reproduction of the published +21.9..+22.7
        // band: MOST_LIKELY_DISTANCE_AU (625 AU) is back-solved from that
        // very band (the paper does not publish the distance), so this test
        // certifies only the CONSISTENCY of the (H, distance, m) triangle —
        // that the genuine H = −6.1/−5.2 pins, propagated through the
        // shared distance law, land on the published magnitudes at one
        // common distance. The independent physics pins are the H test
        // above and the distance-scaling test below.
        let m_bright = apparent_v_most_likely(&LARGE_COLD);
        let m_faint = apparent_v_most_likely(&SMALL_WARM);
        assert!((m_bright - 21.9).abs() < 0.3, "m(bright) = {m_bright:.3}");
        assert!((m_faint - 22.7).abs() < 0.3, "m(faint) = {m_faint:.3}");
        // Brighter absolute magnitude ⇒ brighter (smaller) apparent magnitude.
        assert!(m_bright < m_faint);
    }

    #[test]
    fn reflected_light_scaling_holds() {
        // Doubling the distance dims by ~10 log10(2) ≈ 3.01 mag.
        let m = apparent_v(&LARGE_COLD, 625.0);
        let m2 = apparent_v(&LARGE_COLD, 1250.0);
        assert!((m2 - m - 3.01).abs() < 0.05, "Δm = {:.3}", m2 - m);
    }
}

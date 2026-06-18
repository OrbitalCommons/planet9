//! Place 2015 BP519 in the clustered-ϖ context of the ETNO sample.
//!
//! Becker et al. (2018) note that BP519 "would represent the first member of
//! the population of high-i, ϖ-shepherded TNOs" predicted by Planet Nine. We
//! quantify that by comparing BP519's longitude of perihelion ϖ = ω + Ω to
//! the circular mean ϖ of the Brown (2017) clustering sample
//! (`p9_core::data::etno::BROWN_2017_SAMPLE`), using the shared circular
//! statistics in `p9_core::analysis::circular`.

use crate::bp519::Bp519;
use p9_core::analysis::circular::{circular_mean, mean_resultant_length, wrap_to_pi};
use p9_core::constants::RAD2DEG;
use p9_core::data::etno::longitudes_of_perihelion;
use uom::si::angle::radian;

/// Circular mean ϖ of the Brown (2017) sample, in radians [0, 2π).
pub fn sample_mean_varpi() -> f64 {
    circular_mean(&longitudes_of_perihelion()).expect("sample has a defined mean direction")
}

/// Mean resultant length R̄ of the Brown (2017) sample's ϖ.
pub fn sample_varpi_concentration() -> f64 {
    mean_resultant_length(&longitudes_of_perihelion())
}

/// Signed angular offset of BP519's ϖ from the sample mean, wrapped to
/// (−180°, 180°] and returned in degrees. Near 0° means apsidally aligned
/// with the cluster.
pub fn varpi_offset_from_cluster_deg(bp: &Bp519) -> f64 {
    let delta =
        wrap_to_pi(bp.longitude_of_perihelion_typed().get::<radian>() - sample_mean_varpi());
    delta * RAD2DEG
}

/// Mean resultant length of the sample's ϖ *with* BP519 appended — a quick
/// check that adding BP519 does not destroy the clustering (it should remain
/// comparable if BP519 is consistent with the population).
pub fn concentration_with_bp519(bp: &Bp519) -> f64 {
    let mut varpis = longitudes_of_perihelion();
    varpis.push(bp.longitude_of_perihelion_typed().get::<radian>());
    mean_resultant_length(&varpis)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bp519::discovery_2018;

    #[test]
    fn sample_is_clustered() {
        // Sanity: the Brown 2017 ϖ sample is clustered (R̄ well above the
        // isotropic ~1/√10 ≈ 0.32 noise floor).
        assert!(
            sample_varpi_concentration() > 0.3,
            "R̄ = {:.3}",
            sample_varpi_concentration()
        );
    }

    #[test]
    fn bp519_varpi_offset_is_finite_and_reported() {
        // BP519's ϖ = ω + Ω = 348° + 135° = 483° ≡ 123°. The sample mean is
        // computed; we report the signed offset honestly rather than assert a
        // tight alignment (BP519 sits in the broad cluster, not its core, and
        // its ω is poorly constrained at discovery).
        let bp = discovery_2018();
        let off = varpi_offset_from_cluster_deg(&bp);
        assert!(off.abs() <= 180.0, "offset out of range: {off:.1}°");
        // It is on the clustered side of the circle, not anti-clustered:
        // within a half-circle of the sample mean.
        assert!(
            off.abs() < 90.0,
            "BP519 ϖ should sit within the clustered half-circle, off = {off:.1}°"
        );
    }

    #[test]
    fn adding_bp519_keeps_clustering() {
        let bp = discovery_2018();
        let with = concentration_with_bp519(&bp);
        // Adding one consistent member should not collapse the concentration.
        assert!(
            with > 0.25,
            "R̄ with BP519 = {with:.3} (was {:.3})",
            sample_varpi_concentration()
        );
    }
}

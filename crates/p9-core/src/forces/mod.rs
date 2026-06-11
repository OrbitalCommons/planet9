pub mod galactic_tide;
pub mod j2_secular;

use std::sync::Arc;

use nalgebra::Vector3;

/// Extra acceleration applied during the kick stage of an integrator.
///
/// These are position-dependent perturbations (potential kicks), so they
/// preserve the symplectic structure of the Wisdom-Holman splitting. They are
/// applied to massive bodies and test particles alike, in heliocentric
/// coordinates.
#[derive(Clone)]
pub enum ExtraForce {
    /// Orbit-averaged Jupiter+Saturn+Uranus: effective J2/J4 ring field plus
    /// the ΣGM_JSU monopole boost (see `j2_secular::combined_j2_jsu`).
    /// Valid for perihelia q ≳ 30 AU; integrate Neptune directly.
    J2Jsu,
    /// Milky Way vertical disk tide, evaluated in the galactic frame.
    GalacticTide,
    /// User-supplied acceleration field a(r) in AU/day², heliocentric.
    Custom(Arc<dyn Fn(&Vector3<f64>) -> Vector3<f64> + Send + Sync>),
}

impl ExtraForce {
    /// Acceleration at heliocentric position `pos` (AU), in AU/day².
    pub fn acceleration(&self, pos: &Vector3<f64>) -> Vector3<f64> {
        match self {
            ExtraForce::J2Jsu => {
                let (j2, j4, gm_boost) = j2_secular::combined_j2_jsu();
                j2_secular::secular_j2_acceleration(
                    pos,
                    crate::constants::GM_SUN,
                    j2,
                    j4,
                    gm_boost,
                    1.0,
                )
            }
            ExtraForce::GalacticTide => galactic_tide::galactic_tide_acceleration(pos),
            ExtraForce::Custom(f) => f(pos),
        }
    }
}

impl std::fmt::Debug for ExtraForce {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExtraForce::J2Jsu => write!(f, "ExtraForce::J2Jsu"),
            ExtraForce::GalacticTide => write!(f, "ExtraForce::GalacticTide"),
            ExtraForce::Custom(_) => write!(f, "ExtraForce::Custom(..)"),
        }
    }
}

/// Total extra acceleration from a list of extra forces.
pub fn total_extra_acceleration(forces: &[ExtraForce], pos: &Vector3<f64>) -> Vector3<f64> {
    let mut accel = Vector3::zeros();
    for force in forces {
        accel += force.acceleration(pos);
    }
    accel
}

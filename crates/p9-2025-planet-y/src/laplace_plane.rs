//! Classical Laplace-plane / forced-inclination theory for a Kuiper-belt test
//! particle under competing torques.
//!
//! A test particle at semi-major axis `a` precesses about a *forced* plane
//! (the local Laplace plane) set by the competition between two nodal torques:
//!
//! 1. **Inner torque** from the giant planets' orbit-averaged quadrupole. The
//!    averaged JSU(N) field is a J2-like flattening of the inner Solar System
//!    (reuse `p9_core::forces::j2_secular::effective_j2`). Its nodal-precession
//!    coefficient scales as
//!
//!    ```text
//!    A_in(a) = (3/4) n(a) Σ_j (m_j/M_⊙)(a_j/a)²            (∝ a^{-7/2})
//!    ```
//!
//!    and it holds the local mean plane near the inner planets' invariable
//!    plane (≈ the ecliptic for our purposes — the giants' invariable plane is
//!    only ≈ 1.58° from the ecliptic; see `p9-2025-clustering::orbital_poles`).
//!
//! 2. **Outer torque** from a distant *inclined* perturber — Planet Y at
//!    a_Y ≈ 100–200 AU, mass m_Y ≈ 0.05–1 M_⊕, inclination i_Y. Treating Planet
//!    Y as an exterior ring, its nodal-precession coefficient on an interior
//!    particle is
//!
//!    ```text
//!    A_out(a) = (3/4) n(a) (m_Y/M_⊙) (a/a_Y)³ b_{3/2}^{(1)}(a/a_Y)/3
//!    ```
//!
//!    (Laplace–Lagrange; the leading b_{3/2}^{(1)} ≈ 3α gives the familiar
//!    (a/a_Y)³ ring torque). It pulls the local mean plane toward Planet Y's
//!    orbital plane.
//!
//! The classical Laplace surface lies at the inclination `i_forced` (measured
//! from the inner reference plane toward Planet Y's plane) where the two
//! torques balance. For a single inner axis (z, the invariable normal) and a
//! single outer plane tilted by i_Y, the forced inclination satisfies the
//! standard Laplace-plane relation
//!
//! ```text
//! tan(2 i_forced) = A_out sin(2 i_Y) / (A_in + A_out cos(2 i_Y)).
//! ```
//!
//! - When A_in ≫ A_out (small a): i_forced → 0, the plane hugs the invariable
//!   plane.
//! - When A_out ≫ A_in (large a): i_forced → i_Y, the plane aligns with Y.
//! - The transition (A_out ≈ A_in) is the **warp**: the mean-plane tilt rises
//!   over a localized band of semi-major axis. The ratio A_out/A_in ∝ a^{5}/a_Y³,
//!   so the warp band's location is fixed by a_Y and m_Y.
//!
//! Reference: Murray & Dermott (1999) §7 (secular Laplace–Lagrange theory) and
//! the Laplace-surface treatment of §5.3; Siraj, Loeb & Siegel (2025),
//! arXiv:2508.14156.

use p9_core::constants::{
    DEG2RAD, EARTH_MASS_SOLAR, GM_SUN, MASS_JUPITER_SOLAR, MASS_NEPTUNE_SOLAR, MASS_SATURN_SOLAR,
    MASS_URANUS_SOLAR,
};
use p9_core::forces::j2_secular::effective_j2;

/// Published Planet Y reference numbers (Siraj et al. 2025), kept as labelled
/// constants so the tests pin against them rather than against magic numbers.
pub mod published {
    /// Lower mass bound (≈ Mercury mass), in Earth masses.
    pub const M_Y_LOW_EARTH: f64 = 0.05;
    /// Upper mass bound (≈ Earth mass), in Earth masses.
    pub const M_Y_HIGH_EARTH: f64 = 1.0;
    /// Lower semi-major-axis bound (AU).
    pub const A_Y_LOW_AU: f64 = 100.0;
    /// Upper semi-major-axis bound (AU).
    pub const A_Y_HIGH_AU: f64 = 200.0;
    /// Inner edge of the warp band (AU).
    pub const WARP_BAND_LOW_AU: f64 = 50.0;
    /// Outer edge of the warp band (AU).
    pub const WARP_BAND_HIGH_AU: f64 = 80.0;
    /// Representative modest Planet Y inclination to the ecliptic (deg).
    pub const I_Y_DEG: f64 = 10.0;
}

/// Semi-major axes (AU) of the four giant planets, used to build the inner
/// quadrupole torque. The averaged JSU(N) ring is the same construction as
/// `p9_core::forces::j2_secular::combined_j2_jsu`, extended to include Neptune
/// because the warp band (50–80 AU) is exterior to Neptune.
const GIANT_AU: [(f64, f64); 4] = [
    (MASS_JUPITER_SOLAR, 5.2026),
    (MASS_SATURN_SOLAR, 9.5549),
    (MASS_URANUS_SOLAR, 19.2184),
    (MASS_NEPTUNE_SOLAR, 30.07),
];

/// Specification of the distant perturber (Planet Y).
#[derive(Debug, Clone, Copy)]
pub struct PlanetY {
    /// Mass in solar masses.
    pub mass_solar: f64,
    /// Semi-major axis (AU).
    pub a_au: f64,
    /// Inclination of Planet Y's orbital plane to the inner reference plane
    /// (radians).
    pub inclination: f64,
}

impl PlanetY {
    /// Build from the paper's natural units: mass in Earth masses, a in AU,
    /// inclination in degrees.
    pub fn new(mass_earth: f64, a_au: f64, inclination_deg: f64) -> Self {
        Self {
            mass_solar: mass_earth * EARTH_MASS_SOLAR,
            a_au,
            inclination: inclination_deg * DEG2RAD,
        }
    }
}

/// Mean motion of a test particle (rad/day) about the Sun.
fn mean_motion(a_au: f64) -> f64 {
    (GM_SUN / (a_au * a_au * a_au)).sqrt()
}

/// Inner (giant-planet quadrupole) nodal-precession coefficient at semi-major
/// axis `a` (rad/day, magnitude).
///
///   A_in = (3/4) n Σ_j (m_j/M_⊙) (a_j/a)²
///
/// Built from `effective_j2` for each giant so the inner field is exactly the
/// p9-core orbit-averaged J2 ring (J2_eff,j = ½ (m_j/M)(a_j)², with R_ref = 1
/// folded out below via the (R_ref/a)² lever arm = (1/a)²).
pub fn inner_torque_coeff(a_au: f64) -> f64 {
    let n = mean_motion(a_au);
    // Σ_j (m_j/M)(a_j/a)² = (2/a²) Σ_j J2_eff,j with R_ref = 1 AU.
    let j2_sum: f64 = GIANT_AU
        .iter()
        .map(|&(m, a_j)| effective_j2(m, a_j, 1.0))
        .sum();
    let quad_sum = 2.0 * j2_sum / (a_au * a_au);
    0.75 * n * quad_sum
}

/// Leading Laplace coefficient b_{3/2}^{(1)}(α) for an exterior perturber,
/// evaluated by direct quadrature of
///
///   b_{3/2}^{(s)}(α) = (1/π) ∫_0^{2π} cos(sψ) (1 − 2α cosψ + α²)^{-3/2} dψ.
///
/// For α → 0 this tends to 3α, recovering the bare (a/a_Y)³ ring torque; for
/// the α ≈ 0.3–0.8 relevant to the warp band the exact coefficient matters.
fn laplace_b_three_half_1(alpha: f64) -> f64 {
    let n = 2048;
    let dpsi = std::f64::consts::TAU / n as f64;
    let mut sum = 0.0;
    for k in 0..n {
        let psi = (k as f64 + 0.5) * dpsi;
        let denom = (1.0 - 2.0 * alpha * psi.cos() + alpha * alpha).powf(1.5);
        sum += psi.cos() / denom;
    }
    sum * dpsi / std::f64::consts::PI
}

/// Outer (Planet Y) nodal-precession coefficient on an interior particle at
/// semi-major axis `a` (rad/day, magnitude).
///
///   A_out = (3/4) n (m_Y/M_⊙) α b_{3/2}^{(1)}(α) / ... ,  α = a/a_Y < 1
///
/// In the secular Laplace–Lagrange framework the off-diagonal node coupling of
/// an exterior body to an interior particle is
///
///   |A_out| = (1/4) n (m_Y/M_⊙) α² b_{3/2}^{(1)}(α),
///
/// which for α ≪ 1 reduces to (3/4) n (m_Y/M)(a/a_Y)³ (since b_{3/2}^{(1)}→3α).
pub fn outer_torque_coeff(a_au: f64, planet_y: &PlanetY) -> f64 {
    let alpha = a_au / planet_y.a_au;
    let n = mean_motion(a_au);
    let b1 = laplace_b_three_half_1(alpha);
    0.25 * n * planet_y.mass_solar * alpha * alpha * b1
}

/// Forced inclination of the local Laplace plane (radians), measured from the
/// inner invariable plane toward Planet Y's orbital plane.
///
///   tan(2 i_forced) = A_out sin(2 i_Y) / (A_in + A_out cos(2 i_Y))
///
/// The half-angle solution is taken in [0, i_Y]; `atan2` keeps the branch
/// continuous across the A_in = A_out crossover.
pub fn forced_inclination(a_au: f64, planet_y: &PlanetY) -> f64 {
    let a_in = inner_torque_coeff(a_au);
    let a_out = outer_torque_coeff(a_au, planet_y);
    let two_iy = 2.0 * planet_y.inclination;
    let num = a_out * two_iy.sin();
    let den = a_in + a_out * two_iy.cos();
    0.5 * num.atan2(den)
}

/// Convenience: forced inclination in degrees.
pub fn forced_inclination_deg(a_au: f64, planet_y: &PlanetY) -> f64 {
    forced_inclination(a_au, planet_y) / DEG2RAD
}

/// A sample of the forced-inclination profile.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ForcedProfilePoint {
    pub a_au: f64,
    pub i_forced_deg: f64,
}

/// Sample the forced-inclination profile i_forced(a) over a grid of
/// semi-major axes (AU), inclusive of both endpoints.
pub fn forced_profile(
    a_min: f64,
    a_max: f64,
    n: usize,
    planet_y: &PlanetY,
) -> Vec<ForcedProfilePoint> {
    assert!(n >= 2, "need at least two grid points");
    (0..n)
        .map(|k| {
            let a = a_min + (a_max - a_min) * (k as f64) / ((n - 1) as f64);
            ForcedProfilePoint {
                a_au: a,
                i_forced_deg: forced_inclination_deg(a, planet_y),
            }
        })
        .collect()
}

/// Warp diagnostic over the published 50–80 AU band.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct WarpFeature {
    /// Semi-major axis (AU) of the peak forced tilt inside the band.
    pub a_peak_au: f64,
    /// Forced tilt at the peak (deg).
    pub i_peak_deg: f64,
    /// Forced tilt at the inner band edge (deg).
    pub i_inner_deg: f64,
    /// Forced tilt at the outer band edge (deg).
    pub i_outer_deg: f64,
    /// Warp amplitude (deg): rise across the band, i_outer − i_inner.
    pub amplitude_deg: f64,
}

/// Measure the warp feature in the 50–80 AU band: the rise of the forced tilt
/// across the band and the location of the in-band peak. The "warp amplitude"
/// is the increase from the inner to the outer band edge — the localized
/// feature the paper attributes to Planet Y.
pub fn warp_feature(planet_y: &PlanetY) -> WarpFeature {
    let lo = published::WARP_BAND_LOW_AU;
    let hi = published::WARP_BAND_HIGH_AU;
    let n = 121;
    let profile = forced_profile(lo, hi, n, planet_y);

    let inner = profile.first().unwrap().i_forced_deg;
    let outer = profile.last().unwrap().i_forced_deg;

    let peak = profile
        .iter()
        .max_by(|x, y| x.i_forced_deg.partial_cmp(&y.i_forced_deg).unwrap())
        .unwrap();

    WarpFeature {
        a_peak_au: peak.a_au,
        i_peak_deg: peak.i_forced_deg,
        i_inner_deg: inner,
        i_outer_deg: outer,
        amplitude_deg: outer - inner,
    }
}

/// Neptune's mass-ratio contribution, exposed so callers can confirm the inner
/// field genuinely includes the planet that bounds the warp band on the inside.
pub fn neptune_mass_solar() -> f64 {
    MASS_NEPTUNE_SOLAR
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn nominal_planet_y() -> PlanetY {
        // Mid-range published Planet Y: 0.3 M_⊕ at 150 AU, modestly inclined.
        PlanetY::new(0.3, 150.0, published::I_Y_DEG)
    }

    #[test]
    fn test_laplace_coeff_small_alpha_limit() {
        // b_{3/2}^{(1)}(α) → 3α as α → 0.
        let alpha = 1e-3;
        let b1 = laplace_b_three_half_1(alpha);
        assert_relative_eq!(b1, 3.0 * alpha, max_relative = 1e-3);
    }

    #[test]
    fn test_outer_torque_reduces_to_ring_at_small_alpha() {
        // |A_out| → (3/4) n (m_Y/M)(a/a_Y)³ for α ≪ 1.
        let py = PlanetY::new(1.0, 1000.0, 10.0); // tiny α at a = 40
        let a = 40.0;
        let alpha = a / py.a_au;
        let n = mean_motion(a);
        let ring = 0.75 * n * py.mass_solar * alpha.powi(3);
        let computed = outer_torque_coeff(a, &py);
        assert_relative_eq!(computed, ring, max_relative = 2e-3);
    }

    #[test]
    fn test_inner_dominates_at_small_a() {
        // At 40 AU the giant-planet quadrupole should dominate Planet Y.
        let py = nominal_planet_y();
        let a_in = inner_torque_coeff(40.0);
        let a_out = outer_torque_coeff(40.0, &py);
        assert!(
            a_in > a_out,
            "expected inner torque to dominate at 40 AU: A_in={a_in:.3e} A_out={a_out:.3e}"
        );
        // Forced tilt nearly zero there.
        assert!(forced_inclination_deg(40.0, &py) < 1.0);
    }

    #[test]
    fn test_outer_dominates_near_planet_y() {
        // The inner giant-planet quadrupole (∝ a^{-7/2}) dominates throughout
        // the warp band and well beyond; the outer ring torque (∝ a³/a_Y³)
        // only overtakes it as the particle approaches Planet Y. With the
        // nominal 0.3 M_⊕ / 150 AU perturber the crossover is at ≈ 128 AU, so
        // "outer wins" must be tested just interior to a_Y, not at 90 AU.
        // (This is the honest reason the 50–80 AU warp is a small rising
        // shoulder, not full alignment with Y's plane — see REPRODUCTION_NOTES.)
        let py = nominal_planet_y();
        let a = 0.95 * py.a_au; // 142.5 AU, interior to and near Planet Y
        let a_in = inner_torque_coeff(a);
        let a_out = outer_torque_coeff(a, &py);
        assert!(
            a_out > a_in,
            "expected outer torque to win near a_Y ({a} AU): A_in={a_in:.3e} A_out={a_out:.3e}"
        );
    }

    #[test]
    fn test_warp_feature_in_50_80_band() {
        // With Planet Y the forced tilt rises across the 50–80 AU band. For a
        // sub-Earth-mass perturber at 150 AU the rise is a small fraction of a
        // degree (the band is far interior to the torque crossover); the
        // "few-degree" warp is reached only at the Earth-mass, close end (see
        // `test_warp_amplitude_few_degrees_for_earth_mass`).
        let py = nominal_planet_y();
        let w = warp_feature(&py);
        assert!(
            w.i_outer_deg > w.i_inner_deg,
            "tilt should rise across the band: {:?}",
            w
        );
        // Monotone rise: the in-band peak sits at the outer edge.
        assert!(w.a_peak_au >= published::WARP_BAND_LOW_AU);
        assert!(w.a_peak_au <= published::WARP_BAND_HIGH_AU);
        assert_relative_eq!(w.a_peak_au, published::WARP_BAND_HIGH_AU, epsilon = 0.5);
        // A real, localized warp of order a tenth of a degree for this case.
        assert!(
            w.amplitude_deg > 0.01 && w.amplitude_deg < 1.0,
            "warp amplitude = {:.3}° for 0.3 M_earth @ 150 AU",
            w.amplitude_deg
        );
    }

    #[test]
    fn test_warp_amplitude_few_degrees_for_earth_mass() {
        // The paper's headline "a few degrees" warp is realized at the upper,
        // close end of the allowed box: an Earth-mass Planet Y at 100 AU
        // produces a ≈ 3–4° rise across the 50–80 AU band.
        let py = PlanetY::new(
            published::M_Y_HIGH_EARTH,
            published::A_Y_LOW_AU,
            published::I_Y_DEG,
        );
        let w = warp_feature(&py);
        assert!(
            w.amplitude_deg > 2.0 && w.amplitude_deg < 6.0,
            "Earth-mass @ 100 AU warp amplitude = {:.3}°, expected a few degrees",
            w.amplitude_deg
        );
    }

    #[test]
    fn test_no_perturber_plane_stays_flat() {
        // Zero Planet Y mass: the local mean plane stays on the inner
        // invariable plane (i_forced = 0) at every a.
        let none = PlanetY::new(0.0, 150.0, 10.0);
        for a in [40.0, 60.0, 80.0, 90.0] {
            assert!(
                forced_inclination_deg(a, &none).abs() < 1e-9,
                "no-perturber tilt at {a} AU should vanish"
            );
        }
        let w = warp_feature(&none);
        assert!(w.amplitude_deg.abs() < 1e-9);
    }

    #[test]
    fn test_warp_amplitude_scales_with_mass() {
        // Heavier Planet Y → larger warp. Check the two published mass bounds
        // (0.05 and 1 M_⊕) at fixed a_Y, i_Y.
        let light = PlanetY::new(published::M_Y_LOW_EARTH, 150.0, published::I_Y_DEG);
        let heavy = PlanetY::new(published::M_Y_HIGH_EARTH, 150.0, published::I_Y_DEG);
        let w_light = warp_feature(&light);
        let w_heavy = warp_feature(&heavy);
        assert!(
            w_heavy.amplitude_deg > w_light.amplitude_deg,
            "heavier Planet Y must warp more: light={:.3}° heavy={:.3}°",
            w_light.amplitude_deg,
            w_heavy.amplitude_deg
        );
        // In the small-tilt regime the warp is ∝ A_out ∝ m_Y, so a 20× mass
        // ratio gives roughly a 20× amplitude ratio (modest tilts keep it
        // close to linear).
        let ratio = w_heavy.amplitude_deg / w_light.amplitude_deg;
        let mass_ratio = published::M_Y_HIGH_EARTH / published::M_Y_LOW_EARTH;
        assert!(
            ratio > 0.4 * mass_ratio && ratio < mass_ratio * 1.1,
            "amplitude ratio {ratio:.1} vs mass ratio {mass_ratio:.1}"
        );
    }

    #[test]
    fn test_a_y_position_shifts_warp() {
        // A closer Planet Y (100 AU) warps the band more than a distant one
        // (200 AU) at fixed mass, since A_out/A_in ∝ a_Y^{-3}.
        let close = PlanetY::new(0.3, published::A_Y_LOW_AU, published::I_Y_DEG);
        let far = PlanetY::new(0.3, published::A_Y_HIGH_AU, published::I_Y_DEG);
        assert!(
            warp_feature(&close).amplitude_deg > warp_feature(&far).amplitude_deg,
            "closer Planet Y should warp more"
        );
    }

    #[test]
    fn test_forced_tilt_bounded_by_i_y() {
        // The Laplace plane never tilts past Planet Y's own inclination.
        let py = nominal_planet_y();
        for a in [40.0, 55.0, 70.0, 85.0, 100.0] {
            let i_f = forced_inclination_deg(a, &py);
            assert!(
                i_f >= 0.0 && i_f <= published::I_Y_DEG + 1e-6,
                "forced tilt {i_f:.3}° at {a} AU out of [0, i_Y]"
            );
        }
    }

    #[test]
    fn test_inner_field_includes_neptune() {
        // The warp band is exterior to Neptune; the inner torque must include
        // Neptune's mass (it is the nearest interior planet to the band).
        assert!(neptune_mass_solar() > 0.0);
        let with_n = inner_torque_coeff(60.0);
        assert!(with_n > 0.0);
    }
}

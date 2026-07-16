//! Orbit-averaged (secular) dynamics of an ETNO under the EFE quadrupole tide.
//!
//! The EFE adds the axisymmetric quadrupole potential (see [`crate::efe`])
//!
//!   Φ_efe(r) = -(A/2) [ (r·n̂)² - r²/3 ],
//!
//! about the Sun→galactic-center axis `n̂`. Averaging over one Kepler orbit at
//! fixed `a` removes the fast mean-anomaly dependence and leaves a secular
//! disturbing function `R̄(e, i, ω, Ω)` that drives slow apsidal and nodal
//! evolution. Brown & Mathur's headline is the *forced apsidal equilibrium*:
//! the orbit-averaged torque pushes the longitude of perihelion `ϖ` toward the
//! in-plane projection of `n̂` (the galactic-center direction), reproducing the
//! observed clustering with no planet.
//!
//! We use the *coplanar* idealization (the ETNOs share, to first order, a
//! common orbital plane; `n̂` is projected into it). The orbit average is done
//! by Kepler quadrature in eccentric anomaly — exact for this smooth integrand
//! and mirroring the Gauss-ring averaging style of `p9_core::analysis::secular`.
//!
//! Lagrange's planetary equation for the apse, with `R̄` the orbit-averaged
//! disturbing function (per unit mass) and `n = sqrt(GM/a³)` the mean motion:
//!
//!   dϖ/dt = (√(1−e²) / (n a² e)) ∂R̄/∂e            (coplanar; ω̇ + Ω̇)
//!
//! Reference: Brown & Mathur (2023), arXiv:2304.00576; companion Migaszewski
//! (2023), arXiv:2303.13339.

use std::f64::consts::{PI, TAU};

use nalgebra::Vector3;

use p9_core::analysis::elements::mean_motion;
use p9_core::units::{days, radians, AngularVelocity};

use crate::efe::efe_potential;

/// In-plane geometry of the secular problem for a coplanar orbit whose plane
/// has unit normal `z_hat` and whose perihelion-reference axis is `x_hat`.
pub struct PlaneGeometry {
    /// In-plane unit vector along the reference (ϖ = 0) direction.
    pub x_hat: Vector3<f64>,
    /// In-plane unit vector 90° ahead of `x_hat` (right-handed with `z_hat`).
    pub y_hat: Vector3<f64>,
    /// Orbit-plane normal.
    pub z_hat: Vector3<f64>,
}

impl PlaneGeometry {
    /// Build an orbital plane from its normal. The reference axis is chosen as
    /// the normalized projection of the ecliptic x-axis into the plane (an
    /// arbitrary but fixed ϖ = 0 zero point); `y_hat = z_hat × x_hat`.
    pub fn from_normal(z_hat: Vector3<f64>) -> Self {
        let z_hat = z_hat.normalize();
        let mut seed = Vector3::new(1.0, 0.0, 0.0);
        if seed.cross(&z_hat).norm() < 1e-6 {
            seed = Vector3::new(0.0, 1.0, 0.0);
        }
        let x_hat = (seed - seed.dot(&z_hat) * z_hat).normalize();
        let y_hat = z_hat.cross(&x_hat);
        Self {
            x_hat,
            y_hat,
            z_hat,
        }
    }

    /// Longitude (radians, [0, 2π)) of an arbitrary direction's in-plane
    /// projection, measured from `x_hat` toward `y_hat` (the ϖ convention).
    pub fn in_plane_longitude(&self, v: &Vector3<f64>) -> f64 {
        v.dot(&self.y_hat).atan2(v.dot(&self.x_hat)).rem_euclid(TAU)
    }
}

/// Orbit-averaged EFE disturbing function R̄ = -⟨Φ_efe⟩ (AU²/day²) for a
/// coplanar orbit of semi-major axis `a`, eccentricity `e`, and longitude of
/// perihelion `varpi` (measured from `geom.x_hat`), under EFE amplitude
/// `a_efe` and symmetry axis `n`.
///
/// Averaged over eccentric anomaly with the dM = (1 − e cos E) dE Jacobian.
pub fn averaged_disturbing(
    a: f64,
    e: f64,
    varpi: f64,
    a_efe: f64,
    n: &Vector3<f64>,
    geom: &PlaneGeometry,
    n_quad: usize,
) -> f64 {
    // Perihelion direction in 3D.
    let (sv, cv) = varpi.sin_cos();
    let p_hat = cv * geom.x_hat + sv * geom.y_hat; // toward perihelion
    let q_hat = -sv * geom.x_hat + cv * geom.y_hat; // in-plane, 90° ahead

    let de = TAU / n_quad as f64;
    let mut sum = 0.0;
    let mut wsum = 0.0;
    for k in 0..n_quad {
        let ea = (k as f64 + 0.5) * de;
        // Position in the orbit plane: x' = a(cosE − e), y' = a√(1−e²) sinE,
        // with x' along perihelion.
        let xp = a * (ea.cos() - e);
        let yp = a * (1.0 - e * e).sqrt() * ea.sin();
        let r = xp * p_hat + yp * q_hat;
        let w = 1.0 - e * ea.cos(); // dM/dE
        sum += w * efe_potential(a_efe, n, &r);
        wsum += w;
    }
    -sum / wsum
}

/// Mean motion n = sqrt(GM_sun / a³) as a typed [`AngularVelocity`] (rad per
/// day).
pub fn mean_motion_typed(a: f64) -> AngularVelocity {
    (radians(mean_motion(a)) / days(1.0)).into()
}

/// Secular apsidal precession rate dϖ/dt (rad/day) from Lagrange's planetary
/// equation, evaluated numerically:
///
///   dϖ/dt = (√(1−e²) / (n a² e)) ∂R̄/∂e.
///
/// (Coplanar: ϖ̇ = ω̇ + Ω̇; the nodal/inclination terms vanish for the in-plane
/// projection used here.)
pub fn precession_rate(
    a: f64,
    e: f64,
    varpi: f64,
    a_efe: f64,
    n: &Vector3<f64>,
    geom: &PlaneGeometry,
    n_quad: usize,
) -> f64 {
    let he = 1e-5;
    let r_plus = averaged_disturbing(a, e + he, varpi, a_efe, n, geom, n_quad);
    let r_minus = averaged_disturbing(a, e - he, varpi, a_efe, n, geom, n_quad);
    let dr_de = (r_plus - r_minus) / (2.0 * he);
    let nm = (mean_motion_typed(a) / (radians(1.0) / days(1.0))).value;
    (1.0 - e * e).sqrt() / (nm * a * a * e) * dr_de
}

/// ∂R̄/∂ϖ (the apsidal torque, AU²/day²/rad), evaluated numerically. The
/// forced equilibrium is where this vanishes.
pub fn apsidal_torque(
    a: f64,
    e: f64,
    varpi: f64,
    a_efe: f64,
    n: &Vector3<f64>,
    geom: &PlaneGeometry,
    n_quad: usize,
) -> f64 {
    let hv = 1e-5;
    let r_plus = averaged_disturbing(a, e, varpi + hv, a_efe, n, geom, n_quad);
    let r_minus = averaged_disturbing(a, e, varpi - hv, a_efe, n, geom, n_quad);
    (r_plus - r_minus) / (2.0 * hv)
}

/// Locate the *forced* longitude of perihelion ϖ_forced: the value of ϖ that
/// extremizes R̄ (zero apsidal torque) and is *stable* (a minimum of the
/// secular potential energy Φ̄ = −R̄, equivalently a maximum of R̄, i.e. the
/// libration center). Returned in [0, 2π).
///
/// Found by scanning ϖ on a coarse grid for the global maximum of R̄, then
/// refining with a few Newton steps on the torque.
pub fn forced_varpi(
    a: f64,
    e: f64,
    a_efe: f64,
    n: &Vector3<f64>,
    geom: &PlaneGeometry,
    n_quad: usize,
) -> f64 {
    // Coarse scan for the global maximum of R̄ over ϖ.
    let n_scan = 720;
    let mut best_v = 0.0;
    let mut best_r = f64::NEG_INFINITY;
    for k in 0..n_scan {
        let v = k as f64 / n_scan as f64 * TAU;
        let r = averaged_disturbing(a, e, v, a_efe, n, geom, n_quad);
        if r > best_r {
            best_r = r;
            best_v = v;
        }
    }
    // Newton refinement on ∂R̄/∂ϖ = 0.
    let mut v = best_v;
    for _ in 0..40 {
        let t = apsidal_torque(a, e, v, a_efe, n, geom, n_quad);
        let hv = 1e-4;
        let tp = apsidal_torque(a, e, v + hv, a_efe, n, geom, n_quad);
        let tm = apsidal_torque(a, e, v - hv, a_efe, n, geom, n_quad);
        let d2 = (tp - tm) / (2.0 * hv);
        if d2.abs() < 1e-300 {
            break;
        }
        let step = t / d2;
        v -= step;
        if step.abs() < 1e-10 {
            break;
        }
    }
    v.rem_euclid(TAU)
}

/// Small-amplitude libration frequency about the forced equilibrium ϖ_forced
/// (rad/day), from the linearized (e, ϖ) secular system
///
///   δϖ̇ = A·δe,   δė = −B·δϖ,   ω_lib² = A·B,
///
/// with A = ∂ϖ̇/∂e and B = C(e)·∂²R̄/∂ϖ² (from the Lagrange pair
/// ϖ̇ = C·∂R̄/∂e, ė = −C·∂R̄/∂ϖ, C = √(1−e²)/(n a² e)). The diagonal terms
/// vanish at the fixed point — for the pure-second-harmonic EFE tide
/// ∂ϖ̇/∂ϖ|_eq ≡ 0 identically, so a previous version that finite-differenced
/// exactly that quantity returned √(quadrature noise), and its "stability"
/// test passed on noise. Returns √|ω²| with stability decided by the SIGN of
/// A·B (positive ⇒ stable libration center).
pub fn libration_frequency(
    a: f64,
    e: f64,
    a_efe: f64,
    n: &Vector3<f64>,
    geom: &PlaneGeometry,
    n_quad: usize,
) -> f64 {
    libration_frequency_squared(a, e, a_efe, n, geom, n_quad)
        .abs()
        .sqrt()
}

/// Signed ω_lib² = (∂ϖ̇/∂e)·(C·∂²R̄/∂ϖ²) at the forced apse (rad²/day²);
/// positive for a stable libration center.
pub fn libration_frequency_squared(
    a: f64,
    e: f64,
    a_efe: f64,
    n: &Vector3<f64>,
    geom: &PlaneGeometry,
    n_quad: usize,
) -> f64 {
    let varpi0 = forced_varpi(a, e, a_efe, n, geom, n_quad);

    // A = ∂ϖ̇/∂e at the fixed point.
    let he = 1e-4;
    let rate_p = precession_rate(a, e + he, varpi0, a_efe, n, geom, n_quad);
    let rate_m = precession_rate(a, e - he, varpi0, a_efe, n, geom, n_quad);
    let a_coef = (rate_p - rate_m) / (2.0 * he);

    // B = C(e)·∂²R̄/∂ϖ² via the torque derivative: ė = −C·∂R̄/∂ϖ, so
    // ∂ė/∂ϖ = −C·∂²R̄/∂ϖ².
    let hv = 1e-3;
    let t_p = apsidal_torque(a, e, varpi0 + hv, a_efe, n, geom, n_quad);
    let t_m = apsidal_torque(a, e, varpi0 - hv, a_efe, n, geom, n_quad);
    let d2r_dv2 = (t_p - t_m) / (2.0 * hv);
    let nm = (mean_motion_typed(a) / (radians(1.0) / days(1.0))).value;
    let c = (1.0 - e * e).sqrt() / (nm * a * a * e);
    let b_coef = c * d2r_dv2;

    // δϖ̈ = A·δė = −A·B·δϖ with B = C·R̄_ϖϖ: ω² = A·C·R̄_ϖϖ. Stability sign:
    // restoring when A and R̄_ϖϖ combine to a positive squared frequency.
    -a_coef * (-b_coef)
}

/// Smallest signed angular separation between two longitudes (radians), in
/// (−π, π].
pub fn angle_diff(a: f64, b: f64) -> f64 {
    let mut d = (a - b).rem_euclid(TAU);
    if d > PI {
        d -= TAU;
    }
    d
}

/// Separation between two *apsidal lines* (radians, in [0, π/2]). The EFE
/// quadrupole is symmetric under ϖ → ϖ + π, so the forced state is an aligned
/// *line of apsides*, not a single direction; this folds the longitude
/// difference into [0, π/2].
pub fn apsidal_line_separation(a: f64, b: f64) -> f64 {
    let d = angle_diff(a, b).abs();
    if d > PI / 2.0 {
        PI - d
    } else {
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::efe::galactic_center_ecliptic;
    use approx::assert_relative_eq;
    use p9_core::constants::GM_SUN;

    /// A representative ETNO orbital plane. We use the ecliptic plane (small
    /// ETNO inclinations) so the GC projection's ecliptic longitude is the
    /// directly comparable quantity.
    fn ecliptic_plane() -> PlaneGeometry {
        PlaneGeometry::from_normal(Vector3::new(0.0, 0.0, 1.0))
    }

    #[test]
    fn typed_mean_motion_matches_formula() {
        use uom::si::angular_velocity::radian_per_second;
        let a = 500.0;
        // rad/day read back as rad/s matches n = sqrt(GM/a³) / seconds-per-day.
        let n_rad_per_day = (GM_SUN / (a * a * a)).sqrt();
        assert_relative_eq!(
            mean_motion_typed(a).get::<radian_per_second>(),
            n_rad_per_day / 86_400.0,
            max_relative = 1e-9
        );
    }

    #[test]
    fn test_forced_perihelion_aligns_with_galactic_center_longitude() {
        // The headline of Brown & Mathur (2023): the EFE forces ϖ toward the
        // galactic-center direction's ecliptic longitude.
        let n = galactic_center_ecliptic();
        let geom = ecliptic_plane();
        let a = 350.0;
        let e = 0.7;
        let a_efe = 5.0e-13;

        let varpi_forced = forced_varpi(a, e, a_efe, &n, &geom, 256);
        let gc_lon = geom.in_plane_longitude(&n); // GC projection longitude

        // Both expressed as ecliptic longitudes (the plane IS the ecliptic).
        // The quadrupole aligns the *line of apsides* with n̂ (ϖ→ϖ+π
        // symmetry), so compare apsidal lines.
        let varpi_deg = varpi_forced.to_degrees();
        let gc_deg = gc_lon.to_degrees();
        let sep_deg = apsidal_line_separation(varpi_forced, gc_lon).to_degrees();

        assert!(
            sep_deg < 3.0,
            "ϖ_forced = {varpi_deg:.2}° vs GC longitude {gc_deg:.2}° (apsidal sep {sep_deg:.2}°)"
        );

        // And the GC ecliptic longitude itself is ~267° (Sgr A* direction).
        assert!(
            (gc_deg - 266.8).abs() < 1.5,
            "GC ecliptic longitude {gc_deg:.2}°"
        );
    }

    #[test]
    fn test_forced_equilibrium_is_stable_libration_center() {
        let n = galactic_center_ecliptic();
        let geom = ecliptic_plane();
        let a = 350.0;
        let e = 0.7;
        let a_efe = 5.0e-13;

        // Zero torque at the equilibrium, relative to the torque scale ~A·a²
        // (the peak |∂R̄/∂ϖ| is of order A·a²·e²; require ≤ 1e-8 of A·a²).
        let v0 = forced_varpi(a, e, a_efe, &n, &geom, 256);
        let torque = apsidal_torque(a, e, v0, a_efe, &n, &geom, 256);
        let torque_scale = a_efe * a * a;
        assert!(
            torque.abs() < 1e-8 * torque_scale,
            "residual torque {torque:e} vs scale {torque_scale:e}"
        );

        // Real libration frequency ⇒ genuine restoring (stable) center.
        let omega_lib = libration_frequency(a, e, a_efe, &n, &geom, 256);
        assert!(omega_lib > 0.0 && omega_lib.is_finite());

        // The equilibrium 90° away is an unstable maximum of Φ̄ (zero torque
        // there too, but R̄ is smaller): confirm R̄ is larger at v0.
        let r_eq = averaged_disturbing(a, e, v0, a_efe, &n, &geom, 256);
        let r_perp = averaged_disturbing(a, e, v0 + PI / 2.0, a_efe, &n, &geom, 256);
        assert!(
            r_eq > r_perp,
            "v0 should maximize R̄: {r_eq:e} vs {r_perp:e}"
        );
    }

    #[test]
    fn test_precession_rate_scales_linearly_with_efe_amplitude() {
        // dϖ/dt ∝ A_efe: double the amplitude, double the rate.
        let n = galactic_center_ecliptic();
        let geom = ecliptic_plane();
        let a = 350.0;
        let e = 0.7;
        // Evaluate the rate at a fixed ϖ offset from the forced point (so e is
        // not stationary and the rate is non-trivial).
        let v0 = forced_varpi(a, e, 1.0e-12, &n, &geom, 256);
        let varpi = v0 + 0.6;

        let a1 = 3.0e-13;
        let a2 = 6.0e-13;
        let rate1 = precession_rate(a, e, varpi, a1, &n, &geom, 256);
        let rate2 = precession_rate(a, e, varpi, a2, &n, &geom, 256);

        assert!(rate1.abs() > 0.0);
        let ratio = rate2 / rate1;
        assert_relative_eq!(ratio, 2.0, epsilon = 1e-6);
    }

    #[test]
    fn test_libration_frequency_scales_linearly_with_amplitude() {
        // In this ISOLATED-tide secular problem both linearized coefficients
        // (A = ∂ϖ̇/∂e and B = C·∂²R̄/∂ϖ²) are linear in A_efe, so
        // ω_lib = √(A·B) ∝ A_efe: quadrupling the amplitude quadruples the
        // frequency (timescale ∝ 1/A_efe). The old test asserted ω ∝ √A —
        // an artifact of the noise-based ω it was testing (√ of a quantity
        // whose numerical residue happened to scale with A). A √A law would
        // hold only with an A-independent restoring term (e.g. the giants'
        // precession), which this crate deliberately excludes.
        let n = galactic_center_ecliptic();
        let geom = ecliptic_plane();
        let a = 350.0;
        let e = 0.7;
        let w1 = libration_frequency(a, e, 2.5e-13, &n, &geom, 256);
        let w4 = libration_frequency(a, e, 1.0e-12, &n, &geom, 256);
        assert_relative_eq!(w4 / w1, 4.0, epsilon = 1e-3);
        // And the equilibrium is a genuinely STABLE center: signed ω² > 0.
        let w2 = libration_frequency_squared(a, e, 5.0e-13, &n, &geom, 256);
        assert!(w2 > 0.0, "omega^2 = {w2:.3e} must be positive (stable)");
    }

    #[test]
    fn test_forced_direction_independent_of_eccentricity() {
        // The forced apsidal direction is set by n̂'s geometry, not e: the
        // alignment toward the GC longitude holds across the ETNO e range.
        let n = galactic_center_ecliptic();
        let geom = ecliptic_plane();
        let a = 350.0;
        let a_efe = 5.0e-13;
        let gc_lon = geom.in_plane_longitude(&n);
        for &e in &[0.5, 0.7, 0.85] {
            let v = forced_varpi(a, e, a_efe, &n, &geom, 256);
            let sep = apsidal_line_separation(v, gc_lon).to_degrees();
            assert!(sep < 3.0, "e={e}: sep {sep:.2}°");
        }
    }

    #[test]
    fn test_averaged_potential_recovers_monopole_free_quadrupole() {
        // ⟨(r·n̂)² − r²/3⟩ for a circular orbit (e=0) in a plane: with n̂ in
        // the plane, the average of (r·n̂)² over a circle of radius a is a²/2,
        // and ⟨r²⟩ = a², so ⟨Φ⟩ = -(A/2)(a²/2 − a²/3) = -(A/2)(a²/6).
        let geom = ecliptic_plane();
        let n = Vector3::new(1.0, 0.0, 0.0); // in-plane axis
        let a = 100.0;
        let a_efe = 1.0e-12;
        let rbar = averaged_disturbing(a, 0.0, 0.3, a_efe, &n, &geom, 512);
        let expected = 0.5 * a_efe * (a * a / 6.0); // R̄ = -⟨Φ⟩
        assert_relative_eq!(rbar, expected, epsilon = 1e-6 * expected.abs());
    }
}

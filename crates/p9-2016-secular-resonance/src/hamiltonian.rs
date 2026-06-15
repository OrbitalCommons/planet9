//! The coplanar 1-degree-of-freedom secular Hamiltonian for a distant test
//! particle perturbed by an eccentric Planet Nine.
//!
//! At fixed semi-major axis `a` the doubly-averaged dynamics has a single
//! degree of freedom. The canonical pair is the apsidal angle and its
//! conjugate Delaunay action:
//!
//! ```text
//!     angle:    Δϖ = ϖ − ϖ₉            (the resonant angle)
//!     action:   Γ  = √(GM·a)·(1 − √(1 − e²))
//! ```
//!
//! `Γ` is a monotone increasing function of `e` at fixed `a` (the "deficit"
//! action; it is the quantity conserved together with `a` when the perturber
//! is axisymmetric). Hamilton's equations are
//!
//! ```text
//!     dΔϖ/dt = ∂H/∂Γ,      dΓ/dt = −∂H/∂Δϖ.
//! ```
//!
//! The Hamiltonian itself is the exact doubly-averaged interaction supplied by
//! [`p9_core::analysis::secular::numerical_secular_hamiltonian`] (Gauss-ring
//! quadrature of 1/Δ; no α expansion). Because the system is autonomous and
//! one-dimensional, level curves `H(Δϖ, Γ) = const` *are* the secular
//! trajectories: closed curves encircling the anti-aligned fixed point are
//! libration; curves spanning all Δϖ are circulation. The libration module
//! integrates the equations of motion to confirm this directly.

use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::{GM_SUN, TWO_PI};
use p9_core::types::P9Params;

/// A coplanar secular model: a fixed test-particle semi-major axis under a
/// fixed (eccentric) Planet Nine. The perturber defines the reference frame,
/// so the particle's apsidal angle relative to it is Δϖ and `i = 0`,
/// `Ω = 0` (coplanar, anti-/aligned dynamics live entirely in Δϖ).
#[derive(Debug, Clone, Copy)]
pub struct SecularModel {
    /// Test-particle semi-major axis (AU).
    pub a: f64,
    /// Planet Nine semi-major axis (AU).
    pub a_p: f64,
    /// Planet Nine eccentricity.
    pub e_p: f64,
    /// Planet Nine GM (AU³/day²).
    pub gm_p: f64,
    /// Quadrature nodes per anomaly for the ring average.
    pub n_quad: usize,
    /// Softening length (AU) regularizing orbit-crossing geometries.
    pub softening: f64,
}

impl SecularModel {
    /// Build from a test-particle `a` and a [`P9Params`] perturber.
    pub fn new(a: f64, p9: &P9Params) -> Self {
        Self {
            a,
            a_p: p9.a,
            e_p: p9.e,
            gm_p: p9.gm(),
            n_quad: 64,
            softening: 0.02 * p9.a,
        }
    }

    /// Override the perturber GM (used to scan Planet Nine's mass).
    pub fn with_gm_p(mut self, gm_p: f64) -> Self {
        self.gm_p = gm_p;
        self
    }

    /// The conserved secular Hamiltonian `H(e, Δϖ)` (AU²/day², up to the
    /// constant monopole term that does not affect the dynamics).
    pub fn hamiltonian(&self, e: f64, dvarpi: f64) -> f64 {
        // Coplanar: i = 0, Ω = 0, ω = Δϖ measured from the perturber apse.
        numerical_secular_hamiltonian(
            self.a,
            e,
            0.0,
            dvarpi,
            0.0,
            self.a_p,
            self.e_p,
            self.gm_p,
            self.n_quad,
            self.softening,
        )
    }

    /// Delaunay deficit action Γ = √(GM·a)·(1 − √(1 − e²)) (AU²/day),
    /// monotone in `e` at fixed `a`.
    pub fn action(&self, e: f64) -> f64 {
        (GM_SUN * self.a).sqrt() * (1.0 - (1.0 - e * e).sqrt())
    }

    /// Invert [`action`](Self::action): recover `e` from Γ at fixed `a`.
    pub fn eccentricity(&self, gamma: f64) -> f64 {
        let l = (GM_SUN * self.a).sqrt();
        let root = 1.0 - gamma / l; // = √(1 − e²)
        let r2 = root.clamp(0.0, 1.0).powi(2);
        (1.0 - r2).max(0.0).sqrt()
    }

    /// ∂H/∂Δϖ by a centered finite difference (drives dΓ/dt).
    pub fn dh_ddvarpi(&self, e: f64, dvarpi: f64) -> f64 {
        let h = 1e-4;
        (self.hamiltonian(e, dvarpi + h) - self.hamiltonian(e, dvarpi - h)) / (2.0 * h)
    }

    /// ∂H/∂Γ by a centered finite difference in the action (drives dΔϖ/dt).
    /// Differencing in Γ keeps the equations canonical.
    pub fn dh_dgamma(&self, e: f64, dvarpi: f64) -> f64 {
        let gamma = self.action(e);
        let dg = 1e-5 * (GM_SUN * self.a).sqrt();
        let e_plus = self.eccentricity(gamma + dg);
        let e_minus = self.eccentricity(gamma - dg);
        (self.hamiltonian(e_plus, dvarpi) - self.hamiltonian(e_minus, dvarpi)) / (2.0 * dg)
    }

    /// Apsidal Fourier decomposition of the secular Hamiltonian at fixed `e`:
    ///
    /// ```text
    ///     H(e, Δϖ) = a₀(e) + c₁(e)·cos Δϖ + c₂(e)·cos 2Δϖ + …
    /// ```
    ///
    /// Returns `(a₀, c₁, c₂)`. The axisymmetric mean `a₀` is the bulk of `H`
    /// (it sets the free precession); the apsidal harmonics `c₁` (octupole-like,
    /// the dominant aligned/anti-aligned forcing) and `c₂` are the small but
    /// cleanly convergent terms that build the libration island. Sine terms
    /// vanish by the perturber's apsidal-line reflection symmetry. Computed by
    /// a 64-point DFT over Δϖ — far better conditioned than finite-differencing
    /// the angle directly, since each harmonic is integrated coherently.
    pub fn apsidal_harmonics(&self, e: f64) -> (f64, f64, f64) {
        let n = 48usize;
        let mut a0 = 0.0;
        let mut c1 = 0.0;
        let mut c2 = 0.0;
        for k in 0..n {
            let dv: f64 = TWO_PI * k as f64 / n as f64;
            let h = self.hamiltonian(e, dv);
            a0 += h;
            c1 += h * dv.cos();
            c2 += h * (2.0 * dv).cos();
        }
        let nf = n as f64;
        (a0 / nf, 2.0 * c1 / nf, 2.0 * c2 / nf)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use approx::assert_relative_eq;

    #[test]
    fn action_is_monotone_and_invertible() {
        let m = SecularModel::new(250.0, &published::nominal_p9());
        let mut last = -1.0;
        for k in 0..20 {
            let e = 0.02 + 0.9 * k as f64 / 19.0;
            let g = m.action(e);
            assert!(g > last, "action not monotone at e = {e}");
            last = g;
            // Round-trip.
            assert_relative_eq!(m.eccentricity(g), e, epsilon = 1e-9);
        }
    }

    #[test]
    fn hamiltonian_finite_and_bound() {
        let m = SecularModel::new(250.0, &published::nominal_p9());
        let h = m.hamiltonian(0.7, std::f64::consts::PI);
        assert!(h.is_finite() && h < 0.0, "H = {h}");
    }

    #[test]
    fn anti_aligned_is_a_stationary_point_in_dvarpi() {
        // By the perturber's reflection symmetry (the eccentric ring is
        // symmetric about its apsidal line), Δϖ = π must be a stationary
        // point of H in the angle: ∂H/∂Δϖ = 0 there.
        let m = SecularModel::new(250.0, &published::nominal_p9());
        let slope = m.dh_ddvarpi(0.7, std::f64::consts::PI);
        // Scale: compare against a generic non-symmetric slope.
        let generic = m.dh_ddvarpi(0.7, 1.0).abs();
        assert!(
            slope.abs() < 0.05 * generic.max(1e-30),
            "Δϖ = π not stationary: slope {slope:.3e} vs generic {generic:.3e}"
        );
    }
}

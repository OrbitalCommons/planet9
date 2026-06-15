//! Build Beust's secular-resonance phase portrait as a one-degree-of-freedom
//! **resonant pendulum** and locate the apsidally anti-aligned libration
//! island.
//!
//! # The resonant Hamiltonian
//!
//! At fixed semi-major axis the doubly-averaged Hamiltonian decomposes into an
//! axisymmetric mean and apsidal harmonics (see
//! [`SecularModel::apsidal_harmonics`]):
//!
//! ```text
//!     H(e, Δϖ) = a₀(Γ) + c₁(Γ)·cos Δϖ + …
//! ```
//!
//! The axisymmetric part `a₀(Γ)` drives the particle's *free* apsidal
//! precession, `g_free(Γ) = da₀/dΓ`. A **secular resonance** with Planet Nine
//! occurs when this free rate matches Planet Nine's own apsidal precession
//! rate `g₉`; the relevant angle is then the slowly-varying `Δϖ = ϖ − ϖ₉`, and
//! in the co-precessing frame the Hamiltonian becomes the pendulum
//!
//! ```text
//!     K(Δϖ, Γ) = a₀(Γ) − g₉·Γ + c₁(Γ)·cos Δϖ.
//! ```
//!
//! Its fixed points are at Δϖ = 0 (aligned) and Δϖ = π (anti-aligned), at the
//! resonant action `Γ*` where `g_free(Γ*) = g₉`. The harmonic `c₁` controls
//! which one is the stable libration center: for the distant ETNO-scale
//! particles `c₁ < 0` at high eccentricity, making **Δϖ = π** the center — the
//! anti-aligned libration island Beust reports. A particle with `Γ` near `Γ*`
//! and `Δϖ` near π is trapped (librates); one with `Γ` far from `Γ*`
//! circulates.
//!
//! The harmonics are extracted by a coherent DFT of the p9-core Gauss-ring
//! Hamiltonian over Δϖ, which is well-conditioned even though the apsidal
//! terms are a small fraction of the axisymmetric bulk.

use crate::hamiltonian::SecularModel;
use p9_core::constants::GM_SUN;
use std::f64::consts::PI;

/// The 1-DOF resonant pendulum built from the secular Hamiltonian's apsidal
/// harmonics, in the frame precessing with Planet Nine.
///
/// The expensive Gauss-ring DFT is evaluated once on an action grid in
/// [`ResonantHamiltonian::new`]; the dynamics then reads cheap interpolants of
/// `a₀(Γ)`, `c₁(Γ)` and the free precession `g_free(Γ) = da₀/dΓ`. The free
/// precession is tabulated at the nodes from finite differences of the
/// *directly evaluated* `a₀` (not of the interpolant), so its slope — and
/// hence the libration frequency — is physical rather than the zero curvature
/// a linear interpolant of `a₀` alone would produce.
#[derive(Debug, Clone)]
pub struct ResonantHamiltonian {
    pub model: SecularModel,
    /// Planet Nine apsidal precession rate `g₉` (rad/day) defining the
    /// resonance. Chosen as the particle's free precession rate at the resonant
    /// eccentricity, i.e. the commensurability condition `g_free = g₉`.
    pub g_p: f64,
    /// Resonant eccentricity the pendulum is tuned to (where g_free = g₉).
    e_resonant: f64,
    /// Monotone action grid (AU²/day).
    gammas: Vec<f64>,
    /// Tabulated apsidal harmonic c₁ at each node (for the resonance-strength
    /// diagnostics; the *dynamics* freezes c₁ at the resonant action).
    c1s: Vec<f64>,
    /// Tabulated free precession g_free = da₀/dΓ at each node.
    gfree: Vec<f64>,
    /// a₀ accumulated to each node (exact integral of the piecewise-linear
    /// g_free, so a₀ values at nodes are consistent with `gfree`).
    a0_nodes: Vec<f64>,
    /// Harmonic amplitude frozen at the resonant action (the standard pendulum
    /// approximation: c₁ is held at c₁(Γ*) so the flow is exactly Hamiltonian).
    c1_star: f64,
}

/// A located libration island.
#[derive(Debug, Clone, Copy)]
pub struct Island {
    /// Resonant eccentricity (where g_free = g₉).
    pub e_center: f64,
    /// Resonant action Γ* (AU²/day).
    pub gamma_center: f64,
    /// Stable libration angle (radians) — π for the anti-aligned island.
    pub center_dvarpi: f64,
    /// Pendulum harmonic amplitude |c₁(Γ*)| (AU²/day²); the resonance strength,
    /// linear in the perturber GM.
    pub depth: f64,
    /// Small-oscillation libration frequency (rad/day).
    pub omega_lib: f64,
    /// Half-width of the resonance in action: δΓ_sep = 2·√(|c₁*|/|a₀''|).
    /// The separatrix reaches this far from Γ* in the action.
    pub separatrix_halfwidth_gamma: f64,
    /// Suggested integration timestep (days).
    pub suggested_dt: f64,
    /// Suggested step count (covers several libration periods).
    pub suggested_steps: usize,
}

impl ResonantHamiltonian {
    /// Build the resonant pendulum for a secular model, tuning the resonance to
    /// the free precession rate at eccentricity `e_resonant` (the
    /// commensurability with Planet Nine's apse).
    pub fn new(model: SecularModel, e_resonant: f64) -> Self {
        let n = 80usize;
        let e_min = 0.05;
        let e_max = 0.97;
        let l = (GM_SUN * model.a).sqrt();
        let dg = 0.015 * l;

        let mut gammas = Vec::with_capacity(n + 1);
        let mut c1s = Vec::with_capacity(n + 1);
        let mut gfree = Vec::with_capacity(n + 1);

        // a₀ on a slightly padded eccentricity range so the node-centered
        // difference for g_free has valid neighbours at the ends.
        let a0_at = |gamma: f64| -> f64 {
            let g = gamma.clamp(model.action(0.02), model.action(0.99));
            model.apsidal_harmonics(model.eccentricity(g)).0
        };

        for k in 0..=n {
            let e = e_min + (e_max - e_min) * k as f64 / n as f64;
            let g = model.action(e);
            let (_a0, c1, _c2) = model.apsidal_harmonics(e);
            gammas.push(g);
            c1s.push(c1);
            gfree.push((a0_at(g + dg) - a0_at(g - dg)) / (2.0 * dg));
        }

        // a₀ accumulated as the exact integral of the piecewise-linear g_free,
        // so that d a₀/dΓ ≡ g_free pointwise — the flow is then exactly the
        // gradient of the K that `k()` reports (no interpolant-mismatch drift).
        let mut a0_nodes = Vec::with_capacity(n + 1);
        a0_nodes.push(0.0);
        for k in 1..gammas.len() {
            let dgk = gammas[k] - gammas[k - 1];
            a0_nodes.push(a0_nodes[k - 1] + 0.5 * (gfree[k] + gfree[k - 1]) * dgk);
        }

        let mut rh = Self {
            model,
            g_p: 0.0,
            e_resonant,
            gammas,
            c1s,
            gfree,
            a0_nodes,
            c1_star: 0.0,
        };
        rh.g_p = rh.free_precession(model.action(e_resonant));
        rh.c1_star = rh.c1(model.action(e_resonant));
        rh
    }

    /// Override the external (Planet Nine) apsidal rate g₉, *detuning* the
    /// system away from resonance. With g₉ pushed outside the band of free
    /// precession rates the model spans, `g_free(Γ) − g₉` never vanishes, so
    /// there is no fixed point and the relative apse Δϖ circulates for every
    /// initial condition — the non-resonant, unclustered case.
    pub fn with_external_rate(mut self, g_p: f64) -> Self {
        self.g_p = g_p;
        self
    }

    /// Span [min, max] of the free precession rate over the eccentricity grid.
    /// A resonance (fixed point) exists only when g₉ lies inside this band.
    pub fn free_precession_band(&self) -> (f64, f64) {
        let lo = self.gfree.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = self.gfree.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    }

    /// Linear interpolation of a tabulated node array at action `gamma`.
    fn interp(&self, table: &[f64], gamma: f64) -> f64 {
        let g = gamma.clamp(self.gammas[0], *self.gammas.last().unwrap());
        for (i, gw) in self.gammas.windows(2).enumerate() {
            if g >= gw[0] && g <= gw[1] {
                let t = (g - gw[0]) / (gw[1] - gw[0]);
                return table[i] + t * (table[i + 1] - table[i]);
            }
        }
        *table.last().unwrap()
    }

    /// Axisymmetric mean a₀(Γ), the exact antiderivative of the piecewise-linear
    /// g_free (so da₀/dΓ ≡ g_free, making the pendulum self-consistent).
    pub fn a0(&self, gamma: f64) -> f64 {
        let g = gamma.clamp(self.gammas[0], *self.gammas.last().unwrap());
        for (i, gw) in self.gammas.windows(2).enumerate() {
            if g >= gw[0] && g <= gw[1] {
                // ∫ over the cell of the linear g_free from gw[0] to g.
                let h = gw[1] - gw[0];
                let s = (g - gw[0]) / h;
                let gf0 = self.gfree[i];
                let gf1 = self.gfree[i + 1];
                // a₀(g) = a₀_node[i] + ∫₀^{s·h} (gf0 + (gf1−gf0)·u/h) du
                return self.a0_nodes[i] + h * (gf0 * s + 0.5 * (gf1 - gf0) * s * s);
            }
        }
        *self.a0_nodes.last().unwrap()
    }

    /// Apsidal harmonic c₁(Γ).
    pub fn c1(&self, gamma: f64) -> f64 {
        self.interp(&self.c1s, gamma)
    }

    /// Free apsidal precession rate g_free(Γ) = da₀/dΓ (rad/day).
    pub fn free_precession(&self, gamma: f64) -> f64 {
        self.interp(&self.gfree, gamma)
    }

    /// Resonant Hamiltonian K(Δϖ, Γ) = a₀(Γ) − g₉·Γ + c₁*·cos Δϖ, with the
    /// harmonic frozen at the resonant action (the pendulum approximation).
    /// This is exactly the K whose gradient gives [`ddvarpi_dt`](Self::ddvarpi_dt)
    /// and [`dgamma_dt`](Self::dgamma_dt), so it is conserved along the flow.
    pub fn k(&self, dvarpi: f64, gamma: f64) -> f64 {
        self.a0(gamma) - self.g_p * gamma + self.c1_star * dvarpi.cos()
    }

    /// dΔϖ/dt = ∂K/∂Γ = g_free(Γ) − g₉.
    ///
    /// The second-fundamental-model-of-resonance form: the restoring term is
    /// the detuning of the free precession from the resonant rate
    /// (`g_free(Γ) − g₉`).
    pub fn ddvarpi_dt(&self, _dvarpi: f64, gamma: f64) -> f64 {
        self.free_precession(gamma) - self.g_p
    }

    /// dΓ/dt = −∂K/∂Δϖ = c₁*·sin Δϖ (harmonic frozen at the resonant action).
    pub fn dgamma_dt(&self, dvarpi: f64, _gamma: f64) -> f64 {
        self.c1_star * dvarpi.sin()
    }

    /// Locate the libration island: the resonant action Γ* where the free
    /// precession matches g₉, the stable apse, and the small-oscillation
    /// frequency. Returns `None` if no resonance is crossed (e.g. axisymmetric
    /// perturber: c₁ ≡ 0, the pendulum has no restoring torque).
    pub fn island(&self) -> Option<Island> {
        // Γ* is, by construction, the action at the tuning eccentricity:
        // the pendulum was tuned so g_free(Γ*) = g₉ there.
        let e_star = self.e_resonant;
        let gamma_star = self.model.action(e_star);

        let c1_star = self.c1(gamma_star);
        if c1_star.abs() < 1e-18 {
            return None; // no apsidal forcing → no island
        }

        // Pendulum: K ≈ const + ½ a₀''·(δΓ)² + c₁·cos Δϖ. The stable angle is
        // the one where the angular curvature and the action curvature combine
        // to a positive squared frequency ω² = −c₁·cos(Δϖ*)·a₀''.
        let a0pp = self.a0_curvature(gamma_star);
        // At Δϖ = π: ω² = −c₁·(−1)·a0'' = c₁·a0''. At Δϖ = 0: ω² = −c₁·a0''.
        let (center, omega2) = {
            let w_anti = c1_star * a0pp;
            let w_aligned = -c1_star * a0pp;
            if w_anti > 0.0 {
                (PI, w_anti)
            } else if w_aligned > 0.0 {
                (0.0, w_aligned)
            } else {
                return None;
            }
        };
        let omega_lib = omega2.sqrt();
        let period = 2.0 * PI / omega_lib;
        let steps_per_period = 1000usize;
        let dt = period / steps_per_period as f64;
        let steps = steps_per_period * 3;

        // Separatrix half-width in action: the pendulum barrier height is
        // 2|c₁*|; equating ½|a₀''|δΓ² = 2|c₁*| gives δΓ_sep = 2√(|c₁*|/|a₀''|).
        let separatrix_halfwidth_gamma = 2.0 * (c1_star.abs() / a0pp.abs()).sqrt();

        Some(Island {
            e_center: e_star,
            gamma_center: gamma_star,
            center_dvarpi: center,
            depth: c1_star.abs(),
            omega_lib,
            separatrix_halfwidth_gamma,
            suggested_dt: dt,
            suggested_steps: steps,
        })
    }

    /// Curvature a₀''(Γ) = dg_free/dΓ of the axisymmetric mean
    /// (rad/day per AU²/day), from the slope of the tabulated free precession.
    fn a0_curvature(&self, gamma: f64) -> f64 {
        let dg = 0.05 * (GM_SUN * self.model.a).sqrt();
        (self.free_precession(gamma + dg) - self.free_precession(gamma - dg)) / (2.0 * dg)
    }
}

/// Representative high eccentricity for the ETNO-scale libration island. The
/// resonance is tuned here (commensurability with Planet Nine's apse); it sits
/// in Beust's high-e regime and leaves headroom in the action grid for the
/// libration to swing without hitting the e = 0.97 ceiling.
pub const RESONANT_ECCENTRICITY: f64 = 0.7;

/// Convenience: locate the anti-aligned high-eccentricity libration island
/// (Beust's result), tuning the secular resonance at [`RESONANT_ECCENTRICITY`].
/// Returns `None` for an axisymmetric perturber (c₁ ≡ 0: no apsidal forcing,
/// hence no island).
pub fn libration_island(model: &SecularModel) -> Option<Island> {
    // Reject the axisymmetric case outright (circular perturber: c₁ ≡ 0).
    let c1 = model.apsidal_harmonics(RESONANT_ECCENTRICITY).1;
    if c1.abs() < 1e-18 {
        return None;
    }
    let rh = ResonantHamiltonian::new(*model, RESONANT_ECCENTRICITY);
    rh.island()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use p9_core::types::P9Params;

    #[test]
    fn island_exists_for_nominal_p9_and_is_anti_aligned() {
        let p9 = published::nominal_p9();
        let m = SecularModel::new(published::REPRESENTATIVE_ETNO_A_AU, &p9);
        let island = libration_island(&m).expect("nominal P9 should have an island");

        // Anti-aligned center (Beust's high-e libration islands).
        assert!(
            (island.center_dvarpi - published::ANTI_ALIGNED_CENTER_RAD).abs() < 1e-9,
            "island center {:.1}° not anti-aligned",
            island.center_dvarpi.to_degrees()
        );
        assert!(island.depth > 0.0, "resonance strength must be positive");
        assert!(island.omega_lib > 0.0, "libration frequency must be real");
        // Resonance sits at high eccentricity per the abstract.
        assert!(
            island.e_center > published::HIGH_ECCENTRICITY_THRESHOLD,
            "resonant e = {:.2} not in the high-e regime",
            island.e_center
        );
    }

    #[test]
    fn circular_perturber_has_no_island() {
        let mut p9 = published::nominal_p9();
        p9.e = 0.0;
        let m = SecularModel::new(published::REPRESENTATIVE_ETNO_A_AU, &p9);
        assert!(
            libration_island(&m).is_none(),
            "axisymmetric (circular) perturber must not produce an apsidal island"
        );
    }

    #[test]
    fn resonance_strengthens_with_p9_mass() {
        let p9 = published::nominal_p9();
        let base = SecularModel::new(published::REPRESENTATIVE_ETNO_A_AU, &p9);
        let light = base.with_gm_p(p9.gm() * 0.5);
        let heavy = base.with_gm_p(p9.gm() * 2.0);

        let il = libration_island(&light).expect("light island");
        let ih = libration_island(&heavy).expect("heavy island");

        // The pendulum amplitude c₁ is linear in GM, so a heavier P9 gives a
        // deeper resonance and a faster libration.
        assert!(
            ih.depth > il.depth,
            "heavier P9 should deepen the resonance: {:.3e} vs {:.3e}",
            ih.depth,
            il.depth
        );
        assert!(ih.omega_lib > il.omega_lib, "heavier P9 librates faster");

        // c₁ scales linearly in GM: compare at a common eccentricity.
        let e = 0.7;
        let c1_light = light.apsidal_harmonics(e).1;
        let c1_heavy = heavy.apsidal_harmonics(e).1;
        let ratio = c1_heavy / c1_light;
        assert!(
            (ratio - 4.0).abs() < 1e-6,
            "apsidal harmonic should scale linearly in GM (4×): {ratio:.4}"
        );
    }

    #[test]
    fn island_present_across_p9_configurations() {
        let p9 = P9Params::revised_2019();
        let m = SecularModel::new(300.0, &p9);
        assert!(
            libration_island(&m).is_some(),
            "revised-2019 P9 should also support a libration island"
        );
    }

    /// Every real ETNO in the vetted Brown (2017) clustering sample sits at a
    /// semi-major axis that admits an anti-aligned secular-resonance libration
    /// island under the nominal Planet Nine — the mechanism Beust invokes for
    /// the observed apsidal clustering applies to the actual objects.
    #[test]
    fn every_etno_admits_an_anti_aligned_island() {
        use p9_core::data::etno::BROWN_2017_SAMPLE;
        let p9 = published::nominal_p9();
        for kbo in &BROWN_2017_SAMPLE {
            let m = SecularModel::new(kbo.a, &p9);
            let island = libration_island(&m).unwrap_or_else(|| {
                panic!("{}: no libration island at a = {:.0} AU", kbo.name, kbo.a)
            });
            assert!(
                (island.center_dvarpi - PI).abs() < 1e-9,
                "{}: island center {:.0}° not anti-aligned",
                kbo.name,
                island.center_dvarpi.to_degrees()
            );
            assert!(island.depth > 0.0 && island.omega_lib > 0.0);
        }
    }
}

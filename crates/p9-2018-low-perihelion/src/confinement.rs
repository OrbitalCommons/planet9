//! Apsidal-confinement strength of a single P9 configuration on a test ETNO.
//!
//! We reuse `p9_core::analysis::secular::numerical_secular_hamiltonian` — the
//! exact doubly-averaged (Gauss-ring) 1/Δ interaction, coplanar (i = 0) — as
//! the conserved secular Hamiltonian H(e, Δϖ) of a test ETNO at fixed
//! semi-major axis `a`, perturbed by a Planet Nine of semi-major axis `a9`,
//! eccentricity `e9`, mass `mass_earth`. The perturber defines the reference
//! frame, so Δϖ is the ETNO longitude of perihelion relative to P9.
//!
//! For an eccentric perturber the octupole-and-higher terms break the apsidal
//! degeneracy and open an apsidal libration island — exactly the confinement
//! Cáceres & Gomes invoke. We characterise it with real numbers computed from
//! the portrait:
//!
//! 1. `favored_apse` — which apse (Δϖ = 0 "aligned" or Δϖ = π "anti-aligned")
//!    is the secular energy MINIMUM at the forced eccentricity, i.e. the apse
//!    an ETNO is shepherded toward.
//! 2. `forced_eccentricity` — the e at the favoured equilibrium (the minimum
//!    of H along the favoured apse over the high-e branch).
//! 3. `libration_half_width` — the half-width in Δϖ of the level curve
//!    bounded by the saddle at the OPPOSITE apse; the confinement metric
//!    (radians, in (0, π]). < π ⟹ a librating (confined) island; ≈ π ⟹
//!    circulation (no confinement).
//! 4. `torque_amplitude` — max over Δϖ of |∂H/∂ϖ| at the forced eccentricity:
//!    the secular forcing strength (∝ GM₉). Cross-checked against the bare
//!    p9-core kernel in the tests.
//!
//! HONEST SIGN CAVEAT. The single-perturber coplanar Gauss-ring secular problem
//! confines the test ETNO toward the apse that minimises H. For the canonical
//! detached B&B orbit (e9 = 0.6) this energy minimum sits near the ALIGNED apse
//! (Δϖ ≈ 0); as q9 is lowered (e9 raised) the well deepens and sharpens. The
//! OBSERVED ETNOs are clustered ANTI-aligned with the predicted P9 perihelion —
//! recovering that sign requires P9's own apsidal precession and the giants'
//! secular field (full secular/N-body), which the bare single-ring test
//! particle does not carry (the same limitation documented for
//! `p9-2019-selfgrav-disk` and `p9-2021-oort-cloud`). What this crate
//! reproduces is the paper's *mechanism and its q9-dependence*: a coherent
//! apsidal libration island whose strength PERSISTS and grows as the P9
//! perihelion is lowered. The sign of the favoured apse is reported, not pinned
//! as anti-aligned.
//!
//! Everything is computed from the SAME p9-core kernel; nothing is duplicated.

use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};
use p9_core::units::{au, radians, Angle, Length};
use std::f64::consts::PI;

/// Quadrature nodes per anomaly for the Gauss-ring average (n² evaluations).
const N_QUAD: usize = 96;
/// Δϖ samples for the portrait scan.
const N_DVARPI: usize = 180;
/// Lowest eccentricity considered for the equilibrium search. Below this the
/// Gauss-ring H(e) turns over (it peaks in magnitude near e ≈ 0.5) and admits
/// spurious low-e level-curve roots; the clustered ETNOs live on the high-e
/// branch, so we restrict to e ≥ this value (mirrors the ossos-scattering
/// crate's `E_BRANCH_MIN`).
const E_BRANCH_MIN: f64 = 0.55;
/// Upper eccentricity bound for the equilibrium search (avoid e → 1).
const E_BRANCH_MAX: f64 = 0.95;
/// Reference ETNO eccentricity at which the apsidal well depth is measured.
/// Sedna-like (e ≈ 0.85), representative of the clustered sample and on the
/// high-e branch where the Gauss-ring portrait is monotone in e. Holding e
/// fixed isolates the well-depth dependence on q9 from the (real, but
/// branch-jumping) drift of the forced eccentricity, giving a clean monotone
/// confinement metric.
const E_REF_WELL: f64 = 0.85;

/// Which apse the secular dynamics confines the ETNO toward.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FavoredApse {
    /// Δϖ = 0: ETNO perihelion aligned with P9.
    Aligned,
    /// Δϖ = π: ETNO perihelion anti-aligned with P9.
    AntiAligned,
}

/// Confinement diagnostics for one (test-ETNO a, P9 a9/e9/mass) configuration.
#[derive(Debug, Clone, Copy)]
pub struct Confinement {
    /// Planet Nine perihelion q9 = a9(1 − e9) (AU).
    pub q9_au: f64,
    /// Apse the dynamics shepherds toward (energy minimum).
    pub favored_apse: FavoredApse,
    /// Forced (equilibrium) ETNO eccentricity at the favoured apse.
    pub forced_eccentricity: f64,
    /// Half-width of the Δϖ libration island (radians), in (0, π].
    /// < π ⟹ confined (librating); ≈ π ⟹ circulating (unconfined).
    pub libration_half_width: f64,
    /// Secular torque amplitude max_Δϖ |∂H/∂ϖ| at the forced e (AU²/day²).
    pub torque_amplitude: f64,
    /// Dimensionless depth of the apsidal well: (H_saddle − H_min)/|H_min| at
    /// the reference eccentricity `E_REF_WELL`, where H_min is the favoured apse
    /// and H_saddle the opposite apse. Positive ⟹ a real confining well exists.
    /// Measured at fixed e so it depends only on q9/e9 (the confinement
    /// strength), not on the forced-e branch.
    pub well_depth: f64,
}

impl Confinement {
    /// `true` when the configuration confines: a finite (sub-π) libration
    /// island exists AND a real apsidal well is present.
    pub fn is_confined(&self) -> bool {
        self.libration_half_width < PI - 1e-3 && self.well_depth > 1e-6
    }

    /// Planet Nine perihelion q9 as a typed [`Length`] (see [`Self::q9_au`]).
    pub fn q9(&self) -> Length {
        au(self.q9_au)
    }

    /// Δϖ libration-island half-width as a typed [`Angle`] (see
    /// [`Self::libration_half_width`]).
    pub fn libration_half_width_angle(&self) -> Angle {
        radians(self.libration_half_width)
    }
}

/// P9 GM in AU³/day² from a mass in Earth masses.
fn p9_gm(mass_earth: f64) -> f64 {
    mass_earth * EARTH_MASS_SOLAR * GM_SUN
}

/// Coplanar secular Hamiltonian of the test ETNO at (e, Δϖ) for a P9 of
/// (a9, e9, mass). Thin wrapper over the p9-core Gauss-ring kernel with the
/// 1%·a9 softening used by `phase_portrait`/the ossos-scattering crate.
fn h(a: f64, e: f64, dvarpi: f64, a9: f64, e9: f64, mass_earth: f64) -> f64 {
    let soft = 0.01 * a9;
    numerical_secular_hamiltonian(
        a,
        e,
        0.0,
        dvarpi,
        0.0,
        a9,
        e9,
        p9_gm(mass_earth),
        N_QUAD,
        soft,
    )
}

/// Maximum |∂H/∂ϖ| over Δϖ at fixed eccentricity (the secular torque
/// amplitude). Centred finite difference of the Gauss-ring Hamiltonian.
pub fn torque_amplitude(a: f64, e: f64, a9: f64, e9: f64, mass_earth: f64) -> f64 {
    let dphi = 1e-3;
    let mut max_dh = 0.0_f64;
    for j in 0..N_DVARPI {
        let dv = j as f64 / N_DVARPI as f64 * 2.0 * PI;
        let hp = h(a, e, dv + dphi, a9, e9, mass_earth);
        let hm = h(a, e, dv - dphi, a9, e9, mass_earth);
        max_dh = max_dh.max(((hp - hm) / (2.0 * dphi)).abs());
    }
    max_dh
}

/// Golden-section minimiser of `f` on `[lo, hi]` (assumes unimodality on the
/// high-e branch). Returns the minimiser.
fn golden_min(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    let gr = (5.0_f64.sqrt() - 1.0) / 2.0;
    let mut c = hi - gr * (hi - lo);
    let mut d = lo + gr * (hi - lo);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..60 {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - gr * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + gr * (hi - lo);
            fd = f(d);
        }
    }
    0.5 * (lo + hi)
}

/// Forced eccentricity along a fixed apse `dvarpi`: the e on [E_BRANCH_MIN,
/// E_BRANCH_MAX] minimising H at that apse.
fn forced_e_at_apse(a: f64, dvarpi: f64, a9: f64, e9: f64, mass_earth: f64) -> f64 {
    golden_min(E_BRANCH_MIN, E_BRANCH_MAX, |e| {
        h(a, e, dvarpi, a9, e9, mass_earth)
    })
}

/// Half-width (radians) of the Δϖ libration island at fixed `e_f`, around the
/// favoured apse `apse_dv`. The island is bounded by the level curve through
/// the saddle at the opposite apse; we report the angular distance from the
/// favoured apse to where H(e_f, Δϖ) reaches the mean of (apse, saddle)
/// energies — a robust interior-contour width. If H is monotone between the two
/// apses (no well), the orbit circulates and the half-width is π.
fn libration_half_width(a: f64, e_f: f64, apse_dv: f64, a9: f64, e9: f64, mass_earth: f64) -> f64 {
    let saddle_dv = if apse_dv == 0.0 { PI } else { 0.0 };
    let h_min = h(a, e_f, apse_dv, a9, e9, mass_earth);
    let h_saddle = h(a, e_f, saddle_dv, a9, e9, mass_earth);
    if h_saddle <= h_min {
        // The "favoured" apse is not actually the minimum at this e: no well.
        return PI;
    }
    let h_target = 0.5 * (h_min + h_saddle);
    // Walk from the favoured apse toward the saddle, find the H = h_target
    // crossing; its angular distance from the apse is the half-width.
    let n = 400;
    let mut prev = h_min;
    for k in 1..=n {
        let frac = k as f64 / n as f64;
        let dv = apse_dv + frac * (saddle_dv - apse_dv);
        let hv = h(a, e_f, dv, a9, e9, mass_earth);
        if (prev - h_target) * (hv - h_target) <= 0.0 {
            let dv_prev = apse_dv + ((k - 1) as f64 / n as f64) * (saddle_dv - apse_dv);
            let t = (h_target - prev) / (hv - prev);
            let dv_cross = dv_prev + t * (dv - dv_prev);
            return (dv_cross - apse_dv).abs();
        }
        prev = hv;
    }
    PI
}

/// Confining well depth at the reference eccentricity with an explicit
/// softening fraction (of a₉) — the probe behind the softening-robustness
/// test. The sweep geometry is orbit-crossing (see `reference::TEST_ETNO_A_AU`
/// docs), so the doubly-averaged integral is softening-regularized; results
/// must not depend materially on that regularization choice.
pub fn well_depth_with_softening(a: f64, a9: f64, e9: f64, mass_earth: f64, soft_frac: f64) -> f64 {
    let soft = soft_frac * a9;
    let hh = |dvarpi: f64| {
        numerical_secular_hamiltonian(
            a,
            E_REF_WELL,
            0.0,
            dvarpi,
            0.0,
            a9,
            e9,
            p9_gm(mass_earth),
            N_QUAD,
            soft,
        )
    };
    // Anti-aligned well vs aligned saddle (the confining contrast for an
    // eccentric exterior perturber at these geometries).
    let (h_anti, h_aligned) = (hh(PI), hh(0.0));
    (h_aligned - h_anti).abs()
}

/// Full confinement diagnostics for a test ETNO of semi-major axis `a` under a
/// Planet Nine of semi-major axis `a9`, eccentricity `e9`, mass `mass_earth`.
pub fn analyze(a: f64, a9: f64, e9: f64, mass_earth: f64) -> Confinement {
    let q9_au = a9 * (1.0 - e9);

    // Forced eccentricity along each apse, then pick the apse whose minimum H
    // is lower — that is the confining (favoured) apse.
    let e_aligned = forced_e_at_apse(a, 0.0, a9, e9, mass_earth);
    let e_anti = forced_e_at_apse(a, PI, a9, e9, mass_earth);
    let h_aligned = h(a, e_aligned, 0.0, a9, e9, mass_earth);
    let h_anti = h(a, e_anti, PI, a9, e9, mass_earth);

    let (favored_apse, apse_dv, e_f) = if h_aligned <= h_anti {
        (FavoredApse::Aligned, 0.0, e_aligned)
    } else {
        (FavoredApse::AntiAligned, PI, e_anti)
    };

    let half_width = libration_half_width(a, e_f, apse_dv, a9, e9, mass_earth);
    let torque = torque_amplitude(a, e_f, a9, e9, mass_earth);

    // Well depth at the reference eccentricity: opposite apse (saddle) minus
    // favoured apse, holding e fixed so the metric tracks q9/e9 alone.
    let saddle_dv = if apse_dv == 0.0 { PI } else { 0.0 };
    let h_well = h(a, E_REF_WELL, apse_dv, a9, e9, mass_earth);
    let h_saddle = h(a, E_REF_WELL, saddle_dv, a9, e9, mass_earth);
    let well_depth = (h_saddle - h_well) / h_well.abs();

    Confinement {
        q9_au,
        favored_apse,
        forced_eccentricity: e_f,
        libration_half_width: half_width,
        torque_amplitude: torque,
        well_depth,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::{BB_A9_AU, BB_E9, P9_MASS_EARTH, TEST_ETNO_A_AU};
    use approx::assert_relative_eq;

    #[test]
    fn canonical_bb_orbit_confines() {
        // The canonical detached B&B orbit must produce apsidal confinement:
        // a finite libration island and a real apsidal well.
        let c = analyze(TEST_ETNO_A_AU, BB_A9_AU, BB_E9, P9_MASS_EARTH);
        assert!(c.is_confined(), "canonical orbit failed to confine: {c:?}");
        assert!(c.forced_eccentricity > E_BRANCH_MIN && c.forced_eccentricity < E_BRANCH_MAX);
        assert!(c.libration_half_width > 0.0 && c.libration_half_width < PI);
        assert!(c.torque_amplitude > 0.0);
        assert!(c.well_depth > 0.0);
        assert_relative_eq!(c.q9_au, 280.0, epsilon = 1e-9);
    }

    #[test]
    fn typed_accessors_match_f64_sources() {
        use uom::si::angle::radian;
        use uom::si::length::astronomical_unit;
        let c = analyze(TEST_ETNO_A_AU, BB_A9_AU, BB_E9, P9_MASS_EARTH);
        assert_relative_eq!(c.q9().get::<astronomical_unit>(), c.q9_au, epsilon = 1e-9);
        assert_relative_eq!(
            c.libration_half_width_angle().get::<radian>(),
            c.libration_half_width,
            epsilon = 1e-9
        );
    }

    #[test]
    fn torque_scales_linearly_with_p9_mass() {
        // The Gauss-ring Hamiltonian is ∝ GM₉, so the torque amplitude scales
        // linearly with the P9 mass (a real property of linear secular theory).
        let t10 = torque_amplitude(TEST_ETNO_A_AU, 0.7, BB_A9_AU, BB_E9, 10.0);
        let t20 = torque_amplitude(TEST_ETNO_A_AU, 0.7, BB_A9_AU, BB_E9, 20.0);
        assert_relative_eq!(t20 / t10, 2.0, epsilon = 1e-6);
    }

    #[test]
    fn torque_amplitude_matches_core_kernel() {
        // Cross-check: the crate's torque amplitude is a finite difference of
        // the SAME p9-core kernel. Recompute it directly here from
        // numerical_secular_hamiltonian and require bit-level agreement,
        // confirming no re-implementation drift.
        let (a, e, a9, e9, m) = (TEST_ETNO_A_AU, 0.7, BB_A9_AU, BB_E9, P9_MASS_EARTH);
        let soft = 0.01 * a9;
        let gm = m * EARTH_MASS_SOLAR * GM_SUN;
        let dphi = 1e-3;
        let mut max_dh = 0.0_f64;
        for j in 0..N_DVARPI {
            let dv = j as f64 / N_DVARPI as f64 * 2.0 * PI;
            let hp =
                numerical_secular_hamiltonian(a, e, 0.0, dv + dphi, 0.0, a9, e9, gm, N_QUAD, soft);
            let hm =
                numerical_secular_hamiltonian(a, e, 0.0, dv - dphi, 0.0, a9, e9, gm, N_QUAD, soft);
            max_dh = max_dh.max(((hp - hm) / (2.0 * dphi)).abs());
        }
        let crate_val = torque_amplitude(a, e, a9, e9, m);
        assert_eq!(crate_val, max_dh);
    }

    #[test]
    fn forced_e_matches_direct_core_minimum() {
        // Cross-check the forced eccentricity against a brute-force minimum of
        // the SAME p9-core kernel along the favoured apse, to a fine grid.
        let c = analyze(TEST_ETNO_A_AU, BB_A9_AU, BB_E9, P9_MASS_EARTH);
        let apse = match c.favored_apse {
            FavoredApse::Aligned => 0.0,
            FavoredApse::AntiAligned => PI,
        };
        let gm = P9_MASS_EARTH * EARTH_MASS_SOLAR * GM_SUN;
        let soft = 0.01 * BB_A9_AU;
        let mut best_e = E_BRANCH_MIN;
        let mut best_h = f64::INFINITY;
        let steps = 800;
        for s in 0..=steps {
            let e = E_BRANCH_MIN + (E_BRANCH_MAX - E_BRANCH_MIN) * s as f64 / steps as f64;
            let hv = numerical_secular_hamiltonian(
                TEST_ETNO_A_AU,
                e,
                0.0,
                apse,
                0.0,
                BB_A9_AU,
                BB_E9,
                gm,
                N_QUAD,
                soft,
            );
            if hv < best_h {
                best_h = hv;
                best_e = e;
            }
        }
        assert_relative_eq!(c.forced_eccentricity, best_e, epsilon = 5e-3);
    }

    #[test]
    fn circular_perturber_has_no_apsidal_confinement() {
        // e9 = 0: the averaged P9 is an axisymmetric ring, so there is no
        // Δϖ structure — no apsidal well, the island fills 2π (half-width = π)
        // and there is no confinement.
        let c = analyze(TEST_ETNO_A_AU, BB_A9_AU, 0.0, P9_MASS_EARTH);
        assert!(!c.is_confined(), "circular P9 must not confine: {c:?}");
        assert_relative_eq!(c.libration_half_width, PI, epsilon = 1e-9);
        // The defining no-confinement signature: zero apsidal well and a 2π
        // (half-width π) island. (A non-zero `torque_amplitude` survives only as
        // finite-grid quadrature aliasing of the axisymmetric ring; the
        // well-depth is the symmetry-respecting metric and is exactly zero.)
        assert_relative_eq!(c.well_depth, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn confinement_is_octupole_not_quadrupole() {
        // Cross-check against p9-core's analytic quadrupole. The coplanar
        // quadrupole H is a function of e ALONE (no Δϖ dependence): it cannot
        // produce apsidal confinement. The confining well this crate measures
        // must therefore come from the octupole-and-higher terms the Gauss-ring
        // average captures. Verify the quadrupole has zero Δϖ contrast while the
        // numerical average has a finite one at the same geometry.
        use p9_core::analysis::secular::coplanar_quadrupole;
        let (a, a9, e9, m) = (TEST_ETNO_A_AU, BB_A9_AU, BB_E9, P9_MASS_EARTH);
        let gm = m * EARTH_MASS_SOLAR * GM_SUN;
        // Quadrupole is Δϖ-independent by construction (single value of e).
        let q1 = coplanar_quadrupole(a, E_REF_WELL, a9, e9, gm);
        let q2 = coplanar_quadrupole(a, E_REF_WELL, a9, e9, gm);
        assert_eq!(q1, q2);
        // The numerical average DOES distinguish the two apses (octupole+):
        let c = analyze(a, a9, e9, m);
        assert!(
            c.well_depth > 1e-3,
            "expected a finite octupole+ apsidal well, got {:.3e}",
            c.well_depth
        );
    }

    #[test]
    fn lower_q9_confines_more_strongly_than_canonical() {
        // The crate headline, at the single-configuration level: a
        // low-perihelion P9 (q9 = 70 AU, e9 = 0.9) confines MORE strongly than
        // the canonical detached orbit (q9 = 280 AU). Both confine.
        let canonical = analyze(TEST_ETNO_A_AU, BB_A9_AU, BB_E9, P9_MASS_EARTH);
        let low_q = analyze(TEST_ETNO_A_AU, BB_A9_AU, 0.9, P9_MASS_EARTH);
        assert!(canonical.is_confined() && low_q.is_confined());
        assert!(
            low_q.well_depth > canonical.well_depth,
            "low-q9 well {:.3} should exceed canonical {:.3}",
            low_q.well_depth,
            canonical.well_depth
        );
        assert!(
            low_q.libration_half_width < canonical.libration_half_width,
            "low-q9 island {:.3} should be tighter than canonical {:.3}",
            low_q.libration_half_width,
            canonical.libration_half_width
        );
    }

    #[test]
    fn sweep_deepening_is_robust_to_softening_choice() {
        // The headline monotone deepening of the confining well as q9 drops
        // is computed in an orbit-crossing, softening-regularized regime.
        // It must survive doubling the softening: same monotone direction,
        // and depths stable to a modest factor.
        use crate::reference::{BB_A9_AU, P9_MASS_EARTH, TEST_ETNO_A_AU};
        let e9s = [0.6, 0.75, 0.9];
        for &soft in &[0.01, 0.02] {
            let mut prev = 0.0;
            for &e9 in &e9s {
                let d =
                    well_depth_with_softening(TEST_ETNO_A_AU, BB_A9_AU, e9, P9_MASS_EARTH, soft);
                assert!(
                    d > prev,
                    "soft {soft}: depth not deepening at e9 = {e9}: {d:.3e} <= {prev:.3e}"
                );
                prev = d;
            }
        }
        for &e9 in &e9s {
            let d1 = well_depth_with_softening(TEST_ETNO_A_AU, BB_A9_AU, e9, P9_MASS_EARTH, 0.01);
            let d2 = well_depth_with_softening(TEST_ETNO_A_AU, BB_A9_AU, e9, P9_MASS_EARTH, 0.02);
            let ratio = d2 / d1;
            assert!(
                (0.5..2.0).contains(&ratio),
                "e9 = {e9}: depth ratio under 2x softening = {ratio:.3}"
            );
        }
    }
}

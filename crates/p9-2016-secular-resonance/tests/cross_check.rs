//! Cross-check the crate's secular Hamiltonian against the p9-core ring
//! average it is built on, and confirm the qualitative Beust (2016) result
//! end to end.

use p9_2016_secular_resonance::{libration_island, published, SecularModel};
use p9_core::analysis::secular::numerical_secular_hamiltonian;
use std::f64::consts::PI;

/// `SecularModel::hamiltonian` must equal `p9_core`'s numerical Gauss-ring
/// average at the same arguments (it is a thin coplanar wrapper). We re-derive
/// the value independently here and require agreement to a tight tolerance.
#[test]
fn hamiltonian_matches_p9_core_ring_average() {
    let p9 = published::nominal_p9();
    let m = SecularModel::new(published::REPRESENTATIVE_ETNO_A_AU, &p9);

    for &e in &[0.3, 0.6, 0.85] {
        for &dv in &[0.0, PI / 2.0, PI] {
            let ours = m.hamiltonian(e, dv);
            let core = numerical_secular_hamiltonian(
                m.a,
                e,
                0.0,
                dv,
                0.0,
                m.a_p,
                m.e_p,
                m.gm_p,
                m.n_quad,
                m.softening,
            );
            let rel = ((ours - core) / core).abs();
            assert!(
                rel < 1e-12,
                "wrapper vs p9-core ring average disagree at e={e}, Δϖ={dv}: rel {rel:.2e}"
            );
        }
    }
}

/// Quadrature-convergence cross-check: the apsidal harmonic c₁ (which builds
/// the libration island) and the resulting island depth must be robust to
/// refining the ring-average quadrature, so the island is a physical feature,
/// not a quadrature artifact.
#[test]
fn island_robust_to_quadrature_refinement() {
    let p9 = published::nominal_p9();
    let mut coarse = SecularModel::new(400.0, &p9);
    coarse.n_quad = 48;
    let mut fine = coarse;
    fine.n_quad = 96;

    // The apsidal harmonic at the resonant eccentricity converges with the
    // ring quadrature.
    let c1_c = coarse.apsidal_harmonics(0.7).1;
    let c1_f = fine.apsidal_harmonics(0.7).1;
    let rel_c1 = ((c1_c - c1_f) / c1_f).abs();
    assert!(rel_c1 < 0.02, "c₁ not converged: rel {rel_c1:.3e}");

    // And the located island depth (= |c₁*|) is correspondingly stable, and
    // both quadratures agree the center is anti-aligned.
    let ic = libration_island(&coarse).expect("coarse island");
    let i_f = libration_island(&fine).expect("fine island");
    assert!((ic.center_dvarpi - std::f64::consts::PI).abs() < 1e-9);
    assert!((i_f.center_dvarpi - std::f64::consts::PI).abs() < 1e-9);
    let rel = ((ic.depth - i_f.depth) / i_f.depth).abs();
    assert!(rel < 0.03, "island depth not converged: rel {rel:.3e}");
}

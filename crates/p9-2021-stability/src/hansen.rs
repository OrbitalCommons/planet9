//! Hansen coefficients for the disturbing function.
//!
//! Thin wrapper over the single workspace implementation in
//! `p9_core::analysis::hansen` (this crate previously carried its own copy
//! with a dimensionally wrong `hansen_x_neg3_0` asymptotic — `j·exp(−(1−e)²)`
//! used q/a where the paper's length scale is q/a_N — and an unguarded
//! Newton Kepler loop).
//!
//! The asymptotic X_j^{-3,2} ≈ (2j/5)·exp(−(q/a_N)²) is the form used by
//! Batygin, Mardling & Nesvorny (2021) for the 2:j chain; the cross-validation
//! test below checks it against the convergence-controlled numerical
//! coefficient on the chain itself.

#[cfg(test)]
mod tests {
    use crate::resonance_chain::resonance_semimajor;
    use p9_core::analysis::hansen::{
        hansen_coefficient, hansen_x0_neg3_0, hansen_x_neg3_2_neptune_chain, mean_to_true_anomaly,
    };

    #[test]
    fn test_asymptotic_trends() {
        // Positive, increasing in j, decreasing in q.
        assert!(hansen_x_neg3_2_neptune_chain(20, 35.0) > 0.0);
        assert!(hansen_x_neg3_2_neptune_chain(30, 35.0) > hansen_x_neg3_2_neptune_chain(20, 35.0));
        assert!(hansen_x_neg3_2_neptune_chain(20, 30.0) > hansen_x_neg3_2_neptune_chain(20, 40.0));
    }

    #[test]
    fn test_numerical_matches_closed_form() {
        // X_0^{-3,0}(e) = (1 − e²)^{−3/2}
        for &e in &[0.0, 0.5, 0.9] {
            let num = hansen_coefficient(0, -3, 0, e, 1e-10);
            let exact = hansen_x0_neg3_0(e);
            assert!(
                ((num - exact) / exact).abs() < 1e-8,
                "e = {e}: {num} vs {exact}"
            );
        }
    }

    /// M4: certify the chain asymptotic against the numerical Hansen
    /// coefficient on the resonances it is actually used for (the previous
    /// local copy was never compared to anything, which is how the
    /// dimensionally wrong X_j^{-3,0} survived).
    #[test]
    fn test_asymptotic_vs_numerical_on_chain() {
        for &j in &[15_i64, 20, 30] {
            for &q in &[32.0, 38.0] {
                let a = resonance_semimajor(j);
                let e = 1.0 - q / a;
                let approx = hansen_x_neg3_2_neptune_chain(j, q);
                let numerical = hansen_coefficient(j, -3, 2, e, 1e-8);
                let ratio = approx / numerical;
                assert!(
                    (0.5..2.0).contains(&ratio),
                    "j = {j}, q = {q}: approx {approx:.4} vs numerical {numerical:.4} (ratio {ratio:.3})"
                );
            }
        }
    }
}

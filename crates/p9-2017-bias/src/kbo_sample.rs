//! KBO sample from Brown (2017).
//!
//! Thin re-export of the single vetted ETNO table in
//! `p9_core::data::etno`. The paper's primary sample is the 10 objects with
//! a ≳ 230 AU; note that 2000 CR105 (a ≈ 228 AU) is included as in the
//! paper even though it sits just below the nominal 230 AU cut — the sample
//! membership follows the paper's Table 1, not a strict numerical threshold.

pub use p9_core::data::etno::{Etno, BROWN_2017_SAMPLE};

/// The 10 KBOs used in the paper's primary analysis (Brown 2017, Table 1).
pub fn paper_sample_a230() -> Vec<Etno> {
    BROWN_2017_SAMPLE.to_vec()
}

/// Extract longitude of perihelion ϖ = ω + Ω for each KBO (in radians).
pub fn longitudes_of_perihelion(kbos: &[Etno]) -> Vec<f64> {
    kbos.iter().map(|k| k.longitude_of_perihelion()).collect()
}

/// Extract arguments of perihelion ω for each KBO (in radians).
pub fn arguments_of_perihelion(kbos: &[Etno]) -> Vec<f64> {
    kbos.iter()
        .map(|k| k.omega_deg * p9_core::constants::DEG2RAD)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::analysis::circular::mean_resultant_length;

    #[test]
    fn test_sample_size() {
        let kbos = paper_sample_a230();
        assert_eq!(kbos.len(), 10);
    }

    #[test]
    fn test_all_distant() {
        // The paper's sample is "a > 230 AU" with 2000 CR105 (a ≈ 228 AU)
        // included; everything else clears 230 AU.
        for kbo in &paper_sample_a230() {
            if kbo.name == "2000 CR105" {
                assert!((kbo.a - 228.0).abs() < 1.0);
            } else {
                assert!(
                    kbo.a > 230.0,
                    "{} has a = {:.0} (expected > 230)",
                    kbo.name,
                    kbo.a
                );
            }
        }
    }

    #[test]
    fn test_all_eccentric() {
        for kbo in &paper_sample_a230() {
            assert!(
                kbo.e > 0.6,
                "{} has e = {:.2} (expected > 0.6)",
                kbo.name,
                kbo.e
            );
        }
    }

    #[test]
    fn test_varpi_clustering() {
        let kbos = paper_sample_a230();
        let varpis = longitudes_of_perihelion(&kbos);
        let r_bar = mean_resultant_length(&varpis);
        assert!(
            r_bar > 0.3,
            "R̄ = {:.3} should indicate clustering (> 0.3)",
            r_bar
        );
    }
}

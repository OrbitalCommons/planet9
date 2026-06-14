//! The two new 2025 objects encoded as documented `Etno`-style entries.
//!
//! Provenance of every element below: NASA/JPL Small-Body Database (SBDB)
//! osculating solutions, fetched 2026-06-14, epoch JD 2461200.5 (TDB),
//! frame ECLIPJ2000 (same heliocentric ecliptic frame as
//! `p9_core::data::etno::BROWN_2017_SAMPLE`). The discovery-paper abstracts
//! quote only rounded headline numbers (830 AU / 45 AU for 2017 OF201;
//! 252 AU / 66 AU / 11° for Ammonite); the full angular elements ω and Ω —
//! which the ϖ clustering analysis needs — come from the JPL fits, which are
//! the same arc the papers report.

use p9_core::data::etno::Etno;

/// A new discovery with its discovery-paper headline reference numbers kept
/// separately from the computed orbit.
#[derive(Debug, Clone, Copy)]
pub struct NewDiscovery {
    /// Vetted orbital elements (reuses the workspace `Etno` type).
    pub etno: Etno,
    /// arXiv identifier of the discovery paper.
    pub arxiv: &'static str,
    /// Semi-major axis (AU) as quoted in the paper headline, kept as a
    /// labelled reference value (the JPL osculating `a` differs slightly).
    pub paper_a_au: f64,
    /// Perihelion distance q (AU) as quoted in the paper headline.
    pub paper_q_au: f64,
    /// One-line note on why this object stresses the clustering hypothesis.
    pub note: &'static str,
}

/// 2017 OF201 — Cheng, Yang, et al. 2025 (arXiv:2505.15806).
///
/// JPL SBDB solution 2 (soln 2025-07-03, 25 obs, arc 2004–2018):
/// a = 823.64 AU, e = 0.94519, i = 16.221°, ω = 337.835°, Ω = 328.760°.
/// Paper headline: a ≈ 830 AU, q ≈ 45 AU, dwarf-planet-sized (H ≈ 3.5).
/// Its ϖ ≈ 307° sits far from the clustered ETNO mean (~52°): an object that
/// Planet Nine should herd, yet is not aligned.
pub fn of201() -> NewDiscovery {
    NewDiscovery {
        etno: Etno {
            name: "2017 OF201",
            a: 823.637,
            e: 0.945188,
            i_deg: 16.221,
            omega_deg: 337.835,
            omega_big_deg: 328.760,
            h_mag: 3.49,
        },
        arxiv: "2505.15806",
        paper_a_au: 830.0,
        paper_q_au: 45.0,
        note: "very distant ETNO with ϖ outside the clustered group",
    }
}

/// Ammonite / 2023 KQ14 — Chen, Kavelaars, et al. 2025 (arXiv:2508.02162).
///
/// JPL SBDB solution 3 (soln 2026-05-26, 88 obs, arc 2005–2025):
/// a = 247.43 AU, e = 0.73357, i = 11.001°, ω = 198.704°, Ω = 72.079°.
/// Paper headline: a ≈ 252 AU, q ≈ 66 AU, i ≈ 11° (H ≈ 6.98). Reported as the
/// FOURTH sednoid, but with ϖ ≈ 271° misaligned with the other three sednoids,
/// reducing the apparent sednoid apsidal alignment.
pub fn ammonite() -> NewDiscovery {
    NewDiscovery {
        etno: Etno {
            name: "Ammonite (2023 KQ14)",
            a: 247.429,
            e: 0.733566,
            i_deg: 11.001,
            omega_deg: 198.704,
            omega_big_deg: 72.079,
            h_mag: 6.98,
        },
        arxiv: "2508.02162",
        paper_a_au: 252.0,
        paper_q_au: 66.0,
        note: "fourth sednoid whose ϖ is misaligned with the other three",
    }
}

/// Published headline semi-major axis for 2017 OF201 (AU), Cheng et al. 2025.
pub const OF201_PAPER_A_AU: f64 = 830.0;

/// Ammonite is reported as the FOURTH (misaligned) sednoid, Chen et al. 2025.
pub const AMMONITE_SEDNOID_RANK: u32 = 4;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn of201_is_very_distant() {
        let d = of201();
        // Paper headline a ≈ 830 AU; JPL osculating a within ~1% of it.
        assert_relative_eq!(d.paper_a_au, OF201_PAPER_A_AU);
        assert!(
            (d.etno.a - OF201_PAPER_A_AU).abs() < 10.0,
            "JPL a = {} AU vs paper 830 AU",
            d.etno.a
        );
        assert!(d.etno.a > 800.0, "a = {} AU", d.etno.a);
        // Perihelion q = a(1-e) matches the published q ≈ 45 AU.
        assert!(
            (d.etno.perihelion() - d.paper_q_au).abs() < 1.0,
            "q = {:.2} AU vs paper {} AU",
            d.etno.perihelion(),
            d.paper_q_au
        );
        assert!(d.etno.e > 0.9, "e = {}", d.etno.e);
    }

    #[test]
    fn ammonite_is_a_sednoid() {
        let d = ammonite();
        assert_eq!(AMMONITE_SEDNOID_RANK, 4);
        // Sednoid hallmark: detached, q well beyond Neptune's reach.
        assert!(
            d.etno.perihelion() > 50.0,
            "q = {:.2} AU",
            d.etno.perihelion()
        );
        assert!(
            (d.etno.perihelion() - d.paper_q_au).abs() < 1.0,
            "q = {:.2} AU vs paper {} AU",
            d.etno.perihelion(),
            d.paper_q_au
        );
        // Paper headline a ≈ 252 AU; JPL osculating within a few AU.
        assert!((d.etno.a - 252.0).abs() < 6.0, "a = {} AU", d.etno.a);
        assert!((d.etno.i_deg - 11.0).abs() < 0.5, "i = {}", d.etno.i_deg);
    }
}

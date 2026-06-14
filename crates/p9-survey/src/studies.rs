//! The catalog of Planet Nine orbit solutions surveyed here.
//!
//! Every base orbit is read from a *reproduction crate's* numerical output,
//! not re-typed: the four Caltech-lineage solutions come straight from
//! `p9_core::types::P9Params` (the paper parameter sets the whole workspace
//! integrates), and the independent Siraj, Chyba & Tremaine (2024) solution
//! from `p9_2024_siraj_orbit::reference`. Only the per-element 1σ spreads are
//! added here, as documented assumptions standing in for each paper's stated
//! uncertainty (the same convention as REPRODUCTION_NOTES.md's
//! "tuned/assumed parameters").

use crate::schema::OrbitSolution;
use p9_core::constants::RAD2DEG;
use p9_core::types::P9Params;

/// Fiducial geometric albedo (Neptune V-band) used for every brightness
/// estimate unless a study motivates otherwise. See p9-core photometry.
pub const FIDUCIAL_ALBEDO: f64 = 0.40;

fn from_p9params(
    name: &str,
    citation: &str,
    arxiv: &str,
    p: &P9Params,
    (sa, se, si, som, sbom): (f64, f64, f64, f64, f64),
    note: &str,
) -> OrbitSolution {
    OrbitSolution {
        name: name.to_string(),
        citation: citation.to_string(),
        arxiv: arxiv.to_string(),
        mass_earth: p.mass_earth,
        albedo: FIDUCIAL_ALBEDO,
        a_au: p.a,
        a_sigma_au: sa,
        e: p.e,
        e_sigma: se,
        i_deg: p.i * RAD2DEG,
        i_sigma_deg: si,
        omega_deg: p.omega * RAD2DEG,
        omega_sigma_deg: som,
        omega_big_deg: p.omega_big * RAD2DEG,
        omega_big_sigma_deg: sbom,
        note: note.to_string(),
    }
}

/// The surveyed solutions, oldest to newest.
pub fn catalog() -> Vec<OrbitSolution> {
    let mut v = vec![
        from_p9params(
            "2016 Batygin & Brown (nominal)",
            "Batygin & Brown 2016, AJ 151, 22",
            "1601.05438",
            &P9Params::nominal_2016(),
            (150.0, 0.10, 10.0, 30.0, 30.0),
            "Discovery-paper nominal (10 M⊕, 700 AU, e=0.6, i=30°). Large 1σ \
             spreads reflect the broad 2016 constraints; ω, Ω from p9-core.",
        ),
        from_p9params(
            "2016 Batygin & Brown (inclined-TNO variant)",
            "Batygin & Brown 2016 (inclined-TNO study)",
            "1610.04992",
            &P9Params::inclined_tnos_2016(),
            (140.0, 0.10, 10.0, 30.0, 30.0),
            "Variant used for the high-inclination TNO production (10 M⊕, \
             600 AU, e=0.5, i=30°).",
        ),
        from_p9params(
            "2019 Batygin et al. (review best-fit)",
            "Batygin, Adams, Brown & Becker 2019, Phys. Rep. 805, 1",
            "1902.10103",
            &P9Params::revised_2019(),
            (120.0, 0.08, 8.0, 22.0, 22.0),
            "Review-paper revised parameters (5 M⊕, 500 AU, e=0.25, i=20°).",
        ),
        from_p9params(
            "2021 Brown & Batygin (MCMC median)",
            "Brown & Batygin 2021, AJ 162, 219",
            "2108.09868",
            &P9Params::mcmc_2021(),
            (80.0, 0.08, 5.0, 18.0, 18.0),
            "Statistical MCMC posterior median (6.2 M⊕, 380 AU, e=0.3, i=16°); \
             the current best estimate. ω, Ω from p9-core.",
        ),
    ];

    // Independent inference: Siraj, Chyba & Tremaine (2024). Mass, a, i and
    // their 1σ are read from the reproduction crate's published constants;
    // e and the apsidal/node angles are not inferred by that paper, so they
    // are taken anti-aligned with the ETNO cluster (same ω, Ω as the Caltech
    // solutions) with a moderate eccentricity — documented assumption.
    use p9_2024_siraj_orbit::reference as siraj;
    v.push(OrbitSolution {
        name: "2024 Siraj, Chyba & Tremaine (independent)".to_string(),
        citation: "Siraj, Chyba & Tremaine 2024".to_string(),
        arxiv: "2410.18170".to_string(),
        mass_earth: siraj::SIRAJ_2024_MASS_EARTH,
        albedo: FIDUCIAL_ALBEDO,
        a_au: siraj::SIRAJ_2024_A_AU,
        a_sigma_au: siraj::SIRAJ_2024_A_SIGMA_AU,
        e: 0.30,
        e_sigma: 0.10,
        i_deg: siraj::SIRAJ_2024_I_DEG,
        i_sigma_deg: siraj::SIRAJ_2024_I_SIGMA_DEG,
        omega_deg: 150.0,
        omega_sigma_deg: 25.0,
        omega_big_deg: 100.0,
        omega_big_sigma_deg: 25.0,
        note: "Mass/a/i and 1σ from p9-2024-siraj-orbit::reference \
               (4.4 M⊕, 290 AU, i=6.8°); markedly closer and lower-inclination \
               than Brown & Batygin 2021. e and ω/Ω assumed (anti-aligned)."
            .to_string(),
    });
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalog_has_expected_studies() {
        let c = catalog();
        assert_eq!(c.len(), 5);
        // Newest Caltech solution is the lightest/closest of that lineage.
        let mcmc = c.iter().find(|s| s.name.contains("2021")).unwrap();
        assert_eq!(mcmc.mass_earth, 6.2);
        assert_eq!(mcmc.a_au, 380.0);
    }

    #[test]
    fn siraj_is_sourced_from_reproduction_crate() {
        let c = catalog();
        let s = c.iter().find(|s| s.name.contains("Siraj")).unwrap();
        assert_eq!(s.mass_earth, 4.4);
        assert_eq!(s.a_au, 290.0);
        assert!((s.i_deg - 6.8).abs() < 1e-9);
    }

    #[test]
    fn all_solutions_physical() {
        for s in catalog() {
            assert!(s.mass_earth > 0.0 && s.a_au > 100.0);
            assert!((0.0..1.0).contains(&s.e));
            assert!(s.albedo > 0.0 && s.albedo <= 1.0);
            assert!(s.a_sigma_au > 0.0 && s.i_sigma_deg > 0.0);
        }
    }
}

//! Class tags and tier rules (design/02).
//!
//! The census cut and every class cut mirror the samples the paper crates
//! actually consume; the tier rules encode the a–e degeneracy honestly
//! (1-opposition / high-uncertainty orbits never enter frozen tables).

use crate::ledger::{Quality, Tier};

/// Census membership: is this object worth tracking at all?
///
/// `(a ≥ 100 AU) ∨ (i ≥ 60° ∧ a ≥ 30 AU ∧ q ≥ 15 AU)` — the first branch is
/// the distant population every analysis draws from; the second keeps the
/// Drac/Niku high-inclination family (a ≈ 35–45 AU) while excluding the
/// ~650 comet-like damocloids that a bare i-cut would haul in.
pub fn in_census(a: f64, q: f64, i_deg: f64) -> bool {
    a >= 100.0 || (i_deg >= 60.0 && a >= 30.0 && q >= 15.0)
}

/// Class tags for an orbit (design/02 table; an object can carry several).
/// The gap band matches the machinery it feeds
/// (`p9_2021_perihelion_gap::published::GAP_Q_LOW_AU..GAP_Q_HIGH_AU` =
/// 50..65 AU).
pub fn class_tags(a: f64, q: f64, i_deg: f64) -> Vec<String> {
    use p9_2021_perihelion_gap::published::{GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
    let mut tags = Vec::new();
    if a >= 250.0 && q >= 40.0 {
        tags.push("etno_vetted_class");
    }
    if a >= 230.0 && q > 30.0 && q < 40.0 {
        tags.push("etno_napier_class");
    }
    if (GAP_Q_LOW_AU..=GAP_Q_HIGH_AU).contains(&q) {
        tags.push("gap_band");
    }
    if q >= 40.0 && (100.0..250.0).contains(&a) {
        tags.push("gap_context");
    }
    if i_deg >= 60.0 {
        tags.push("high_inclination");
    }
    if i_deg > 90.0 {
        tags.push("retrograde");
    }
    if a > 100.0 && i_deg < 40.0 && q < 30.0 {
        tags.push("neptune_crossing_class");
    }
    if q >= 60.0 {
        tags.push("sedna_like");
    }
    if tags.is_empty() && a >= 150.0 {
        tags.push("watch_only");
    }
    tags.into_iter().map(String::from).collect()
}

/// Tier assignment (design/02): `vetted` requires multi-opposition, a long
/// arc, AND a known-good condition code; anything less certain — including
/// unknown quality — stays `provisional`.
pub fn tier(quality: &Quality, e: f64) -> Tier {
    if e >= 1.0 {
        return Tier::Provisional;
    }
    let multi_opp = quality.n_oppositions.is_some_and(|n| n >= 2);
    let long_arc = quality.arc_days.is_some_and(|d| d >= 400.0);
    let good_cc = quality.condition_code.is_some_and(|u| u <= 5);
    if multi_opp && long_arc && good_cc {
        Tier::Vetted
    } else {
        Tier::Provisional
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n_opp: u32, arc: f64, cc: u32) -> Quality {
        Quality {
            arc_days: Some(arc),
            n_oppositions: Some(n_opp),
            condition_code: Some(cc),
            n_obs: None,
        }
    }

    #[test]
    fn census_cut_keeps_drac_drops_damocloids_and_classicals() {
        assert!(in_census(41.7, 20.2, 103.4)); // Drac 2008 KV42
        assert!(in_census(523.0, 43.2, 12.6)); // 2025 LS2
        assert!(!in_census(44.1, 41.0, 2.2)); // classical KBO (Albion)
        assert!(!in_census(50.0, 2.0, 120.0)); // damocloid: retrograde, q ~ comet
    }

    #[test]
    fn class_table_matches_design_examples() {
        // Sedna: ETNO + sedna_like, not gap (q = 76).
        let t = class_tags(543.7, 76.2, 11.9);
        assert!(t.contains(&"etno_vetted_class".into()));
        assert!(t.contains(&"sedna_like".into()));
        assert!(!t.contains(&"gap_band".into()));
        // 2025 LS2: ETNO only (q = 43.2 below gap band, below sedna_like).
        let t = class_tags(523.0, 43.2, 12.6);
        assert_eq!(t, vec!["etno_vetted_class".to_string()]);
        // Gap-band object at a = 180: gap_band + gap_context.
        let t = class_tags(180.0, 55.0, 20.0);
        assert!(t.contains(&"gap_band".into()) && t.contains(&"gap_context".into()));
        // Neptune-crossing class.
        let t = class_tags(150.0, 25.0, 30.0);
        assert_eq!(t, vec!["neptune_crossing_class".to_string()]);
        // Distant nothing-special: watch_only.
        let t = class_tags(160.0, 35.0, 10.0);
        assert_eq!(t, vec!["watch_only".to_string()]);
        // Napier class.
        let t = class_tags(240.0, 33.0, 20.0);
        assert_eq!(t, vec!["etno_napier_class".to_string()]);
    }

    #[test]
    fn tier_rules_encode_the_degeneracy() {
        // Sedna-grade: vetted.
        assert_eq!(tier(&q(26, 13000.0, 5), 0.86), Tier::Vetted);
        // 2025 LS2 today: 2 oppositions but U = 9 → provisional.
        assert_eq!(tier(&q(2, 350.0, 9), 0.919), Tier::Provisional);
        // Multi-opp, long arc, unknown condition code → provisional.
        let unknown_cc = Quality {
            arc_days: Some(1000.0),
            n_oppositions: Some(3),
            condition_code: None,
            n_obs: None,
        };
        assert_eq!(tier(&unknown_cc, 0.5), Tier::Provisional);
        // Hyperbolic solution can never vet.
        assert_eq!(tier(&q(5, 2000.0, 1), 1.01), Tier::Provisional);
    }
}

//! The statistics battery (design/06): every number recomputed on a delta,
//! reusing the paper crates as libraries. Deterministic (seeded MCs).

use crate::ledger::{Ledger, ObjectRecord, Tier};
use p9_2026_apsidal_clustering::estimator::fit;
use p9_2026_apsidal_clustering::significance::sigma;
use p9_core::analysis::circular::{kuiper_p_value, mean_resultant_length, rayleigh_p_value};
use p9_core::analysis::selection::VarpiSelection;
use p9_core::constants::{DEG2RAD, RAD2DEG, TWO_PI};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct ClusterStats {
    pub n: usize,
    pub r_bar: f64,
    pub kappa: f64,
    pub mu_deg: f64,
    /// Gaussian-equivalent σ of the von Mises likelihood ratio (Wilks).
    pub sigma_gauss: f64,
    pub rayleigh_p: f64,
    pub kuiper_p: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct SelectionSpread {
    /// p under the cluster-aligned CosineLobes stand-in (napier-critique).
    pub p_cosine_lobes: f64,
    /// p under the CrossingDips stand-in (2021-orbit survey-bias null).
    pub p_crossing_dips: f64,
    /// p under the FlooredCrossingDips stand-in (2017-bias longitude factor).
    pub p_floored_dips: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct GapStats {
    /// Frozen paper-epoch sample baseline (p9-2021-perihelion-gap).
    pub p_paper_epoch: f64,
    /// Frozen current-epoch sample baseline.
    pub p_current_frozen: f64,
    /// Ledger objects in the 50–65 AU gap band not in the frozen tables —
    /// candidates for human curation into the gap sample.
    pub new_gap_band: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct BatteryResult {
    /// B1/B2 on the ledger's vetted-tier ETNO sample.
    pub vetted: Option<ClusterStats>,
    /// Same including provisional-tier objects.
    pub with_provisional: Option<ClusterStats>,
    /// Same as `with_provisional` but excluding this run's new objects
    /// (the "without" column; present only when the delta had new objects).
    pub without_new: Option<ClusterStats>,
    /// Frozen-crate baseline (p9-2026-apsidal-clustering vetted sample).
    pub baseline_frozen: ClusterStats,
    /// B3 on the ledger vetted sample.
    pub selection_spread: Option<SelectionSpread>,
    /// B4.
    pub gap: GapStats,
    /// B5: ledger Neptune-crossing objects absent from the frozen sample.
    pub neptune_crossing_new: Vec<String>,
    /// B6: ledger high-i/retrograde objects absent from the paper table.
    pub high_inclination_new: Vec<String>,
    /// B8.
    pub posterior_stale: bool,
    pub posterior_sample_diff: Vec<String>,
}

fn cluster_stats(varpi_rad: &[f64]) -> Option<ClusterStats> {
    if varpi_rad.len() < 3 {
        return None;
    }
    let f = fit(varpi_rad);
    Some(ClusterStats {
        n: f.n,
        r_bar: f.r_bar,
        kappa: f.kappa,
        mu_deg: f.mu.rem_euclid(TWO_PI) * RAD2DEG,
        sigma_gauss: sigma(f.lambda),
        rayleigh_p: rayleigh_p_value(varpi_rad),
        kuiper_p: kuiper_p_value(varpi_rad),
    })
}

/// ETNO ϖ sample from the ledger (radians) for the given tiers, e < 1 only.
fn etno_varpi(ledger: &Ledger, tiers: &[Tier], exclude_ids: &[String]) -> Vec<f64> {
    ledger
        .objects
        .iter()
        .filter(|o| {
            tiers.contains(&o.tier)
                && o.class_tags.iter().any(|t| t == "etno_vetted_class")
                && o.elements.e < 1.0
                && !exclude_ids.contains(&o.id)
        })
        .map(|o| o.elements.varpi_deg * DEG2RAD)
        .collect()
}

/// Selection-aware MC p (fraction of null resamples at least as concentrated
/// as observed), drawing ϖ from the density ∝ w(ϖ) by rejection. Mirrors the
/// core verdict-sensitivity machinery.
fn selection_p(sel: &VarpiSelection, varpi_obs: &[f64], seed: u64, iters: usize) -> f64 {
    let r_obs = mean_resultant_length(varpi_obs);
    let w_max = (0..3600)
        .map(|k| sel.weight(k as f64 * TWO_PI / 3600.0))
        .fold(0.0_f64, f64::max);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut hits = 0usize;
    for _ in 0..iters {
        let draw: Vec<f64> = (0..varpi_obs.len())
            .map(|_| loop {
                let v = rng.gen::<f64>() * TWO_PI;
                if rng.gen::<f64>() * w_max < sel.weight(v) {
                    break v;
                }
            })
            .collect();
        if mean_resultant_length(&draw) >= r_obs {
            hits += 1;
        }
    }
    hits as f64 / iters as f64
}

fn selection_spread(varpi: &[f64]) -> SelectionSpread {
    let lobe = VarpiSelection::CosineLobes {
        a1: 0.90,
        phi1: 52.0 * DEG2RAD,
        a2: 0.09,
        phi2: 52.0 * DEG2RAD,
    };
    let dips = VarpiSelection::CrossingDips {
        depth: 0.85,
        sigma: 20.0 * DEG2RAD,
    };
    let floored = VarpiSelection::FlooredCrossingDips {
        sigmas: [15.0 * DEG2RAD, 40.0 * DEG2RAD],
        floors: [0.5, 0.02],
    };
    SelectionSpread {
        p_cosine_lobes: selection_p(&lobe, varpi, 11, 4000),
        p_crossing_dips: selection_p(&dips, varpi, 13, 4000),
        p_floored_dips: selection_p(&floored, varpi, 17, 4000),
    }
}

/// Does this ledger object match a frozen-table name? Frozen names come in
/// "Sedna", "2012 VP113" and "Name (prov)" shapes; match against the alias
/// set and any parenthetical.
fn matches_frozen(obj: &ObjectRecord, frozen_name: &str) -> bool {
    let inner = frozen_name
        .rsplit_once('(')
        .map(|(_, t)| t.trim_end_matches(')').trim().to_string());
    let candidates: Vec<&str> = std::iter::once(frozen_name)
        .chain(inner.as_deref())
        .collect();
    candidates
        .iter()
        .any(|c| obj.primary_desig == *c || obj.aliases.iter().any(|a| a == c))
}

pub fn run(ledger: &Ledger, new_ids: &[String]) -> BatteryResult {
    use p9_2021_perihelion_gap::published::{GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
    use p9_2021_perihelion_gap::sample::{observed_perihelia, paper_epoch_perihelia};
    use p9_2021_perihelion_gap::synthetic::continuous_null_dip_p_value;

    let vetted_varpi = etno_varpi(ledger, &[Tier::Vetted], &[]);
    let both_varpi = etno_varpi(ledger, &[Tier::Vetted, Tier::Provisional], &[]);
    let without_new_varpi = if new_ids.is_empty() {
        Vec::new()
    } else {
        etno_varpi(ledger, &[Tier::Vetted, Tier::Provisional], new_ids)
    };

    // B4 baselines on the frozen crate samples (deterministic).
    let gap = GapStats {
        p_paper_epoch: continuous_null_dip_p_value(
            &paper_epoch_perihelia(),
            GAP_Q_LOW_AU,
            GAP_Q_HIGH_AU,
            12.0,
            20_000,
            2021,
        ),
        p_current_frozen: continuous_null_dip_p_value(
            &observed_perihelia(),
            GAP_Q_LOW_AU,
            GAP_Q_HIGH_AU,
            12.0,
            20_000,
            2021,
        ),
        new_gap_band: ledger
            .objects
            .iter()
            .filter(|o| o.tier != Tier::Retired && o.class_tags.iter().any(|t| t == "gap_band"))
            .map(|o| o.primary_desig.clone())
            .collect(),
    };

    // B5: Neptune-crossing membership diff vs the frozen sample.
    let frozen_nc = p9_2024_neptune_crossing::observed_tnos::observed_sample();
    let neptune_crossing_new = ledger
        .objects
        .iter()
        .filter(|o| {
            o.tier != Tier::Retired
                && o.class_tags.iter().any(|t| t == "neptune_crossing_class")
                && !frozen_nc.iter().any(|t| matches_frozen(o, t.name))
        })
        .map(|o| o.primary_desig.clone())
        .collect();

    // B6: high-inclination membership diff vs the paper table.
    let frozen_hi = p9_2016_inclined_tnos::known_objects::paper_tnos();
    let high_inclination_new = ledger
        .objects
        .iter()
        .filter(|o| {
            o.tier != Tier::Retired
                && o.class_tags.iter().any(|t| t == "high_inclination")
                && !frozen_hi
                    .iter()
                    .any(|t| matches_frozen(o, t.name) || matches_frozen(o, t.designation))
        })
        .map(|o| o.primary_desig.clone())
        .collect();

    // B8: vetted-tier ETNO sample vs BROWN_2017_SAMPLE.
    let brown = &p9_core::data::etno::BROWN_2017_SAMPLE;
    let ledger_etnos: Vec<&ObjectRecord> = ledger
        .objects
        .iter()
        .filter(|o| o.tier == Tier::Vetted && o.class_tags.iter().any(|t| t == "etno_vetted_class"))
        .collect();
    let mut diff: Vec<String> = Vec::new();
    for b in brown.iter() {
        if !ledger_etnos.iter().any(|o| matches_frozen(o, b.name)) {
            diff.push(format!("frozen-only: {}", b.name));
        }
    }
    for o in &ledger_etnos {
        if !brown.iter().any(|b| matches_frozen(o, b.name)) {
            diff.push(format!("ledger-only: {}", o.primary_desig));
        }
    }
    let posterior_stale = diff.len() >= 3;

    BatteryResult {
        vetted: cluster_stats(&vetted_varpi),
        with_provisional: cluster_stats(&both_varpi),
        without_new: if new_ids.is_empty() {
            None
        } else {
            cluster_stats(&without_new_varpi)
        },
        baseline_frozen: cluster_stats(&p9_2026_apsidal_clustering::samples::vetted_etno_varpi())
            .expect("frozen baseline sample nonempty"),
        selection_spread: if vetted_varpi.len() >= 3 {
            Some(selection_spread(&vetted_varpi))
        } else {
            None
        },
        gap,
        neptune_crossing_new,
        high_inclination_new,
        posterior_stale,
        posterior_sample_diff: diff,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ledger::{Elements, HistoryEvent, Quality};

    fn etno(id: &str, desig: &str, varpi: f64, tier: Tier) -> ObjectRecord {
        ObjectRecord {
            id: id.into(),
            primary_desig: desig.into(),
            aliases: vec![desig.into()],
            tier,
            class_tags: vec!["etno_vetted_class".into()],
            elements: Elements {
                epoch_jd: 0.0,
                a_au: 400.0,
                e: 0.85,
                q_au: 60.0,
                i_deg: 15.0,
                omega_deg: varpi,
                omega_big_deg: 0.0,
                varpi_deg: varpi,
                h_mag: None,
                source: "sbdb".into(),
                fetched_utc: "t".into(),
            },
            elements_alt: None,
            quality: Quality::default(),
            flags: vec![],
            history: vec![HistoryEvent {
                utc: "t".into(),
                event: "ingested".into(),
                detail: "test".into(),
            }],
        }
    }

    #[test]
    fn battery_runs_and_separates_tiers() {
        let mut ledger = Ledger::empty("t");
        // Tight vetted cluster at ~50 deg, one provisional outlier at 200.
        for (k, v) in [40.0, 50.0, 55.0, 60.0, 45.0].iter().enumerate() {
            ledger.objects.push(etno(
                &format!("rw-{:04}", k + 1),
                &format!("V{k}"),
                *v,
                Tier::Vetted,
            ));
        }
        ledger
            .objects
            .push(etno("rw-0006", "P0", 200.0, Tier::Provisional));
        let out = run(&ledger, &["rw-0006".to_string()]);
        let vetted = out.vetted.unwrap();
        let both = out.with_provisional.unwrap();
        assert_eq!(vetted.n, 5);
        assert_eq!(both.n, 6);
        assert!(vetted.r_bar > both.r_bar, "outlier must dilute R-bar");
        // without_new == vetted sample here (the provisional IS the new one).
        let without = out.without_new.unwrap();
        assert_eq!(without.n, 5);
        assert!((without.r_bar - vetted.r_bar).abs() < 1e-12);
        // Gap baselines reproduce the crate's documented values.
        assert!(out.gap.p_paper_epoch < 0.10);
        assert!(out.gap.p_current_frozen > out.gap.p_paper_epoch);
        // Selection spread exists and is ordered the way the core test pins.
        let s = out.selection_spread.unwrap();
        assert!(s.p_cosine_lobes > s.p_crossing_dips);
        // Synthetic 5-object sample shares nothing with BROWN_2017 -> stale.
        assert!(out.posterior_stale);
    }

    #[test]
    fn frozen_name_matching_handles_parentheticals() {
        let o = etno("rw-0001", "2004 VN112", 0.0, Tier::Vetted);
        assert!(matches_frozen(&o, "2004 VN112"));
        assert!(matches_frozen(&o, "Alicanto (2004 VN112)"));
        assert!(!matches_frozen(&o, "Sedna"));
    }
}

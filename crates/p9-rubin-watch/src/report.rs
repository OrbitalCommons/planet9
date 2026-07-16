//! Delta-report rendering (design/06 body structure, design/07 frontmatter).
//! The frontmatter duplicates every number machine-readably (as JSON, which
//! is valid YAML flow style) so downstream tooling never parses prose.

use crate::battery::{BatteryResult, ClusterStats};
use crate::ledger::Ledger;
use crate::sync::Delta;

pub fn render(
    run_utc: &str,
    ledger: &Ledger,
    delta: &Delta,
    battery: &BatteryResult,
    degraded: &[String],
) -> String {
    let mut s = String::new();
    let battery_json = serde_json::to_string(battery).expect("battery serializes");
    let delta_json = serde_json::to_string(delta).expect("delta serializes");
    s.push_str(&format!(
        "---\nschema_version: 1\nrun_utc: \"{run_utc}\"\nlayer: mpc\ntrigger: {delta_json}\ndegraded: {deg}\nbattery: {battery_json}\n---\n\n",
        deg = serde_json::to_string(degraded).unwrap(),
    ));

    s.push_str(&format!("# Rubin watch — MPC sync {run_utc}\n\n"));

    // 1. Trigger
    s.push_str("## Trigger\n\n");
    if delta.is_empty() {
        s.push_str("No census changes.\n\n");
    } else {
        for (label, list) in [
            ("New", &delta.new),
            ("Promoted", &delta.promoted),
            ("Demoted", &delta.demoted),
            ("Reclassified", &delta.reclassified),
            ("Renamed", &delta.renamed),
            ("Retired", &delta.retired),
            ("Discrepant (MPC vs SBDB)", &delta.discrepant),
        ] {
            if !list.is_empty() {
                s.push_str(&format!("- **{label}**: {}\n", list.join(", ")));
            }
        }
        s.push_str(&format!("- Elements refreshed: {}\n\n", delta.refreshed));
        const MAX_DETAIL: usize = 40;
        if delta.new.len() > MAX_DETAIL {
            s.push_str(&format!(
                "  (showing the first {MAX_DETAIL} of {} new objects; the ledger has all)\n",
                delta.new.len()
            ));
        }
        for desig in delta.new.iter().take(MAX_DETAIL) {
            if let Some(idx) = ledger.find_by_desig(desig) {
                let o = &ledger.objects[idx];
                s.push_str(&format!(
                    "  - {} — a = {:.1} AU, e = {:.3}, q = {:.1} AU, i = {:.1}°, ϖ = {:.1}°; \
                     tier {:?}, tags {:?}, arc {} d, opps {}, U {}\n",
                    o.primary_desig,
                    o.elements.a_au,
                    o.elements.e,
                    o.elements.q_au,
                    o.elements.i_deg,
                    o.elements.varpi_deg,
                    o.tier,
                    o.class_tags,
                    o.quality.arc_days.map_or("?".into(), |v| format!("{v:.0}")),
                    o.quality
                        .n_oppositions
                        .map_or("?".into(), |v| v.to_string()),
                    o.quality
                        .condition_code
                        .map_or("?".into(), |v| v.to_string()),
                ));
            }
        }
        s.push('\n');
    }

    // 2. Battery
    s.push_str("## Battery\n\n");
    s.push_str("### B1/B2 — ETNO apsidal clustering (ledger sample: a ≥ 250, q ≥ 40)\n\n");
    s.push_str("| sample | n | R̄ | κ | µ (°) | σ_gauss | Rayleigh p | Kuiper p |\n");
    s.push_str("|---|---|---|---|---|---|---|---|\n");
    push_cluster_row(&mut s, "frozen baseline", Some(&battery.baseline_frozen));
    push_cluster_row(&mut s, "ledger vetted", battery.vetted.as_ref());
    push_cluster_row(
        &mut s,
        "ledger vetted+prov",
        battery.with_provisional.as_ref(),
    );
    push_cluster_row(
        &mut s,
        "…without this run's new",
        battery.without_new.as_ref(),
    );
    s.push('\n');

    if let Some(sp) = &battery.selection_spread {
        s.push_str("### B3 — selection-conditioned clustering p (vetted sample)\n\n");
        s.push_str(&format!(
            "| CosineLobes (cluster-aligned) | CrossingDips | FlooredCrossingDips |\n|---|---|---|\n| {:.3} | {:.3} | {:.3} |\n\n",
            sp.p_cosine_lobes, sp.p_crossing_dips, sp.p_floored_dips
        ));
    }

    s.push_str("### B4 — perihelion gap (frozen baselines)\n\n");
    s.push_str(&format!(
        "paper-epoch p = {:.3}; current-epoch p = {:.3}. Ledger objects in the 50–65 AU band: {}\n\n",
        battery.gap.p_paper_epoch,
        battery.gap.p_current_frozen,
        if battery.gap.new_gap_band.is_empty() {
            "none".to_string()
        } else {
            battery.gap.new_gap_band.join(", ")
        }
    ));

    if !battery.neptune_crossing_new.is_empty() {
        s.push_str(&format!(
            "### B5 — Neptune-crossing candidates not in the frozen sample\n\n{}\n\n",
            battery.neptune_crossing_new.join(", ")
        ));
    }
    if !battery.high_inclination_new.is_empty() {
        s.push_str(&format!(
            "### B6 — high-inclination objects not in the paper table\n\n{}\n\n",
            battery.high_inclination_new.join(", ")
        ));
    }
    s.push_str(&format!(
        "### B8 — orbit-posterior staleness\n\nstale: **{}** (sample diff vs BROWN_2017: {})\n\n",
        battery.posterior_stale,
        if battery.posterior_sample_diff.is_empty() {
            "none".to_string()
        } else {
            battery.posterior_sample_diff.join("; ")
        }
    ));

    // 3. Guardrails
    s.push_str("## Interpretation guardrails\n\n");
    s.push_str(
        "- The clustering verdict is selection-model dependent (B3): quote the spread, \
         never one number.\n\
         - Provisional-tier orbits ride the a–e fit degeneracy; they never enter frozen \
         tables and their column is advisory.\n\
         - Gap-band candidates require human curation into the gap sample before the \
         B4 statistic moves (the frozen baselines above are unchanged by ingestion).\n\n",
    );

    // 4. Actions
    s.push_str("## Actions\n\n");
    if battery.posterior_stale {
        s.push_str(
            "- **posterior_stale**: open a re-derivation issue for the orbit-solution inputs.\n",
        );
    }
    if !battery.gap.new_gap_band.is_empty()
        || !battery.neptune_crossing_new.is_empty()
        || !battery.high_inclination_new.is_empty()
    {
        s.push_str(
            "- Frozen-table adoption candidates listed above (human-reviewed PR required).\n",
        );
    }
    if delta.is_empty() {
        s.push_str("- None.\n");
    }
    s
}

fn push_cluster_row(s: &mut String, label: &str, c: Option<&ClusterStats>) {
    match c {
        Some(c) => s.push_str(&format!(
            "| {label} | {} | {:.3} | {:.2} | {:.1} | {:.2} | {:.4} | {:.4} |\n",
            c.n, c.r_bar, c.kappa, c.mu_deg, c.sigma_gauss, c.rayleigh_p, c.kuiper_p
        )),
        None => s.push_str(&format!("| {label} | — | | | | | | |\n")),
    }
}

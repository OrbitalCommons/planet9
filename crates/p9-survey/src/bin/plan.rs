//! Evaluate observing strategies for a fixed JBT/SPENCER time budget.
//!
//! Usage: cargo run -p p9-survey --bin plan -- [--hours H] [--duty D] [--samples N]
//!
//! For several depth/area/cadence strategies it reports the coverable area, the
//! stacked depth, and the captured Planet Nine probability (≈ detection chance)
//! both for an equal-weight blend of the orbit solutions and for the current
//! best (2021 MCMC) and the tightest/closest (2024 Siraj) solutions.

use p9_survey::plan::{blended_prob, capture, stacked_depth, TimeBudget};
use p9_survey::{default_grid, run_survey};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut hours = 96.0; // 4 days
    let mut duty = 0.5;
    let mut samples = 300_000usize;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--hours" => hours = it.next().and_then(|v| v.parse().ok()).unwrap_or(hours),
            "--duty" => duty = it.next().and_then(|v| v.parse().ok()).unwrap_or(duty),
            "--samples" => samples = it.next().and_then(|v| v.parse().ok()).unwrap_or(samples),
            _ => {}
        }
    }

    let grid = default_grid();
    let d = run_survey(samples, 2026);
    let blend = blended_prob(&d.studies);
    let s2021 = d
        .studies
        .iter()
        .find(|s| s.solution.name.contains("2021 Brown"))
        .unwrap();
    let s2024 = d
        .studies
        .iter()
        .find(|s| s.solution.name.contains("Siraj"))
        .unwrap();

    // (label, exposure_s, overhead_s, epochs)
    let strategies = [
        ("wide-shallow  (120s ×3 epochs)", 120.0, 60.0, 3u32),
        ("balanced      (300s ×3 epochs)", 300.0, 60.0, 3),
        ("deep-narrow   (600s ×4 epochs)", 600.0, 90.0, 4),
    ];

    println!(
        "JBT 0.5m + SPENCER observing plan — {:.0} h wall ({:.1} days), duty {:.0}%, FoV 0.32 deg²\n",
        hours,
        hours / 24.0,
        duty * 100.0,
    );
    println!(
        "  {:32}  {:>7}  {:>6}  {:>6}   {:>8}  {:>8}  {:>8}",
        "strategy", "area", "depth", "Vmed?", "blend", "2021", "2024"
    );
    println!("  {}", "-".repeat(86));

    for (label, exp, oh, epochs) in strategies {
        let b = TimeBudget {
            hours_total: hours,
            duty_cycle: duty,
            exposure_s: exp,
            overhead_s: oh,
            epochs,
            fov_deg2: 0.32,
        };
        let area = b.area_deg2();
        let depth = stacked_depth(b.integration_per_field_s());

        // Depth-sufficiency factor per study: fraction of the V distribution
        // brighter than the stacked depth, estimated from the p16/med/p84 points.
        let eff = |vmed: f64, vp84: f64| -> f64 {
            if depth >= vp84 {
                1.0
            } else if depth >= vmed {
                0.75
            } else {
                0.4
            }
        };
        let cap_blend = capture(&blend, &grid, area).captured_prob;
        let c21 = capture(&s2021.prob, &grid, area);
        let c24 = capture(&s2024.prob, &grid, area);
        let det21 = c21.captured_prob * eff(s2021.v_median, s2021.v_p84);
        let det24 = c24.captured_prob * eff(s2024.v_median, s2024.v_p84);
        let det_blend = cap_blend * 0.85; // generic linking/gap efficiency

        let vmed_ok = if depth >= s2021.v_median {
            "ok"
        } else {
            "FAINT"
        };
        println!(
            "  {label:32}  {area:6.0}°  {depth:5.1}  {vmed_ok:>6}   {:7.1}%  {:7.1}%  {:7.1}%",
            100.0 * det_blend,
            100.0 * det21,
            100.0 * det24,
        );
    }

    // Recommended pointing for the wide-shallow strategy on the blended prior.
    let b = TimeBudget {
        hours_total: hours,
        duty_cycle: duty,
        exposure_s: 120.0,
        overhead_s: 60.0,
        epochs: 3,
        fov_deg2: 0.32,
    };
    let cov = capture(&blend, &grid, b.area_deg2());
    println!(
        "\n  Recommended field (wide-shallow, blended prior): {:.0} cells, {:.0} deg²\n  \
         RA {:.1}h–{:.1}h, Dec {:+.0}°…{:+.0}°, centroid RA {:.1}h Dec {:+.0}° — captures {:.1}% of the prior.",
        cov.n_cells as f64,
        cov.area_used_deg2,
        cov.ra_lo_deg / 15.0,
        cov.ra_hi_deg / 15.0,
        cov.dec_lo_deg,
        cov.dec_hi_deg,
        cov.ra_centroid_deg / 15.0,
        cov.dec_centroid_deg,
        100.0 * cov.captured_prob,
    );

    // Betting on the current best + independent modern solutions concentrates
    // the field at a brighter, tighter lobe and yields more.
    let cov21 = capture(&s2021.prob, &grid, b.area_deg2());
    let cov24 = capture(&s2024.prob, &grid, b.area_deg2());
    println!(
        "  If you back the modern orbit (recommended):\n  \
         · 2021 MCMC  → RA {:.1}h Dec {:+.0}°, captures {:.1}% (source V≈{:.1}, easy)\n  \
         · 2024 Siraj → RA {:.1}h Dec {:+.0}°, captures {:.1}% (source V≈{:.1}, easy)",
        cov21.ra_centroid_deg / 15.0,
        cov21.dec_centroid_deg,
        100.0 * cov21.captured_prob,
        s2021.v_median,
        cov24.ra_centroid_deg / 15.0,
        cov24.dec_centroid_deg,
        100.0 * cov24.captured_prob,
        s2024.v_median,
    );
}

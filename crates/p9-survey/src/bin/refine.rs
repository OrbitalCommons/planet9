//! Refine JBT/SPENCER's target: cede Rubin's southern sky and remove what ZTF,
//! Pan-STARRS and DES have already covered deep enough, then re-run the 4-day
//! observing plan on the residual.
//!
//! Usage: cargo run -p p9-survey --bin refine -- [--hours H] [--duty D]

use p9_survey::plan::{capture, TimeBudget};
use p9_survey::refine::{prior_surveys, refine, residual_mass, RUBIN_DEC_MAX_DEG};
use p9_survey::{default_grid, run_survey};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut hours = 96.0;
    let mut duty = 0.5;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--hours" => hours = it.next().and_then(|v| v.parse().ok()).unwrap_or(hours),
            "--duty" => duty = it.next().and_then(|v| v.parse().ok()).unwrap_or(duty),
            _ => {}
        }
    }

    let grid = default_grid();
    let d = run_survey(300_000, 2026);
    let surveys = prior_surveys();
    let budget = TimeBudget {
        hours_total: hours,
        duty_cycle: duty,
        exposure_s: 120.0,
        overhead_s: 60.0,
        epochs: 3,
        fov_deg2: 0.32,
    };
    let area = budget.area_deg2();

    println!(
        "Refined JBT target — cede Rubin (Dec ≤ {:+.0}°), survive ZTF/PS1/DES.\n\
         {hours:.0} h / duty {:.0}% → coverable area ≈ {area:.0} deg² (wide-shallow, depth ~V24).\n",
        RUBIN_DEC_MAX_DEG,
        duty * 100.0,
    );
    println!(
        "  {:30} {:>6}  {:>9}  {:>9}  {:>10}   {}",
        "solution (source V)", "rawcap", "residual", "JBT 4-day", "vs raw", "residual field"
    );
    println!("  {}", "-".repeat(98));

    for s in &d.studies {
        let v = s.v_median;
        let raw_cap = capture(&s.prob, &grid, area).captured_prob; // full prior, no cuts
        let resid = residual_mass(&s.prob, &grid, v, &surveys); // JBT's niche size
        let refined = refine(&s.prob, &grid, v, &surveys);
        let cov = capture(&refined, &grid, area); // 4-day reach on the residual
        let gain = if raw_cap > 0.0 {
            cov.captured_prob / raw_cap
        } else {
            0.0
        };
        let field = if cov.captured_prob > 1e-4 {
            format!(
                "RA {:.1}–{:.1}h Dec {:+.0}…{:+.0}°",
                cov.ra_lo_deg / 15.0,
                cov.ra_hi_deg / 15.0,
                cov.dec_lo_deg,
                cov.dec_hi_deg
            )
        } else {
            "(little left — others cover it)".to_string()
        };
        println!(
            "  {:30} {:5.1}%  {:8.1}%  {:8.1}%  {:9.2}×   {}",
            format!(
                "{} (V{:.1})",
                s.solution.name.split('(').next().unwrap().trim(),
                v
            ),
            100.0 * raw_cap,
            100.0 * resid,
            100.0 * cov.captured_prob,
            gain,
            field,
        );
    }

    println!(
        "\n  rawcap   = fraction of the full prior JBT could cover in 4 days (no cuts)\n  \
         residual = prior mass left for JBT after Rubin-cede + ZTF/PS1/DES survival (its niche)\n  \
         JBT 4-day= fraction actually caught in 4 days pointed at the residual\n  \
         Read-off: bright/close solutions → tiny residual (Rubin+ZTF+PS1 mop them up);\n  \
         faint/distant solutions → JBT owns the deep northern sky nobody else reaches."
    );
}

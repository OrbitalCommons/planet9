//! Refine JBT/SPENCER's target: cede Rubin's southern sky, remove what ZTF,
//! Pan-STARRS and DES have already covered (real footprints + plane depth
//! degradation), and fold in the Cassini/Iorio ephemeris favored-ν phase
//! prior. Then re-run the 4-day observing plan on the residual.
//!
//! Usage: cargo run -p p9-survey --bin refine -- [--hours H] [--duty D]

use p9_survey::ephemeris::nu_weight;
use p9_survey::plan::{capture, TimeBudget};
use p9_survey::refine::{prior_surveys, refine, RUBIN_DEC_MAX_DEG};
use p9_survey::skymap::ephemeris_weighted_grid;
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
    let (nu_lo, nu_hi) = (
        d.constraints.favored_nu_lo_deg,
        d.constraints.favored_nu_hi_deg,
    );
    let arc_dec_lo = d
        .constraints
        .favored_arc
        .iter()
        .map(|p| p[1])
        .fold(f64::MAX, f64::min);
    let arc_dec_hi = d
        .constraints
        .favored_arc
        .iter()
        .map(|p| p[1])
        .fold(f64::MIN, f64::max);

    println!(
        "Refined JBT target — cede Rubin (Dec ≤ {:+.0}°) + ZTF/PS1/DES survival + Cassini ν∈[{:.0},{:.0}]°.\n\
         {hours:.0} h / duty {:.0}% → coverable area ≈ {area:.0} deg² (wide-shallow, depth ~V24).\n\
         NOTE: the Cassini favored-ν arc lands at Dec {:+.0}…{:+.0}° — {} Rubin's +12° line, so the\n\
         ephemeris phase prior and the Rubin-cede pull against each other (see resid+eph below).\n",
        RUBIN_DEC_MAX_DEG, nu_lo, nu_hi, duty * 100.0, arc_dec_lo, arc_dec_hi,
        if arc_dec_hi <= RUBIN_DEC_MAX_DEG { "entirely south of" } else { "straddling" },
    );
    println!(
        "  {:26} {:>9} {:>9}   {:>8} {:>8}   {}",
        "solution (V)", "resid R+S", "resid+eph", "4d R+S", "4d+eph", "ephemeris-cut field"
    );
    println!("  {}", "-".repeat(96));

    for s in &d.studies {
        let v = s.v_median;

        // (a) Rubin-cede + survey survival on the uniform-phase prior.
        let refined_rs = refine(&s.prob, &grid, v, &surveys);
        let resid_rs: f64 = refined_rs.iter().sum();
        let cov_rs = capture(&refined_rs, &grid, area);

        // (b) Same, but on the ephemeris-phase-weighted prior.
        let ephem = ephemeris_weighted_grid(&s.solution, &grid, 200_000, 4242, &nu_weight);
        let refined_eph = refine(&ephem, &grid, v, &surveys);
        let resid_eph: f64 = refined_eph.iter().sum();
        let cov_eph = capture(&refined_eph, &grid, area);

        let field = if cov_eph.captured_prob > 1e-4 {
            format!(
                "RA {:.1}–{:.1}h Dec {:+.0}…{:+.0}°",
                cov_eph.ra_lo_deg / 15.0,
                cov_eph.ra_hi_deg / 15.0,
                cov_eph.dec_lo_deg,
                cov_eph.dec_hi_deg
            )
        } else {
            "(little left)".to_string()
        };
        let _ = cov_rs;
        println!(
            "  {:26} {:8.1}% {:8.1}%   {:7.1}% {:7.1}%   {}",
            format!(
                "{} V{:.1}",
                s.solution.name.split('(').next().unwrap().trim(),
                v
            ),
            100.0 * resid_rs,
            100.0 * resid_eph,
            100.0 * cov_rs.captured_prob,
            100.0 * cov_eph.captured_prob,
            field,
        );
    }

    println!(
        "\n  resid R+S  = prior mass left after Rubin-cede + ZTF/PS1/DES survival\n  \
         resid+eph  = …and after the Cassini favored-ν phase prior\n  \
         4d R+S / 4d+eph = fraction caught in 4 days, without / with the ν cut\n  \
         Result: the Cassini favored-ν zone sits SOUTH of +12° — in Rubin's sky — so the\n  \
         ν prior and the Rubin-cede oppose each other: folding it in *shrinks* JBT's residual\n  \
         (4d+eph < 4d R+S). Read another way, the ν prior says 'if Cassini is right, let Rubin\n  \
         have it'; JBT's northern niche is the case where the (weak) Cassini fit does NOT hold."
    );
}

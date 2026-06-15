//! Generate the Planet Nine survey dataset as JSON for the Python plotter.
//!
//! Usage:
//!   cargo run -p p9-survey --bin survey -- [--out PATH] [--samples N] [--seed S]
//!
//! Defaults: --out p9_survey_data.json, --samples 300000, --seed 2026.

use std::fs;
use std::process::ExitCode;

fn main() -> ExitCode {
    let mut out = "p9_survey_data.json".to_string();
    let mut samples: usize = 300_000;
    let mut seed: u64 = 2026;

    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--out" => match it.next() {
                Some(v) => out = v.clone(),
                None => return fail("--out needs a path"),
            },
            "--samples" => match it.next().and_then(|v| v.parse().ok()) {
                Some(v) => samples = v,
                None => return fail("--samples needs a positive integer"),
            },
            "--seed" => match it.next().and_then(|v| v.parse().ok()) {
                Some(v) => seed = v,
                None => return fail("--seed needs an integer"),
            },
            "-h" | "--help" => {
                println!(
                    "survey [--out PATH] [--samples N] [--seed S]\n\
                     writes a p9-survey SurveyDataset as JSON"
                );
                return ExitCode::SUCCESS;
            }
            other => return fail(&format!("unknown argument: {other}")),
        }
    }

    eprintln!("p9-survey: {samples} samples/study, seed {seed} -> {out}");
    let dataset = p9_survey::run_survey(samples, seed);

    // Compact human-readable summary to stderr (the JSON is the real output).
    eprintln!("\n  study                                          peak RA    Dec    V(med)  d(med)   95% area");
    for s in &dataset.studies {
        eprintln!(
            "  {:44}  {:5.1}h  {:+5.1}°  {:5.1}   {:4.0}AU  {:7.0} deg²",
            s.solution.name,
            s.peak_ra_deg / 15.0,
            s.peak_dec_deg,
            s.v_median,
            s.dist_median_au,
            s.area95_deg2,
        );
    }
    eprintln!("\n  detection probability (in footprint AND bright enough):");
    for t in &dataset.telescopes {
        let best = t
            .per_study
            .iter()
            .map(|d| d.detection_prob)
            .fold(0.0_f64, f64::max);
        eprintln!(
            "  {:36}  depth {:4.1}  best P(detect) over studies = {:4.1}%",
            t.telescope.name,
            t.telescope.limiting_mag,
            100.0 * best
        );
    }

    // Parallax channel: LEO (per-orbit) vs Earth's annual baseline, per study.
    use p9_core::coords::parallax::{annual_parallax_arcsec, leo_parallax_arcsec, LEO_ALTITUDE_KM};
    use p9_survey::telescope::JBT_PLATE_SCALE_ARCSEC_PER_PX as PXSCALE;
    eprintln!(
        "\n  parallax of P9 at its median distance — LEO ({:.0} km, per-orbit) vs annual (Earth 2 AU):",
        LEO_ALTITUDE_KM
    );
    for s in &dataset.studies {
        let d = s.dist_median_au;
        let leo = leo_parallax_arcsec(LEO_ALTITUDE_KM, d);
        let ann = annual_parallax_arcsec(d);
        eprintln!(
            "  {:44}  d={:4.0}AU  LEO {:.3}″ ({:.2} px)  | annual {:5.0}″ ({:.1}′)  | annual/LEO {:.0}×",
            s.solution.name,
            d,
            leo,
            leo / PXSCALE,
            ann,
            ann / 60.0,
            ann / leo,
        );
    }

    let json = match serde_json::to_string_pretty(&dataset) {
        Ok(j) => j,
        Err(e) => return fail(&format!("serialize: {e}")),
    };
    if let Err(e) = fs::write(&out, json) {
        return fail(&format!("write {out}: {e}"));
    }
    eprintln!("\nwrote {out}");
    ExitCode::SUCCESS
}

fn fail(msg: &str) -> ExitCode {
    eprintln!("error: {msg}");
    ExitCode::FAILURE
}

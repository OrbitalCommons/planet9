//! Film export for `p9-2025-new-discoveries`: the numbers its scene and ledger entry draw.

use p9_2025_new_discoveries::discoveries::{NewDiscovery, ammonite, of201};
use p9_2025_new_discoveries::stress::{ClusteringStat, run_stress};
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::{GM_SUN, RAD2DEG, TWO_PI};
use p9_core::data::etno::{BROWN_2017_SAMPLE, Etno};
use p9_core::types::{OrbitalElements, elements_to_cartesian, true_to_mean_anomaly};
use serde_json::{Value, json};

/// Points along each drawn orbit.
const N_PATH: usize = 160;
/// Objects in the Brown & Batygin (2021) clustering sample: the state of the
/// field that these discoveries extend.
const FIELD_SAMPLE_BEFORE: usize = 11;

/// The orbit projected on the ecliptic plane, as (x, y) in AU.
fn orbit_path(elements: &OrbitalElements) -> Vec<(f64, f64)> {
    (0..=N_PATH)
        .map(|k| {
            let nu = TWO_PI * k as f64 / N_PATH as f64;
            let at = OrbitalElements {
                mean_anomaly: true_to_mean_anomaly(elements.e, nu),
                ..*elements
            };
            let pos = elements_to_cartesian(&at, GM_SUN).pos;
            (pos.x, pos.y)
        })
        .collect()
}

fn etno_json(o: &Etno) -> Value {
    json!({
        "name": o.name,
        "a": o.a,
        "e": o.e,
        "q": o.perihelion(),
        "aphelion": o.a * (1.0 + o.e),
        "i_deg": o.i_deg,
        "h_mag": o.h_mag,
        "varpi_deg": o.longitude_of_perihelion() * RAD2DEG,
        "path": orbit_path(&o.elements()),
    })
}

fn discovery_json(d: &NewDiscovery, offset_deg: f64) -> Value {
    let mut v = etno_json(&d.etno);
    let extra = json!({
        "arxiv": d.arxiv,
        "paper_a": d.paper_a_au,
        "paper_q": d.paper_q_au,
        "offset_deg": offset_deg,
    });
    v.as_object_mut()
        .unwrap()
        .extend(extra.as_object().unwrap().clone());
    v
}

fn stat_json(label: &str, s: &ClusteringStat) -> Value {
    json!({
        "label": label,
        "n": s.n,
        "r_bar": s.r_bar,
        "rayleigh_p": s.rayleigh_p,
        "sigma": p_value_to_sigma(s.rayleigh_p),
        "mean_varpi_deg": s.mean_varpi_deg.rem_euclid(360.0),
    })
}

pub fn export() -> Value {
    let (of, am) = (of201(), ammonite());
    let stress = run_stress(&of, &am);

    json!({
        "baseline": BROWN_2017_SAMPLE.iter().map(etno_json).collect::<Vec<_>>(),
        "of201": discovery_json(&of, stress.of201_offset_deg),
        "ammonite": discovery_json(&am, stress.ammonite_offset_deg),
        "stress": [
            stat_json("before", &stress.baseline),
            stat_json("+ 2017 OF201", &stress.with_of201),
            stat_json("+ Ammonite", &stress.with_ammonite),
            stat_json("+ both", &stress.with_both),
        ],
        // This paper adds 2017 OF201 only; Ammonite is a later paper.
        "n_field": FIELD_SAMPLE_BEFORE + (stress.with_of201.n - stress.baseline.n),
        "of201_a": of.etno.a,
    })
}

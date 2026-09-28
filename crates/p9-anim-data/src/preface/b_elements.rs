//! Film export for `scenes/preface/preface_b_elements.py`: the numbers its scenes draw.
//!
//! * `demo`: the real elements of 2012 VP113, the orbit the elements lesson builds.
//! * `geography`: every orbit on the outer-system zoom (planets, Pluto, the
//!   distant objects with full elements, the published Planet Nine orbits),
//!   the (a, q) map of the distant objects, and light-time anchors.

use std::collections::BTreeMap;

use serde_json::{Value, json};

use p9_core::analysis::thermal::{AU_M, C_LIGHT};
use p9_core::constants::{A_NEPTUNE_AU, DEG2RAD};
use p9_core::data::etno::{BROWN_2017_SAMPLE, Etno};
use p9_core::data::stable_kbos::stable_kbos;
use p9_core::initial_conditions::giant_planets::GIANT_PLANETS;
use p9_core::types::{OrbitalElements, P9Params};

use p9_2021_des_catalog::catalog::DES_EXTREME_TNOS;
use p9_2021_detached_inclinations::reference::DETACHED_Q_MIN_AU;
use p9_2021_perihelion_gap::published::{GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
use p9_2021_perihelion_gap::sample::EXTENDED_DISTANT_TNOS;
use p9_2025_new_discoveries::{ammonite, of201};

use super::a_foundations::{
    PLUTO_A_AU, PLUTO_ARGP_DEG, PLUTO_E, PLUTO_I_DEG, PLUTO_NODE_DEG, orbit_xyz, p9_elements,
    period_yr,
};

/// Classical Kuiper belt: from Neptune out to the ~50 AU edge (the 2:1
/// resonance with Neptune sits at 47.8 AU).
const KUIPER_OUTER_AU: f64 = 50.0;

/// Voyager 1: launched 1977-09-05; 162 AU from the Sun on 2024-01-01,
/// receding at 3.57 AU/yr (JPL).
const VOYAGER1_LAUNCH_YEAR: f64 = 1977.68;
const VOYAGER1_AU_2024: f64 = 162.0;
const VOYAGER1_AU_PER_YR: f64 = 3.57;
const FILM_YEAR: f64 = 2026.75;

fn demo() -> Value {
    let o = BROWN_2017_SAMPLE
        .iter()
        .find(|o| o.name == "2012 VP113")
        .expect("2012 VP113 in the Brown (2017) sample");
    json!({
        "name": o.name,
        "a": o.a,
        "e": o.e,
        "i_deg": o.i_deg,
        "node_deg": o.omega_big_deg,
        "argp_deg": o.omega_deg,
        "varpi_deg": o.longitude_of_perihelion().to_degrees(),
        "q": o.perihelion(),
        "big_q": o.a * (1.0 + o.e),
        "period_yr": period_yr(o.a),
    })
}

fn object(name: &str, el: &OrbitalElements, source: &str) -> Value {
    json!({
        "name": name,
        "a": el.a,
        "e": el.e,
        "q": el.a * (1.0 - el.e),
        "i_deg": el.i.to_degrees(),
        "source": source,
        "xyz": orbit_xyz(el, 200),
    })
}

fn geography() -> Value {
    // Distant objects with full elements, one entry per name (first source wins).
    let mut full: BTreeMap<String, Value> = BTreeMap::new();
    let etno = |o: &Etno, src: &str| (o.name.to_string(), object(o.name, &o.elements(), src));
    for o in BROWN_2017_SAMPLE.iter() {
        let (k, v) = etno(o, "Brown 2017");
        full.entry(k).or_insert(v);
    }
    for k in stable_kbos() {
        full.entry(k.name.to_string())
            .or_insert_with(|| object(k.name, &k.elements, "Batygin & Brown 2016"));
    }
    for t in DES_EXTREME_TNOS.iter() {
        let el = OrbitalElements {
            a: t.a,
            e: t.e,
            i: t.i_deg * DEG2RAD,
            omega_big: t.node_deg * DEG2RAD,
            omega: t.argp_deg * DEG2RAD,
            mean_anomaly: 0.0,
        };
        full.entry(t.name.to_string())
            .or_insert_with(|| object(t.name, &el, "DES"));
    }
    for d in [of201(), ammonite()] {
        let (k, v) = etno(&d.etno, "2025 discoveries");
        full.entry(k).or_insert(v);
    }
    // Objects known only by (a, e): they join the (a, q) map, not the zoom.
    let mut aq: Vec<Value> = full
        .values()
        .map(|v| json!({"name": v["name"], "a": v["a"], "q": v["q"]}))
        .collect();
    for t in EXTENDED_DISTANT_TNOS.iter() {
        if !full.contains_key(t.name) {
            aq.push(json!({"name": t.name, "a": t.a, "q": t.a * (1.0 - t.e)}));
        }
    }

    let names = ["Jupiter", "Saturn", "Uranus", "Neptune"];
    let mut planets = vec![json!({"name": "Earth", "a": 1.0})];
    for (name, &(_, a)) in names.iter().zip(GIANT_PLANETS.iter()) {
        planets.push(json!({"name": name, "a": a}));
    }
    let pluto = OrbitalElements {
        a: PLUTO_A_AU,
        e: PLUTO_E,
        i: PLUTO_I_DEG * DEG2RAD,
        omega_big: PLUTO_NODE_DEG * DEG2RAD,
        omega: PLUTO_ARGP_DEG * DEG2RAD,
        mean_anomaly: 0.0,
    };
    let p9: Vec<Value> = [
        ("Batygin & Brown 2016", P9Params::nominal_2016()),
        ("Batygin et al. 2019", P9Params::revised_2019()),
        ("Brown & Batygin 2021", P9Params::mcmc_2021()),
    ]
    .iter()
    .map(|(label, p)| {
        json!({
            "label": label,
            "a": p.a,
            "e": p.e,
            "q": p.a * (1.0 - p.e),
            "mass_earth": p.mass_earth,
            "xyz": orbit_xyz(&p9_elements(p), 200),
        })
    })
    .collect();

    let light_hr_per_au = AU_M / C_LIGHT / 3600.0;
    let p9_ref = P9Params::revised_2019();
    json!({
        "planets": planets,
        "pluto": object("Pluto", &pluto, "JPL"),
        "kuiper_belt_au": [A_NEPTUNE_AU, KUIPER_OUTER_AU],
        "distant": full.into_values().collect::<Vec<_>>(),
        "aq": aq,
        "p9": p9,
        "detached_q_min": DETACHED_Q_MIN_AU,
        "gap_q": [GAP_Q_LOW_AU, GAP_Q_HIGH_AU],
        "neptune_a": A_NEPTUNE_AU,
        "light_hr_per_au": light_hr_per_au,
        "light_hr_neptune": A_NEPTUNE_AU * light_hr_per_au,
        "light_day_p9": p9_ref.a * light_hr_per_au / 24.0,
        "sunlight_p9_vs_neptune": (A_NEPTUNE_AU / p9_ref.a).powi(2),
        "voyager1_au": VOYAGER1_AU_2024 + VOYAGER1_AU_PER_YR * (FILM_YEAR - 2024.0),
        "voyager1_years_flown": FILM_YEAR - VOYAGER1_LAUNCH_YEAR,
        "voyager1_yr_to_p9": p9_ref.a / VOYAGER1_AU_PER_YR,
    })
}

pub fn export() -> Value {
    json!({
        "demo": demo(),
        "geography": geography(),
    })
}

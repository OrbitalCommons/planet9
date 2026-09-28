//! Film export for `scenes/preface/preface_a_foundations.py`: the numbers its scenes draw.
//!
//! * `hook`: the real Brown (2017) ETNO orbits projected on the ecliptic, their
//!   longitudes of perihelion, and a Planet Nine orbit for the opening image.
//! * `launch`: bodies launched sideways from Neptune's distance at several
//!   speeds, propagated with the core universal-variable Kepler solver.
//! * `sedna`: Sedna's shape, and its distance over one full orbit in time.
//! * `kepler3`: periods of real bodies and of the published Planet Nine orbits.

use nalgebra::Vector3;
use serde_json::{Value, json};

use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::analysis::elements::mean_motion;
use p9_core::constants::{A_NEPTUNE_AU, GM_SUN, KMS_TO_AUDAY, TWO_PI, YEAR_DAYS};
use p9_core::data::etno::{BROWN_2017_SAMPLE, Etno};
use p9_core::initial_conditions::giant_planets::GIANT_PLANETS;
use p9_core::integrator::kepler_step::kepler_drift;
use p9_core::types::{OrbitalElements, P9Params, StateVector, cartesian_to_elements};

/// Pluto's heliocentric J2000 elements (JPL): the familiar Kuiper-belt anchor.
pub const PLUTO_A_AU: f64 = 39.48;
pub const PLUTO_E: f64 = 0.2488;
pub const PLUTO_I_DEG: f64 = 17.14;
pub const PLUTO_NODE_DEG: f64 = 110.30;
pub const PLUTO_ARGP_DEG: f64 = 113.76;

/// Publication of Newton's Principia, and the year Neptune was found.
const PRINCIPIA_YEAR: f64 = 1687.0;
const NEPTUNE_FOUND_YEAR: f64 = 1846.0;
const FILM_YEAR: f64 = 2026.0;

fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

/// Orbital period in years for a heliocentric semi-major axis in AU.
pub fn period_yr(a_au: f64) -> f64 {
    TWO_PI / mean_motion(a_au) / YEAR_DAYS
}

/// Heliocentric ecliptic positions around a full orbit, sampled uniformly in
/// eccentric anomaly (dense at both apsides, so high-e ellipses stay smooth).
pub fn orbit_xyz(el: &OrbitalElements, n: usize) -> Vec<[f64; 3]> {
    (0..=n)
        .map(|k| {
            let ea = TWO_PI * k as f64 / n as f64;
            let el = OrbitalElements {
                mean_anomaly: ea - el.e * ea.sin(),
                ..*el
            };
            let p = el.to_state_vector(GM_SUN).pos;
            [round2(p.x), round2(p.y), round2(p.z)]
        })
        .collect()
}

pub fn p9_elements(p: &P9Params) -> OrbitalElements {
    OrbitalElements {
        a: p.a,
        e: p.e,
        i: p.i,
        omega_big: p.omega_big,
        omega: p.omega,
        mean_anomaly: 0.0,
    }
}

fn etno_json(o: &Etno) -> Value {
    json!({
        "name": o.name,
        "a": o.a,
        "e": o.e,
        "q": o.perihelion(),
        "big_q": o.a * (1.0 + o.e),
        "i_deg": o.i_deg,
        "varpi_deg": o.longitude_of_perihelion().to_degrees(),
        "xyz": orbit_xyz(&o.elements(), 160),
    })
}

fn hook() -> Value {
    let varpis: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| o.longitude_of_perihelion())
        .collect();
    let p9 = P9Params::revised_2019();
    json!({
        "etnos": BROWN_2017_SAMPLE.iter().map(etno_json).collect::<Vec<_>>(),
        "r_bar": mean_resultant_length(&varpis),
        "mean_varpi_deg": circular_mean(&varpis).unwrap_or(0.0).to_degrees().rem_euclid(360.0),
        "neptune_a": A_NEPTUNE_AU,
        "p9": {
            "a": p9.a,
            "e": p9.e,
            "mass_earth": p9.mass_earth,
            "varpi_deg": (p9.omega + p9.omega_big).to_degrees().rem_euclid(360.0),
            "xyz": orbit_xyz(&p9_elements(&p9), 160),
        },
    })
}

/// Launch a body sideways from Neptune's distance at `f` times the circular
/// speed and follow it: one full period when bound, out to 260 AU when not.
fn launch_one(f: f64) -> Value {
    let r0 = A_NEPTUNE_AU;
    let vc = (GM_SUN / r0).sqrt();
    let s0 = StateVector::new(Vector3::new(r0, 0.0, 0.0), Vector3::new(0.0, f * vc, 0.0));
    let el = cartesian_to_elements(&s0, GM_SUN);
    let bound = f * f < 2.0;
    let n = 360;
    let t_end = if bound {
        TWO_PI * (el.a.powi(3) / GM_SUN).sqrt()
    } else {
        // step outward until the body passes 260 AU
        let mut t = 0.0;
        while kepler_drift(&s0, t, GM_SUN).pos.norm() < 260.0 {
            t += 365.25;
        }
        t
    };
    let xy: Vec<[f64; 2]> = (0..=n)
        .map(|k| {
            let s = kepler_drift(&s0, t_end * k as f64 / n as f64, GM_SUN);
            [round2(s.pos.x), round2(s.pos.y)]
        })
        .collect();
    json!({
        "factor": f,
        "speed_kms": f * vc / KMS_TO_AUDAY,
        "bound": bound,
        "a": el.a,
        "e": el.e,
        "q": el.a * (1.0 - el.e),
        "big_q": if bound { json!(el.a * (1.0 + el.e)) } else { Value::Null },
        "duration_yr": t_end / YEAR_DAYS,
        "xy": xy,
    })
}

fn launch() -> Value {
    let r0 = A_NEPTUNE_AU;
    let vc = (GM_SUN / r0).sqrt() / KMS_TO_AUDAY;
    let paths: Vec<Value> = [0.6, 0.8, 1.0, 1.2, 1.5]
        .iter()
        .map(|&f| launch_one(f))
        .collect();
    json!({
        "r0_au": r0,
        "v_circ_kms": vc,
        "v_esc_kms": vc * 2f64.sqrt(),
        "paths": paths,
    })
}

fn sedna() -> Value {
    let s = BROWN_2017_SAMPLE
        .iter()
        .find(|o| o.name == "Sedna")
        .expect("Sedna in the Brown (2017) sample");
    let el = s.elements();
    let period_days = TWO_PI / mean_motion(s.a);
    let s0 = el.to_state_vector(GM_SUN); // mean anomaly 0: at perihelion
    let n = 600;
    let (mut t_yr, mut r_au) = (Vec::new(), Vec::new());
    let mut inside_100 = 0usize;
    for k in 0..=n {
        let t = period_days * k as f64 / n as f64;
        let r = kepler_drift(&s0, t, GM_SUN).pos.norm();
        if k < n && r < 100.0 {
            inside_100 += 1;
        }
        t_yr.push(round2(t / YEAR_DAYS));
        r_au.push(round2(r));
    }
    let frac = inside_100 as f64 / n as f64;
    json!({
        "a": s.a,
        "e": s.e,
        "q": s.perihelion(),
        "big_q": s.a * (1.0 + s.e),
        "period_yr": period_days / YEAR_DAYS,
        "t_yr": t_yr,
        "r_au": r_au,
        "frac_inside_100": frac,
        "yr_inside_100": frac * period_days / YEAR_DAYS,
    })
}

fn kepler3() -> Value {
    let names = ["Jupiter", "Saturn", "Uranus", "Neptune"];
    let mut bodies = vec![json!({"name": "Earth", "a": 1.0, "period_yr": period_yr(1.0)})];
    for (name, &(_, a)) in names.iter().zip(GIANT_PLANETS.iter()) {
        bodies.push(json!({"name": name, "a": a, "period_yr": period_yr(a)}));
    }
    bodies.push(json!({"name": "Pluto", "a": PLUTO_A_AU, "period_yr": period_yr(PLUTO_A_AU)}));

    let estimates = [
        ("Batygin & Brown 2016", P9Params::nominal_2016()),
        ("Batygin et al. 2019", P9Params::revised_2019()),
        ("Brown & Batygin 2021", P9Params::mcmc_2021()),
    ];
    let p9: Vec<Value> = estimates
        .iter()
        .map(|(label, p)| json!({"label": label, "a": p.a, "period_yr": period_yr(p.a)}))
        .collect();
    let p_ref = period_yr(P9Params::revised_2019().a);
    json!({
        "bodies": bodies,
        "p9": p9,
        "years_since_principia": FILM_YEAR - PRINCIPIA_YEAR,
        "p9_deg_since_principia": 360.0 * (FILM_YEAR - PRINCIPIA_YEAR) / p_ref,
        "neptune_orbits_since_found": (FILM_YEAR - NEPTUNE_FOUND_YEAR) / period_yr(A_NEPTUNE_AU),
    })
}

pub fn export() -> Value {
    json!({
        "hook": hook(),
        "launch": launch(),
        "sedna": sedna(),
        "kepler3": kepler3(),
    })
}

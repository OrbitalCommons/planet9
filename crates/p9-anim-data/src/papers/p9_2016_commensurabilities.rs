//! Film export for `p9-2016-commensurabilities`: the numbers its scene and ledger entry draw.
//!
//! The ledger entry is de la Fuente Marcos, de la Fuente Marcos & Aarseth
//! (2016), "Dynamical impact of the Planet Nine scenario: N-body experiments":
//! the six objects first linked to the hypothesis are integrated under the
//! nominal Planet Nine to see whether their orbits stay confined. The
//! reproduction crate holds the grouping statistics; the integration itself is
//! run here with the `p9-core` hybrid integrator at reduced scale (one clone
//! per object, tens of Myr).

use std::thread;

use p9_2016_commensurabilities::clustering::grouping_of;
use p9_core::constants::{DEG2RAD, GM_SUN, RAD2DEG, YEAR_DAYS};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets::neptune_j2000;
use p9_core::integrator::hybrid::hybrid_step_with_forces;
use p9_core::types::{P9Params, SimConfig, cartesian_to_elements};
use serde_json::{Value, json};

/// The six objects of Batygin & Brown (2016) that the paper integrates.
const SIX: [&str; 6] = [
    "Sedna",
    "2012 VP113",
    "2010 GB174",
    "2004 VN112",
    "2013 RF98",
    "2007 TG422",
];

/// Reduced-scale run length and sampling.
const T_MYR: f64 = 30.0;
const DT_DAYS: f64 = 3000.0;
const SNAPSHOT_MYR: f64 = 0.2;

/// An orbit counts as unstable once its semi-major axis has wandered by more
/// than this fraction of its starting value (or it has been removed).
const UNSTABLE_EXCURSION: f64 = 0.05;

/// The nominal Planet Nine of Batygin & Brown (2016), placed at aphelion.
fn nominal_planet() -> P9Params {
    P9Params {
        mass_earth: 10.0,
        a: 700.0,
        e: 0.6,
        i: 30.0 * DEG2RAD,
        omega: 150.0 * DEG2RAD,
        omega_big: 113.0 * DEG2RAD,
        mean_anomaly: 180.0 * DEG2RAD,
    }
}

struct Track {
    name: &'static str,
    t_myr: Vec<f64>,
    a_au: Vec<f64>,
    q_au: Vec<f64>,
    i_deg: Vec<f64>,
    omega_rad: Vec<f64>,
    /// Time the object was removed from the run or became unbound, if it did.
    lost_myr: Option<f64>,
    /// Longitude of perihelion and eccentricity at the start (for the top view).
    varpi0_deg: f64,
    e0: f64,
}

impl Track {
    /// Largest fractional departure of the semi-major axis from its start.
    fn excursion(&self) -> f64 {
        let a0 = self.a_au[0];
        self.a_au
            .iter()
            .map(|a| (a - a0).abs() / a0)
            .fold(0.0, f64::max)
    }

    /// First time the semi-major axis excursion passed the threshold.
    fn unstable_myr(&self) -> Option<f64> {
        let a0 = self.a_au[0];
        self.a_au
            .iter()
            .position(|a| (a - a0).abs() / a0 > UNSTABLE_EXCURSION)
            .map(|k| self.t_myr[k])
            .or(self.lost_myr)
    }
}

fn integrate(name: &'static str) -> Track {
    let etno = BROWN_2017_SAMPLE
        .iter()
        .find(|e| e.name == name)
        .unwrap_or_else(|| panic!("{name} missing from BROWN_2017_SAMPLE"));
    let mut particles = vec![etno.elements().to_state_vector(GM_SUN)];
    let mut active = vec![true];
    let mut bodies = vec![neptune_j2000(), nominal_planet().to_body()];
    let config = SimConfig {
        dt: DT_DAYS,
        t_start: 0.0,
        t_end: T_MYR * 1e6 * YEAR_DAYS,
        removal_inner_au: 5.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: SNAPSHOT_MYR * 1e6 * YEAR_DAYS,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };
    let extra = [ExtraForce::J2Jsu];
    let n_steps = (config.t_end / DT_DAYS).ceil() as usize;
    let snap_every = (config.snapshot_interval_days / DT_DAYS).ceil() as usize;

    let mut track = Track {
        name,
        t_myr: Vec::new(),
        a_au: Vec::new(),
        q_au: Vec::new(),
        i_deg: Vec::new(),
        omega_rad: Vec::new(),
        lost_myr: None,
        varpi0_deg: etno.longitude_of_perihelion() * RAD2DEG,
        e0: etno.e,
    };
    for step in 0..=n_steps {
        if step > 0 {
            hybrid_step_with_forces(
                &mut bodies,
                &mut particles,
                &mut active,
                DT_DAYS,
                &config,
                &extra,
            );
        }
        if step % snap_every != 0 {
            continue;
        }
        let t = step as f64 * DT_DAYS / YEAR_DAYS / 1e6;
        if !active[0] {
            track.lost_myr.get_or_insert(t);
            break;
        }
        let el = cartesian_to_elements(&particles[0], GM_SUN);
        if el.e >= 1.0 || el.a <= 0.0 {
            track.lost_myr.get_or_insert(t);
            break;
        }
        track.t_myr.push(t);
        track.a_au.push(el.a);
        track.q_au.push(el.a * (1.0 - el.e));
        track.i_deg.push(el.i * RAD2DEG);
        track.omega_rad.push(el.omega);
    }
    track
}

pub fn export() -> Value {
    let tracks: Vec<Track> = thread::scope(|s| {
        let handles: Vec<_> = SIX
            .iter()
            .map(|&name| s.spawn(move || integrate(name)))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    // Grouping of the arguments of perihelion of whichever objects are still
    // in the run, snapshot by snapshot.
    let n_snap = tracks.iter().map(|t| t.t_myr.len()).max().unwrap_or(0);
    let mut group_t = Vec::new();
    let mut group_r = Vec::new();
    let mut group_n = Vec::new();
    for k in 0..n_snap {
        let angles: Vec<f64> = tracks
            .iter()
            .filter(|t| k < t.omega_rad.len())
            .map(|t| t.omega_rad[k])
            .collect();
        if angles.len() < 2 {
            break;
        }
        group_t.push(k as f64 * SNAPSHOT_MYR);
        group_r.push(grouping_of(&angles).r_bar);
        group_n.push(angles.len());
    }

    let n_unstable = tracks.iter().filter(|t| t.unstable_myr().is_some()).count();
    let objects: Vec<Value> = tracks
        .iter()
        .map(|t| {
            json!({
                "name": t.name,
                "t_myr": t.t_myr,
                "a_au": t.a_au,
                "q_au": t.q_au,
                "i_deg": t.i_deg,
                "omega_deg": t.omega_rad.iter().map(|w| w * RAD2DEG).collect::<Vec<_>>(),
                "lost_myr": t.lost_myr,
                "unstable_myr": t.unstable_myr(),
                "a_excursion": t.excursion(),
                "varpi0_deg": t.varpi0_deg,
                "e0": t.e0,
            })
        })
        .collect();

    let p9 = nominal_planet();
    json!({
        "t_myr": T_MYR,
        "planet": {
            "mass_earth": p9.mass_earth,
            "a_au": p9.a,
            "e": p9.e,
            "i_deg": p9.i * RAD2DEG,
            "varpi_deg": ((p9.omega + p9.omega_big) * RAD2DEG).rem_euclid(360.0),
        },
        "objects": objects,
        "n_objects": SIX.len(),
        "n_unstable": n_unstable,
        "n_stable": SIX.len() - n_unstable,
        "unstable_excursion": UNSTABLE_EXCURSION,
        "omega_grouping": {"t_myr": group_t, "r_bar": group_r, "n": group_n},
        "omega_r_bar_start": group_r.first(),
        "omega_r_bar_end": group_r.last(),
    })
}

//! Film export for `p9-2016-inclination-instability`: the numbers its scene and ledger entry draw.

use p9_2016_inclination_instability::growth::{
    C_GROWTH, efolding_time_myr, orbital_period_typed, secular_time_typed,
};
use p9_2016_inclination_instability::suppression::{
    Stability, classify, critical_mass_earth, giant_planet_j2_eff,
};
use p9_core::constants::{EARTH_MASS_SOLAR, TWO_PI, YEAR_DAYS};
use p9_core::units::days;
use serde_json::{Value, json};

/// Disc used for the headline numbers: the upper end of the paper's 1-10 Earth
/// mass range, at a characteristic 250 AU with the crate's fiducial e = 0.6.
const FIDUCIAL_MASS_EARTH: f64 = 10.0;
const FIDUCIAL_A_AU: f64 = 250.0;
const FIDUCIAL_E: f64 = 0.6;

/// Semi-major axes for the e-folding curves.
const CURVE_A_AU: [f64; 3] = [100.0, 250.0, 500.0];
/// Disc eccentricities for the critical-mass curves.
const CURVE_E: [f64; 3] = [0.2, 0.6, 0.8];

fn log_grid(lo: f64, hi: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|k| lo * (hi / lo).powf(k as f64 / (n - 1) as f64))
        .collect()
}

pub fn export() -> Value {
    let masses = log_grid(0.1, 100.0, 61);
    let efold: Vec<Value> = CURVE_A_AU
        .iter()
        .map(|&a| {
            let tau: Vec<f64> = masses
                .iter()
                .map(|&m| efolding_time_myr(m * EARTH_MASS_SOLAR, a))
                .collect();
            json!({"a_au": a, "tau_myr": tau})
        })
        .collect();

    let radii = log_grid(50.0, 1000.0, 61);
    let critical: Vec<Value> = CURVE_E
        .iter()
        .map(|&e| {
            let m: Vec<f64> = radii.iter().map(|&a| critical_mass_earth(a, e)).collect();
            json!({"e": e, "m_crit_earth": m})
        })
        .collect();

    let m_fid = FIDUCIAL_MASS_EARTH * EARTH_MASS_SOLAR;
    let period_yr = (orbital_period_typed(FIDUCIAL_A_AU) / days(1.0)).value / YEAR_DAYS;
    let t_sec_myr = (secular_time_typed(m_fid, FIDUCIAL_A_AU) / days(1.0)).value / YEAR_DAYS / 1e6;
    let tau_myr = efolding_time_myr(m_fid, FIDUCIAL_A_AU);

    // Radius beyond which each of the paper's bracketing disc masses wins
    // against the giant planets' differential precession.
    let unstable_beyond = |m_earth: f64| -> Option<f64> {
        radii
            .iter()
            .copied()
            .find(|&a| classify(m_earth * EARTH_MASS_SOLAR, a, FIDUCIAL_E) == Stability::Unstable)
    };

    json!({
        "fiducial": {
            "mass_earth": FIDUCIAL_MASS_EARTH,
            "a_au": FIDUCIAL_A_AU,
            "e": FIDUCIAL_E,
            "period_yr": period_yr,
            "t_sec_myr": t_sec_myr,
            "tau_myr": tau_myr,
        },
        "tau_over_t_sec": C_GROWTH / TWO_PI,
        "tau_fiducial_myr": tau_myr,
        "tau_low_mass_myr": efolding_time_myr(EARTH_MASS_SOLAR, FIDUCIAL_A_AU),
        "m_crit_fiducial_earth": critical_mass_earth(FIDUCIAL_A_AU, FIDUCIAL_E),
        "m_crit_100au_earth": critical_mass_earth(100.0, FIDUCIAL_E),
        "j2_eff_au2": giant_planet_j2_eff(),
        "unstable_beyond_au": {
            "m1": unstable_beyond(1.0),
            "m10": unstable_beyond(10.0),
        },
        "paper_mass_range_earth": [1.0, 10.0],
        "masses_earth": masses,
        "efolding": efold,
        "radii_au": radii,
        "critical_mass": critical,
    })
}

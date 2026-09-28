//! Film export for `p9-2016-obliquity-lai`: the numbers its scene and ledger entry draw.

use p9_2016_obliquity::secular_hamiltonian::{SecularParams, SpinOrbitState, integrate_obliquity};
use p9_2016_obliquity_lai::planets::{l_jupiter, total_orbital_angular_momentum};
use p9_2016_obliquity_lai::precession::published::{
    A_TILDE_RANGE_AU, OBSERVED_OBLIQUITY_DEG, OMEGA_L_PERIOD_GYR, OMEGA_SPIN_PERIOD_GYR,
};
use p9_2016_obliquity_lai::precession::{
    PlanetNine, corotating_frequencies_typed, omega_l_typed, omega_spin_precession_typed,
    solar_obliquity_deg, solar_obliquity_typed,
};
use p9_2016_obliquity_lai::survey::{
    NodeGeometry, required_a_tilde_typed, solve_a_tilde_for_obliquity,
};
use p9_core::constants::{EARTH_MASS_SOLAR, GYR_DAYS, YEAR_DAYS};
use p9_core::initial_conditions::giant_planets::{
    giant_planet_angular_momentum, p9_angular_momentum,
};
use p9_core::units::{au, days, radians};
use serde_json::{Value, json};
use std::f64::consts::PI;

/// Age of the Solar System (Gyr).
const AGE_GYR: f64 = 4.5;

/// Solar rotation parameter P/lambda of Lai (2016) Fig. 1 (days).
const SPIN_PERIOD_DAYS: f64 = 20.0;

/// Reference perturber of Lai (2016) Eqs. 5 and 8: 10 Earth masses at an
/// effective semi-major axis of 400 AU, solar rotation period 10 days.
const REF_MASS_EARTH: f64 = 10.0;
const REF_A_TILDE_AU: f64 = 400.0;
const REF_SPIN_PERIOD_DAYS: f64 = 10.0;

/// Perturber compared against the numerical integration of Bailey et al.
const COMPARE: (f64, f64, f64, f64) = (15.0, 500.0, 0.5, 20.0);

/// Node offsets between Planet Nine and the solar equator spanned by Lai's
/// Fig. 1 (degrees).
const DELTA_NODE_DEG: [f64; 3] = [12.0, 45.0, 52.0];

/// Inclinations whose tilt-versus-distance curves are drawn (degrees).
const INCLINATIONS_DEG: [f64; 3] = [20.0, 30.0, 40.0];

/// Planet Nine masses of Lai's Fig. 1 (Earth masses).
const CURVE_MASSES_EARTH: [f64; 2] = [10.0, 20.0];

fn period_gyr(rate_rad_per_day: f64) -> f64 {
    2.0 * PI / rate_rad_per_day.abs() / GYR_DAYS
}

fn circular(mass_earth: f64, a_tilde_au: f64, inclination_deg: f64) -> PlanetNine {
    PlanetNine {
        mass_earth,
        a_au: a_tilde_au,
        e: 0.0,
        inclination_rad: inclination_deg.to_radians(),
    }
}

pub fn export() -> Value {
    let rate = radians(1.0) / days(1.0);

    // The two frequencies of the closed form, at Lai's reference scalings.
    let reference = circular(REF_MASS_EARTH, REF_A_TILDE_AU, 20.0);
    let omega_l = (omega_l_typed(&reference) / rate).value;
    let omega_spin = (omega_spin_precession_typed(REF_SPIN_PERIOD_DAYS) / rate).value;

    // Tilt after 4.5 Gyr against effective semi-major axis.
    let a_tilde: Vec<f64> = (0..=200).map(|k| 250.0 + 2.5 * k as f64).collect();
    let tilt_curves: Vec<Value> = CURVE_MASSES_EARTH
        .iter()
        .flat_map(|&mass| INCLINATIONS_DEG.iter().map(move |&inc| (mass, inc)))
        .map(|(mass, inc)| {
            let tilt: Vec<f64> = a_tilde
                .iter()
                .map(|&a| solar_obliquity_deg(&circular(mass, a, inc), SPIN_PERIOD_DAYS))
                .collect();
            json!({
                "mass_earth": mass,
                "inclination_deg": inc,
                "tilt_deg": tilt,
                "a_tilde_for_observed_au": solve_a_tilde_for_obliquity(
                    mass,
                    inc.to_radians(),
                    0.0,
                    SPIN_PERIOD_DAYS,
                    OBSERVED_OBLIQUITY_DEG,
                ),
            })
        })
        .collect();

    // Lai's Eq. 17: the effective semi-major axis that yields the observed
    // tilt, against inclination, for each node offset.
    let theta_deg: Vec<f64> = (0..=70).map(|k| 5.0 + 0.5 * k as f64).collect();
    let contours: Vec<Value> = DELTA_NODE_DEG
        .iter()
        .map(|&dn| {
            let geom = NodeGeometry::new(OBSERVED_OBLIQUITY_DEG, dn.to_radians());
            let a: Vec<f64> = theta_deg
                .iter()
                .map(|t| {
                    (required_a_tilde_typed(REF_MASS_EARTH, t.to_radians(), &geom) / au(1.0)).value
                })
                .collect();
            json!({"delta_node_deg": dn, "a_tilde_au": a})
        })
        .collect();
    let geom45 = NodeGeometry::new(OBSERVED_OBLIQUITY_DEG, 45f64.to_radians());
    let a_tilde_eq17 =
        (required_a_tilde_typed(REF_MASS_EARTH, 20f64.to_radians(), &geom45) / au(1.0)).value;

    // Closed form against the numerical integration of Bailey et al. for the
    // same perturber.
    let (mass, a9, e9, inc) = COMPARE;
    let p9 = PlanetNine {
        mass_earth: mass,
        a_au: a9,
        e: e9,
        inclination_rad: inc.to_radians(),
    };
    let (omega_y, omega_z) = corotating_frequencies_typed(&p9, SPIN_PERIOD_DAYS);
    let t_gyr: Vec<f64> = (0..=90).map(|k| AGE_GYR * k as f64 / 90.0).collect();
    let analytic: Vec<f64> = t_gyr
        .iter()
        .map(|&t| (solar_obliquity_typed(&p9, SPIN_PERIOD_DAYS, t) / radians(1.0)).value)
        .map(f64::to_degrees)
        .collect();
    let m9_solar = mass * EARTH_MASS_SOLAR;
    let snaps = integrate_obliquity(
        SpinOrbitState::from_inclinations(
            inc.to_radians(),
            PI,
            giant_planet_angular_momentum(),
            p9_angular_momentum(m9_solar, a9, e9),
        ),
        &SecularParams {
            m9_solar,
            a9,
            e9,
            t_total: AGE_GYR * GYR_DAYS,
            dt: 5e4 * YEAR_DAYS,
        },
        0.25 * GYR_DAYS,
    );

    json!({
        "observed_obliquity_deg": OBSERVED_OBLIQUITY_DEG,
        "age_gyr": AGE_GYR,
        "spin_period_days": SPIN_PERIOD_DAYS,
        "reference": {
            "mass_earth": REF_MASS_EARTH,
            "a_tilde_au": REF_A_TILDE_AU,
            "spin_period_days": REF_SPIN_PERIOD_DAYS,
        },
        "plane_precession_period_gyr": period_gyr(omega_l),
        "spin_precession_period_gyr": period_gyr(omega_spin),
        "l_total_over_l_jupiter": total_orbital_angular_momentum() / l_jupiter(),
        "a_tilde_au": a_tilde,
        "tilt_curves": tilt_curves,
        "theta_deg": theta_deg,
        "contours": contours,
        "a_tilde_eq17_au": a_tilde_eq17,
        "compare": {
            "mass_earth": mass,
            "a_au": a9,
            "e": e9,
            "inclination_deg": inc,
            "a_tilde_au": p9.a_tilde_au(),
            "drive_period_gyr": period_gyr((omega_y / rate).value),
            "mismatch_period_gyr": period_gyr((omega_z / rate).value),
            "t_gyr": t_gyr,
            "analytic_deg": analytic,
            "analytic_final_deg": solar_obliquity_deg(&p9, SPIN_PERIOD_DAYS),
            "numerical_t_gyr": snaps.iter().map(|s| s.t / GYR_DAYS).collect::<Vec<_>>(),
            "numerical_deg": snaps.iter().map(|s| s.obliquity.to_degrees()).collect::<Vec<_>>(),
            "numerical_final_deg": snaps.last().map(|s| s.obliquity.to_degrees()),
        },
        "published": {
            "plane_precession_period_gyr": OMEGA_L_PERIOD_GYR,
            "spin_precession_period_gyr": OMEGA_SPIN_PERIOD_GYR,
            "a_tilde_range_au": A_TILDE_RANGE_AU,
        },
    })
}

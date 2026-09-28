//! Film export for `p9-2026-stellar-companions`: the numbers its scene and ledger entry draw.

use p9_2024_panstarrs::combined_exclusion::UpdatedParameters;
use p9_2026_stellar_companions::envelope::mass_envelope_earth;
use p9_2026_stellar_companions::halo::{dm_halo_mass_within_kg, dm_halo_mass_within_solar};
use p9_2026_stellar_companions::{
    ANCHOR_DISTANCE_AU, ANCHOR_MASS_EARTH, M_PLUTO_KG, RHO_LOCAL_DM_MSUN_PC3,
};
use p9_core::constants::{MASS_JUPITER_SOLAR, MASS_SATURN_SOLAR};
use p9_core::units::{EARTH_MASS_KG, SOLAR_MASS_KG};
use serde_json::{Value, json};

/// Benakli (2026), Eq. 32: the envelope in Earth masses at 1000 AU, and the
/// rows of the paper's table of envelope crossings (distance AU, mass M_earth,
/// label).
const PUBLISHED_ENVELOPE_1000: f64 = 38.0;
const PUBLISHED_TABLE: [(f64, f64, &str); 8] = [
    (300.0, 1.0, "Earth mass"),
    (500.0, 4.8, "super-Earth"),
    (650.0, 10.5, "INPOP19a anchor"),
    (766.0, 17.0, "Neptune mass"),
    (1000.0, 38.0, "sub-Saturn"),
    (1355.0, 95.0, "Saturn mass"),
    (1500.0, 129.0, "sub-Jovian"),
    (2026.0, 318.0, "Jupiter mass"),
];

/// Distance (AU) at which the envelope admits `mass_earth`.
fn distance_for_mass(mass_earth: f64) -> f64 {
    let (mut lo, mut hi) = (10.0_f64, 1.0e5_f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if mass_envelope_earth(mid) < mass_earth {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn export() -> Value {
    let earth_per_solar = SOLAR_MASS_KG / EARTH_MASS_KG;
    let jupiter = MASS_JUPITER_SOLAR * earth_per_solar;
    let saturn = MASS_SATURN_SOLAR * earth_per_solar;

    let log_d: Vec<f64> = (0..=60).map(|k| 2.0 + 1.6 * k as f64 / 60.0).collect();
    let distances: Vec<f64> = log_d.iter().map(|l| 10f64.powf(*l)).collect();
    let envelope: Vec<f64> = distances.iter().map(|&d| mass_envelope_earth(d)).collect();
    let halo: Vec<f64> = distances
        .iter()
        .map(|&d| dm_halo_mass_within_kg(d) / EARTH_MASS_KG)
        .collect();

    let p9 = UpdatedParameters::paper_values();
    let table: Vec<Value> = PUBLISHED_TABLE
        .iter()
        .map(|(d, m, label)| {
            json!({
                "distance_au": d,
                "published_mass_earth": m,
                "mass_earth": mass_envelope_earth(*d),
                "label": label,
            })
        })
        .collect();

    json!({
        "anchor": {"mass_earth": ANCHOR_MASS_EARTH, "distance_au": ANCHOR_DISTANCE_AU},
        "distance_au": distances,
        "envelope_earth": envelope,
        "halo_earth": halo,
        "envelope_300": mass_envelope_earth(300.0),
        "envelope_1000": mass_envelope_earth(1000.0),
        "envelope_2000": mass_envelope_earth(2000.0),
        "distance_for_saturn": distance_for_mass(saturn),
        "distance_for_jupiter": distance_for_mass(jupiter),
        "jupiter_earth": jupiter,
        "saturn_earth": saturn,
        "pluto_earth": M_PLUTO_KG / EARTH_MASS_KG,
        "halo_1000_kg": dm_halo_mass_within_kg(1000.0),
        "halo_1000_pluto": dm_halo_mass_within_kg(1000.0) / M_PLUTO_KG,
        "halo_1000_earth": dm_halo_mass_within_solar(1000.0) * earth_per_solar,
        "halo_shortfall_1000": mass_envelope_earth(1000.0)
            / (dm_halo_mass_within_solar(1000.0) * earth_per_solar),
        "rho_dm_msun_pc3": RHO_LOCAL_DM_MSUN_PC3,
        "planet_nine": {"mass_earth": p9.mass_earth_median, "a_au": p9.a_median},
        "table": table,
        "published": {"envelope_1000": PUBLISHED_ENVELOPE_1000},
    })
}

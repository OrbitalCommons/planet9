//! Film export for `p9-2022-iras-candidate`: the numbers its scene and ledger entry draw.

use p9_2022_iras_candidate::chance::{
    N_REJECT_SINGLE_HCON, N_UNIDENTIFIED_60UM, expected_chance_associations,
    surface_density_per_arcmin2,
};
use p9_2022_iras_candidate::distance::{expected_flux_jy, implied_distance_au};
use p9_2022_iras_candidate::{LAMBDA_60UM_M, REF_CANDIDATE};
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::{
    ecliptic_to_equatorial_deg, ecliptic_vec_to_equatorial_deg, equatorial_to_galactic,
};
use p9_core::data::reference_population::generate_reference_population;
use p9_core::types::OrbitalElements;
use rand::SeedableRng;
use serde_json::{Value, json};

/// The candidate as published by Rowan-Robinson (arXiv:2111.03831, Sec. 6.5):
/// centre of the fitted parallactic ellipse (ecliptic J2000, degrees), 60 µm
/// flux, motion between the first and last IRAS passes, the quoted galactic
/// latitude, and the follow-up annulus.
const CANDIDATE_ECLIPTIC_LON_DEG: f64 = 13.62;
const CANDIDATE_ECLIPTIC_LAT_DEG: f64 = 71.62;
const CANDIDATE_FLUX_60UM_JY: f64 = 0.57;
const CANDIDATE_MOTION_ARCMIN: f64 = 20.3;
const CANDIDATE_MOTION_WEEKS: f64 = 11.9;
const PUBLISHED_GALACTIC_LAT_DEG: f64 = 9.0;
const FOLLOW_UP_ANNULUS_DEG: (f64, f64) = (2.5, 4.0);

/// 90%-completeness limit of the single-HCON search at 60 µm (Jy).
const COMPLETENESS_60UM_JY: f64 = 0.45;

/// Associations the paper examined by eye: 532 pairs at 2-25 arcmin,
/// 75 triplets and 76 close pairs.
const PUBLISHED_PAIRS: f64 = 532.0;
const PUBLISHED_TRIPLETS: f64 = 75.0;
const PUBLISHED_CLOSE_PAIRS: f64 = 76.0;
const PAIR_ANNULUS_ARCMIN: (f64, f64) = (2.0, 25.0);

/// Fraction of the sky IRAS surveyed with hours-confirmed passes.
const IRAS_SKY_COVERAGE: f64 = 0.96;

/// Predicted Planet Nines drawn for the sky comparison.
const N_POPULATION: usize = 500;

/// Temperature at which a body of `mass_earth` at `distance_au` emits
/// `flux_jy` at 60 µm, by bisection on the crate's forward model.
fn temperature_for_flux(mass_earth: f64, distance_au: f64, flux_jy: f64) -> f64 {
    let (mut lo, mut hi) = (15.0_f64, 150.0_f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if expected_flux_jy(mass_earth, distance_au, mid, LAMBDA_60UM_M) < flux_jy {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn export() -> Value {
    let c = REF_CANDIDATE;
    let masses = [c.mass_lo_earth, c.mass_earth, c.mass_hi_earth];

    let distances: Vec<f64> = (0..=70).map(|k| 100.0 + 5.0 * k as f64).collect();
    let curves: Vec<Value> = masses
        .iter()
        .map(|&m| {
            let flux: Vec<f64> = distances
                .iter()
                .map(|&d| expected_flux_jy(m, d, c.t_eff_k, LAMBDA_60UM_M))
                .collect();
            json!({
                "mass_earth": m,
                "flux_60um_jy": flux,
                "implied_distance_au":
                    implied_distance_au(CANDIDATE_FLUX_60UM_JY, c.t_eff_k, m, LAMBDA_60UM_M),
                "temperature_at_published_distance_k":
                    temperature_for_flux(m, c.distance_au, CANDIDATE_FLUX_60UM_JY),
            })
        })
        .collect();

    let (ra, dec) =
        ecliptic_to_equatorial_deg(CANDIDATE_ECLIPTIC_LON_DEG, CANDIDATE_ECLIPTIC_LAT_DEG);
    let (_, gal_b) = equatorial_to_galactic(ra.to_radians(), dec.to_radians());

    let mut rng = rand::rngs::StdRng::seed_from_u64(1983);
    let mut max_abs_lat = 0.0_f64;
    let population: Vec<Value> = generate_reference_population(N_POPULATION, &mut rng)
        .iter()
        .map(|obj| {
            let pos = OrbitalElements {
                a: obj.a,
                e: obj.e,
                i: obj.i,
                omega: obj.omega,
                omega_big: obj.omega_big,
                mean_anomaly: obj.mean_anomaly,
            }
            .to_state_vector(GM_SUN)
            .pos;
            let (ra, dec) = ecliptic_vec_to_equatorial_deg(&pos);
            let lat = (pos.z / pos.norm()).asin().to_degrees();
            max_abs_lat = max_abs_lat.max(lat.abs());
            json!({
                "ra_deg": ra,
                "dec_deg": dec,
                "dist_au": pos.norm(),
                "mass_earth": obj.mass,
                "ecliptic_lat_deg": lat,
            })
        })
        .collect();

    let sigma = surface_density_per_arcmin2(N_REJECT_SINGLE_HCON, IRAS_SKY_COVERAGE);
    let chance = expected_chance_associations(
        N_UNIDENTIFIED_60UM,
        sigma,
        PAIR_ANNULUS_ARCMIN.0,
        PAIR_ANNULUS_ARCMIN.1,
    );

    json!({
        "candidate": {
            "ra_deg": ra,
            "dec_deg": dec,
            "ecliptic_lon_deg": CANDIDATE_ECLIPTIC_LON_DEG,
            "ecliptic_lat_deg": CANDIDATE_ECLIPTIC_LAT_DEG,
            "galactic_lat_deg": gal_b.to_degrees(),
            "published_galactic_lat_deg": PUBLISHED_GALACTIC_LAT_DEG,
            "flux_60um_jy": CANDIDATE_FLUX_60UM_JY,
            "motion_arcmin": CANDIDATE_MOTION_ARCMIN,
            "motion_weeks": CANDIDATE_MOTION_WEEKS,
            "follow_up_annulus_deg": [FOLLOW_UP_ANNULUS_DEG.0, FOLLOW_UP_ANNULUS_DEG.1],
        },
        "published_distance_au": c.distance_au,
        "published_distance_err_au": c.distance_err_au,
        "published_mass_lo": c.mass_lo_earth,
        "published_mass_hi": c.mass_hi_earth,
        "model_temp_k": c.t_eff_k,
        "crate_reference_flux_jy": c.flux_60um_jy,
        "completeness_60um_jy": COMPLETENESS_60UM_JY,
        "implied_distance_au":
            implied_distance_au(CANDIDATE_FLUX_60UM_JY, c.t_eff_k, c.mass_earth, LAMBDA_60UM_M),
        "parallax_semi_major_arcmin": (1.0 / c.distance_au).atan().to_degrees() * 60.0,
        "distance_au": distances,
        "curves": curves,
        "population": population,
        "population_max_ecliptic_lat_deg": max_abs_lat,
        "chance_associations": chance,
        "published_associations": PUBLISHED_PAIRS + PUBLISHED_TRIPLETS + PUBLISHED_CLOSE_PAIRS,
        "published_pairs": PUBLISHED_PAIRS,
        "published_triplets": PUBLISHED_TRIPLETS,
        "published_close_pairs": PUBLISHED_CLOSE_PAIRS,
        "survivors": 1,
    })
}

//! Combined Planet Nine "still-viable" parameter-space map.
//!
//! For a grid of masses we compute, for each survey, the maximum heliocentric
//! distance at which a Planet Nine of that mass would have been detected (its
//! *reach*). A (mass, distance) point is RULED OUT if it sits below the reach
//! of any survey that has actually taken data; the surviving region is the set
//! of (mass, distance) pairs that are too faint (too small and/or too far) to
//! have been caught yet.
//!
//! All reaches are computed from the workspace's own physical models — reflected
//! optical light via [`p9_core::analysis::photometry`], thermal/mm flux via
//! [`p9_core::analysis::thermal`], and the WISE W1 model via the dedicated
//! `p9-2018-wise-search` crate — so this map tracks the reproductions, not a
//! hand-drawn cartoon. Output is JSON consumed by `scripts/plot_viability.py`.

use p9_core::analysis::photometry::mass_radius_neptunian;
use p9_core::analysis::photometry::{ALBEDO_NEPTUNE, planet_apparent_magnitude};
use p9_core::analysis::thermal::{
    C_LIGHT, effective_temp, flux_to_magnitude, max_detectable_distance, thermal_flux_jy,
};
use p9_core::constants::EARTH_RADIUS_KM;
use p9_core::types::{P9Params, helio_distance_at_true_anomaly};
use p9_core::units::au;

use p9_2018_wise_search::detectability::max_detectable_distance as wise_reach;
use p9_2018_wise_search::survey_model::WiseSurvey;

use p9_2025_simons_forecast::bands::{ACT_SENSITIVITY, SO_SENSITIVITY};
use p9_2025_simons_forecast::forecast::max_detectable_distance as mm_reach;

use serde::Serialize;

/// Distance bracket for all reach bisections (AU).
const D_MIN: f64 = 60.0;
const D_MAX: f64 = 20_000.0;
/// Internal-heat temperature floor for a cold giant (K); matches the value used
/// across the thermal reproductions.
const INTERNAL_TEMP_K: f64 = 40.0;

/// Body radius in metres from the Neptune-anchored mass–radius relation.
fn radius_m(mass_earth: f64) -> f64 {
    mass_radius_neptunian(mass_earth) * EARTH_RADIUS_KM * 1_000.0
}

/// Reflected-light reach (AU) of an optical survey at limiting magnitude `depth`.
fn optical_reach(mass_earth: f64, depth: f64) -> Option<f64> {
    max_detectable_distance(D_MIN, D_MAX, depth, |d| {
        planet_apparent_magnitude(mass_earth, ALBEDO_NEPTUNE, d)
    })
}

/// Thermal reach (AU) for a flux-limited survey at frequency `nu_hz` with a
/// point-source limit `flux_limit_jy`. Uses the same monotone-bisection helper
/// by mapping flux through a (zero-point-cancelling) magnitude.
fn thermal_reach(mass_earth: f64, nu_hz: f64, flux_limit_jy: f64) -> Option<f64> {
    let r_m = radius_m(mass_earth);
    let depth = flux_to_magnitude(flux_limit_jy, 1.0);
    max_detectable_distance(D_MIN, D_MAX, depth, |d| {
        let t = effective_temp(d, ALBEDO_NEPTUNE, INTERNAL_TEMP_K);
        flux_to_magnitude(thermal_flux_jy(t, r_m, d, nu_hz), 1.0)
    })
}

#[derive(Serialize)]
struct SurveyReach {
    name: &'static str,
    kind: &'static str, // "optical" | "thermal" | "mm"
    /// Whether this is real archival data ("data") or a forecast ("forecast").
    status: &'static str,
    /// "allsky" (≳3/4 sky — excludes regardless of where P9 sits) or "partial"
    /// (deep but small footprint — only excludes in-footprint objects).
    coverage: &'static str,
    /// Approximate fraction of sky searched (for the footprint caveat).
    sky_fraction: f64,
    reference: &'static str,
    /// reach[i] = max detectable heliocentric distance (AU) at mass[i]; null if
    /// the body is undetectable anywhere in the bracket at that mass.
    reach: Vec<Option<f64>>,
}

#[derive(Serialize)]
struct NominalOrbit {
    name: &'static str,
    mass_earth: f64,
    a_au: f64,
    e: f64,
    perihelion_au: f64,
    aphelion_au: f64,
    /// Heliocentric distance at the Cassini-favored current true anomaly (AU).
    current_distance_au: f64,
}

#[derive(Serialize)]
struct MagCurve {
    mass_earth: f64,
    /// Apparent reflected V magnitude at each distance in `distance_au`.
    v_mag: Vec<f64>,
}

#[derive(Serialize)]
struct Dataset {
    generated_by: &'static str,
    mass_earth: Vec<f64>,
    surveys: Vec<SurveyReach>,
    /// Upper envelope over ALL real-data surveys (footprint ignored): the
    /// deepest anyone has seen at each mass. Below this, P9 is excluded *if* it
    /// fell in the relevant footprint.
    detection_envelope_au: Vec<Option<f64>>,
    /// Upper envelope over near-all-sky real-data surveys only: the robust
    /// exclusion floor — below this, P9 is ruled out *regardless* of sky
    /// position. This is the boundary of the genuinely-viable region.
    allsky_envelope_au: Vec<Option<f64>>,
    nominal_orbits: Vec<NominalOrbit>,
    /// Favoured current-distance band (AU) from the Cassini true-anomaly window.
    favored_distance_band_au: [f64; 2],
    // Companion magnitude-vs-distance panel.
    distance_au: Vec<f64>,
    mag_curves: Vec<MagCurve>,
    survey_optical_depths: Vec<(&'static str, f64)>,
}

fn nominal(name: &'static str, p: P9Params) -> NominalOrbit {
    // Cassini-favoured current true anomaly window ~108-129 deg; use the
    // midpoint to place the body's present heliocentric distance.
    let nu_mid = 118.5_f64.to_radians();
    NominalOrbit {
        name,
        mass_earth: p.mass_earth,
        a_au: p.a,
        e: p.e,
        perihelion_au: p.a * (1.0 - p.e),
        aphelion_au: p.a * (1.0 + p.e),
        current_distance_au: helio_distance_at_true_anomaly(&p, nu_mid),
    }
}

fn main() {
    // Mass grid: 1-20 Earth masses, fine enough for smooth curves.
    let mass_earth: Vec<f64> = (0..=190).map(|k| 1.0 + k as f64 * 0.1).collect();

    // Optical reflected-light surveys (limiting r/V magnitude, real archival
    // data). Fields: name, status, coverage, sky_fraction, reference, depth.
    let optical: &[(&str, &str, &str, f64, &str, f64)] = &[
        ("CRTS", "data", "allsky", 0.75, "Catalina RTS, V~19.5", 19.5),
        (
            "ZTF",
            "data",
            "allsky",
            0.70,
            "Brown & Batygin 2022, r~20.5",
            20.5,
        ),
        (
            "Pan-STARRS1 3pi",
            "data",
            "allsky",
            0.75,
            "Brown/Holman/Batygin 2024, r~21.5",
            21.5,
        ),
        (
            "TESS (stack)",
            "data",
            "partial",
            0.40,
            "Payne/Holman/Pal 2019, I_C~22.0",
            22.0,
        ),
        (
            "DES",
            "data",
            "partial",
            0.12,
            "Belyakov et al. 2022/2024, r~23.8",
            23.8,
        ),
        (
            "Rubin/LSST (visit)",
            "forecast",
            "partial",
            0.45,
            "Schwamb et al. 2023, r~24.5",
            24.5,
        ),
        (
            "Rubin/LSST (10yr)",
            "forecast",
            "partial",
            0.45,
            "LSST coadd, r~27.0",
            27.0,
        ),
    ];

    let mut surveys: Vec<SurveyReach> = Vec::new();

    for &(name, status, coverage, sky_fraction, reference, depth) in optical {
        surveys.push(SurveyReach {
            name,
            kind: "optical",
            status,
            coverage,
            sky_fraction,
            reference,
            reach: mass_earth
                .iter()
                .map(|&m| optical_reach(m, depth))
                .collect(),
        });
    }

    // WISE W1 thermal-IR, via the dedicated reproduction (cold-body W1 model).
    let wise = WiseSurvey::default();
    surveys.push(SurveyReach {
        name: "WISE W1",
        kind: "thermal",
        status: "data",
        coverage: "allsky",
        sky_fraction: 0.95,
        reference: "Meisner et al. 2018, W1~16.5",
        reach: mass_earth.iter().map(|&m| wise_reach(&wise, m)).collect(),
    });

    // Far-IR archival surveys (thermal blackbody near the ~40 K peak).
    let iras_nu = C_LIGHT / 60e-6;
    let akari_nu = C_LIGHT / 90e-6;
    surveys.push(SurveyReach {
        name: "IRAS 60um",
        kind: "thermal",
        status: "data",
        coverage: "allsky",
        sky_fraction: 0.96,
        reference: "Phan et al. 2025, FSC ~0.2 Jy",
        reach: mass_earth
            .iter()
            .map(|&m| thermal_reach(m, iras_nu, 0.2))
            .collect(),
    });
    surveys.push(SurveyReach {
        name: "AKARI 90um",
        kind: "thermal",
        status: "data",
        coverage: "allsky",
        sky_fraction: 0.94,
        reference: "Phan et al. 2025, MU ~0.55 Jy",
        reach: mass_earth
            .iter()
            .map(|&m| thermal_reach(m, akari_nu, 0.55))
            .collect(),
    });

    // Millimetre surveys via the Simons forecast reproduction.
    let mm: &[(
        &str,
        &str,
        &str,
        f64,
        &str,
        &p9_2025_simons_forecast::bands::MmSensitivity,
    )] = &[
        (
            "ACT (mm)",
            "data",
            "partial",
            0.40,
            "ACT 229 GHz, 8 mJy",
            &ACT_SENSITIVITY,
        ),
        (
            "Simons Obs. (mm)",
            "forecast",
            "partial",
            0.40,
            "SO 280 GHz, 4.4 mJy",
            &SO_SENSITIVITY,
        ),
    ];
    for &(name, status, coverage, sky_fraction, reference, sens) in mm {
        surveys.push(SurveyReach {
            name,
            kind: "mm",
            status,
            coverage,
            sky_fraction,
            reference,
            reach: mass_earth
                .iter()
                .map(|&m| Some((mm_reach(m, sens) / au(1.0)).value))
                .collect(),
        });
    }

    // Envelope helper: max reach at each mass over surveys passing `keep`.
    let envelope = |keep: &dyn Fn(&SurveyReach) -> bool| -> Vec<Option<f64>> {
        (0..mass_earth.len())
            .map(|i| {
                surveys
                    .iter()
                    .filter(|s| keep(s))
                    .filter_map(|s| s.reach[i])
                    .fold(None::<f64>, |acc, d| Some(acc.map_or(d, |a| a.max(d))))
            })
            .collect()
    };
    // All real data (footprint ignored) vs near-all-sky real data (robust floor).
    let detection_envelope_au = envelope(&|s| s.status == "data");
    let allsky_envelope_au = envelope(&|s| s.status == "data" && s.coverage == "allsky");

    let nominal_orbits = vec![
        nominal("Batygin & Brown 2016", P9Params::nominal_2016()),
        nominal("Batygin et al. 2019 (review)", P9Params::revised_2019()),
        nominal("Brown & Batygin 2021 (MCMC)", P9Params::mcmc_2021()),
    ];

    // Companion: apparent V magnitude vs heliocentric distance for fiducial masses.
    let distance_au: Vec<f64> = (0..=120).map(|k| 100.0 + k as f64 * 25.0).collect();
    let mag_curves: Vec<MagCurve> = [3.0, 5.0, 7.0, 10.0, 15.0]
        .iter()
        .map(|&m| MagCurve {
            mass_earth: m,
            v_mag: distance_au
                .iter()
                .map(|&d| planet_apparent_magnitude(m, ALBEDO_NEPTUNE, d))
                .collect(),
        })
        .collect();

    let dataset = Dataset {
        generated_by: "p9-viability",
        mass_earth,
        surveys,
        detection_envelope_au,
        allsky_envelope_au,
        nominal_orbits,
        favored_distance_band_au: [400.0, 800.0],
        distance_au,
        mag_curves,
        survey_optical_depths: vec![
            ("ZTF", 20.5),
            ("PS1", 21.5),
            ("TESS", 22.0),
            ("DES", 23.8),
            ("LSST 10yr", 27.0),
        ],
    };

    let json = serde_json::to_string_pretty(&dataset).expect("serialize");
    std::fs::write("figures/viability.json", json).expect("write figures/viability.json");
    eprintln!("wrote figures/viability.json");
}

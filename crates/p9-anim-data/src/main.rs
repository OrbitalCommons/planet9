//! Export real, computed values from the reproduction crates to JSON for the
//! manim film, so the visualisations are derived from data rather than
//! hand-tuned constants. Writes `animations/data/anim.json`.

use std::collections::BTreeMap;

use serde::Serialize;
use serde_json::{json, Value};

use p9_core::analysis::circular::{circular_mean, mean_resultant_length, rayleigh_p_value};
use p9_core::analysis::resonance::{critical_perihelion, resonance_semi_major_axis};
use p9_core::analysis::secular::phase_portrait;
use p9_core::analysis::photometry::{planet_apparent_magnitude, ALBEDO_NEPTUNE};
use p9_core::analysis::surveys::SURVEY_DEPTHS;
use p9_core::analysis::thermal::{effective_temp, planck_bnu, solar_equilibrium_temp, C_LIGHT};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::data::stable_kbos::{longitude_of_perihelion, stable_kbos};
use p9_core::types::{helio_distance_at_true_anomaly, P9Params};

use p9_2016_cassini_ranging::perturbation::{favored_true_anomaly, range_perturbation_curve};
use p9_2024_panstarrs::combined_exclusion::compute_combined_from_population;

#[derive(Serialize)]
struct Orbit {
    name: &'static str,
    mass_earth: f64,
    a: f64,
    e: f64,
    i_deg: f64,
    varpi_deg: f64,
    perihelion_au: f64,
    aphelion_au: f64,
    current_distance_au: f64,
}

fn nominal(name: &'static str, p: P9Params) -> Orbit {
    let nu = 118.5_f64.to_radians(); // Cassini-favoured current true anomaly
    Orbit {
        name,
        mass_earth: p.mass_earth,
        a: p.a,
        e: p.e,
        i_deg: p.i.to_degrees(),
        varpi_deg: (p.omega + p.omega_big).to_degrees().rem_euclid(360.0),
        perihelion_au: p.a * (1.0 - p.e),
        aphelion_au: p.a * (1.0 + p.e),
        current_distance_au: helio_distance_at_true_anomaly(&p, nu),
    }
}

/// Real ETNO sample (Brown 2017) with its circular-statistics clustering.
fn etno_sample() -> Value {
    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "i_deg": o.i_deg,
                "varpi_deg": o.longitude_of_perihelion().to_degrees().rem_euclid(360.0),
            })
        })
        .collect();
    let varpis: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| o.longitude_of_perihelion())
        .collect();
    json!({
        "objects": objects,
        "r_bar_varpi": mean_resultant_length(&varpis),
        "rayleigh_p": rayleigh_p_value(&varpis),
        "mean_varpi_deg": circular_mean(&varpis).unwrap_or(0.0).to_degrees(),
    })
}

/// The six stable KBOs (Batygin & Brown 2016) with their clustering.
fn kbo_sample() -> Value {
    let kbos = stable_kbos();
    let objects: Vec<Value> = kbos
        .iter()
        .map(|k| {
            json!({
                "name": k.name,
                "a": k.elements.a,
                "e": k.elements.e,
                "i_deg": k.elements.i.to_degrees(),
                "varpi_deg": longitude_of_perihelion(&k.elements).to_degrees().rem_euclid(360.0),
            })
        })
        .collect();
    let varpis: Vec<f64> = kbos
        .iter()
        .map(|k| longitude_of_perihelion(&k.elements))
        .collect();
    json!({
        "objects": objects,
        "r_bar_varpi": mean_resultant_length(&varpis),
        "rayleigh_p": rayleigh_p_value(&varpis),
    })
}

/// Normalised Planck spectrum of a cold (40 K) body vs wavelength, plus the
/// solar-equilibrium and effective temperatures at 600 AU.
fn thermal() -> Value {
    let t_cold = 40.0;
    let logs: Vec<f64> = (0..=120).map(|k| -1.0 + k as f64 * (4.0 / 120.0)).collect(); // log10(um), 0.1..1000
    let wl_um: Vec<f64> = logs.iter().map(|l| 10f64.powf(*l)).collect();
    let raw: Vec<f64> = wl_um
        .iter()
        .map(|um| planck_bnu(t_cold, C_LIGHT / (um * 1e-6)))
        .collect();
    let peak = raw.iter().cloned().fold(0.0_f64, f64::max).max(1e-300);
    let norm: Vec<f64> = raw.iter().map(|b| b / peak).collect();
    json!({
        "wavelength_um": wl_um,
        "planck_norm_40k": norm,
        "t_eq_600": solar_equilibrium_temp(600.0, 0.41),
        "t_eff_600": effective_temp(600.0, 0.41, 40.0),
        "wien_peak_um": 2.897_771_955e3 / t_cold, // b(µm·K)/T
    })
}

/// Mean-motion-resonance semimajor axes interior to Planet Nine.
fn resonance() -> Value {
    let a9 = P9Params::revised_2019().a; // 500 AU
    // interior p:q (test particle interior -> p < q), real locations via core.
    let pairs = [(2u32, 3u32), (1, 2), (2, 5), (1, 3)];
    let items: Vec<Value> = pairs
        .iter()
        .map(|&(p, q)| {
            json!({
                "label": format!("{}:{}", q, p),
                "a_res": resonance_semi_major_axis(p, q, a9),
            })
        })
        .collect();
    json!({ "a9": a9, "items": items })
}

/// Secular phase portrait H(e, Δϖ) for a test particle under Planet Nine.
fn portrait() -> Value {
    let p9 = P9Params::revised_2019();
    let (e_vals, dvarpi, grid) = phase_portrait(250.0, p9.a, p9.e, p9.gm(), 36, 60);
    json!({
        "e": e_vals,
        "dvarpi_deg": dvarpi.iter().map(|d| d.to_degrees()).collect::<Vec<_>>(),
        "h": grid,
    })
}

/// Combined ZTF + DES + PS1 exclusion, computed per-orbit over the BB21 posterior.
fn exclusion() -> Value {
    let ex = compute_combined_from_population(3000, 2024);
    json!({
        "ztf": ex.ztf_frac,
        "des_unique": ex.des_unique,
        "ps1_unique": ex.ps1_unique,
        "combined": ex.combined,
        "remaining": 1.0 - ex.combined,
    })
}

/// Cassini Earth–Saturn range residual vs Planet Nine true anomaly.
fn cassini() -> Value {
    let p = P9Params::nominal_2016();
    let n = 360;
    let curve = range_perturbation_curve(&p, n);
    json!({
        "nu_deg": curve.true_anomaly.iter().map(|r| r.to_degrees()).collect::<Vec<_>>(),
        "amplitude_km": curve.amplitude_km,
        "favored_nu_deg": favored_true_anomaly(&p, n).to_degrees(),
    })
}

/// Scattered-disk stability boundary: critical perihelion vs semimajor axis.
fn stability() -> Value {
    let a: Vec<f64> = (0..=60).map(|k| 200.0 + k as f64 * (800.0 / 60.0)).collect();
    let q_crit: Vec<f64> = a.iter().map(|&av| critical_perihelion(av)).collect();
    json!({ "a_au": a, "q_crit_au": q_crit })
}

/// Physically-realistic solar-system reference scales for the background scenes:
/// real distances/eccentricities with derived periods & speeds, plus reference
/// magnitudes and temperatures.
fn solar_system() -> Value {
    // (name, a_AU, e, i_deg) — real measured values.
    let bodies: [(&str, f64, f64, f64); 6] = [
        ("Earth", 1.000, 0.017, 0.0),
        ("Jupiter", 5.203, 0.048, 1.3),
        ("Neptune", 30.07, 0.009, 1.8),
        ("Pluto", 39.48, 0.249, 17.2),
        ("Sedna", 506.0, 0.85, 11.9),
        ("Planet Nine (nominal)", 500.0, 0.25, 20.0),
    ];
    let body_json: Vec<Value> = bodies
        .iter()
        .map(|&(name, a, e, i)| {
            json!({
                "name": name,
                "a_au": a,
                "e": e,
                "i_deg": i,
                "q_au": a * (1.0 - e),
                "Q_au": a * (1.0 + e),
                "period_yr": a.powf(1.5),          // Kepler III in (yr, AU)
                "speed_kms": 29.78 / a.sqrt(),       // circular-orbit speed (Earth = 29.78)
            })
        })
        .collect();

    // Apparent magnitudes at opposition: real observed for Neptune/Pluto;
    // computed for Planet Nine from the Neptune-anchored photometry.
    let magnitudes = json!([
        {"name": "Neptune", "dist_au": 30.0, "app_mag": 7.8},     // observed
        {"name": "Pluto", "dist_au": 39.5, "app_mag": 14.4},      // observed
        {"name": "Planet Nine", "dist_au": 600.0,
         "app_mag": planet_apparent_magnitude(6.0, ALBEDO_NEPTUNE, 600.0)},
    ]);

    // Effective temperatures: Neptune observed; Planet Nine computed; CMB fact.
    let temperatures = json!([
        {"name": "Neptune", "t_k": 59.0},
        {"name": "Planet Nine", "t_k": effective_temp(600.0, 0.41, 40.0)},
        {"name": "CMB", "t_k": 2.725},
    ]);

    let depths: Vec<Value> = SURVEY_DEPTHS
        .iter()
        .map(|d| json!({"survey": d.survey, "band": d.band, "limiting_mag": d.limiting_mag}))
        .collect();

    json!({
        "bodies": body_json,
        "magnitudes": magnitudes,
        "temperatures": temperatures,
        "survey_depths": depths,
    })
}

fn main() {
    let mut out: BTreeMap<&str, Value> = BTreeMap::new();
    out.insert(
        "nominal_orbits",
        serde_json::to_value(vec![
            nominal("Batygin & Brown 2016", P9Params::nominal_2016()),
            nominal("Batygin et al. 2019 (review)", P9Params::revised_2019()),
            nominal("Brown & Batygin 2021 (MCMC)", P9Params::mcmc_2021()),
        ])
        .unwrap(),
    );
    out.insert("etno_sample", etno_sample());
    out.insert("kbo_sample", kbo_sample());
    out.insert("thermal", thermal());
    out.insert("resonance", resonance());
    out.insert("phase_portrait", portrait());
    out.insert("exclusion", exclusion());
    out.insert("cassini", cassini());
    out.insert("stability", stability());
    out.insert("solar_system", solar_system());

    let json = serde_json::to_string_pretty(&out).unwrap();
    std::fs::create_dir_all("animations/data").unwrap();
    std::fs::write("animations/data/anim.json", json).unwrap();
    eprintln!("wrote animations/data/anim.json ({} sections)", out.len());
}

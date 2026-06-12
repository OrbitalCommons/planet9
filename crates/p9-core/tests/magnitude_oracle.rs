//! Cross-validation of p9-core's generic-body photometry against
//! starfield's Mallama & Hilton (2018) planetary magnitude model, using
//! Neptune — the closest physical analog to Planet Nine and the anchor of
//! the mass-radius relation — as the common target.
//!
//! The two paths share no constants:
//! - starfield `magnitudelib::planetary_magnitude`: empirical Neptune
//!   model (year-dependent baseline from secular brightening, fitted phase
//!   polynomial above 1.9°), real DE-kernel observer geometry.
//! - p9-core `analysis::photometry`: H from radius + geometric albedo via
//!   the 1329 km diameter relation, the reflected-light 5 log10(rΔ)
//!   distance law, and the IAU H-G phase function (G = 0.15) — evaluated
//!   at the *same* r, Δ, α taken from the observed kernel geometry.
//!
//! Kernel-gated: skips when the starfield cache holds neither de421.bsp
//! nor de440.bsp (never downloads).

use starfield::jplephem::SpiceKernel;
use starfield::jplephem_ext::SpiceKernelExt;
use starfield::magnitudelib::planetary_magnitude;
use starfield::time::Timescale;

use p9_core::analysis::photometry::{
    absolute_magnitude, apparent_magnitude_with_phase, mass_radius_neptunian, phase_angle,
    ALBEDO_NEPTUNE, DEFAULT_SLOPE_G, NEPTUNE_MASS_EARTH,
};
use p9_core::constants::EARTH_RADIUS_KM;

/// Open a JPL ephemeris from the starfield cache (de421 preferred, de440
/// accepted), or `None` when absent so callers can skip hermetically.
fn cached_kernel() -> Option<SpiceKernel> {
    let cache = starfield::data::get_cache_dir();
    for name in ["de421.bsp", "de440.bsp"] {
        let path = cache.join(name);
        let populated = std::fs::metadata(&path)
            .map(|m| m.len() > 0)
            .unwrap_or(false);
        if populated {
            return SpiceKernel::open(&path).ok();
        }
    }
    None
}

/// Neptune's V magnitude computed (a) by starfield's Mallama-Hilton model
/// on real ephemeris geometry and (b) by p9-core's generic photometry
/// (Neptune mass through the mass-radius relation, p_V = 0.41, H-G phase
/// term) at the same r, Δ, α — pinning the two implementations against
/// each other independently of shared constants.
///
/// Agreement is not expected to be tight: Mallama-Hilton carries Neptune's
/// secular brightening (V(0) clamped at −7.00 after ~2000 vs the static
/// H ≈ −6.88 from R = 24 622 km, p_V = 0.41) and no phase term below
/// α = 1.9°, while the H-G law's opposition surge adds ~0.02 mag at
/// Neptune's near-opposition α. Measured at the 2020-09-11 opposition
/// epoch (r = 29.92 AU, Δ = 28.92 AU, α = 0.042°): V_MH = 7.686,
/// V_p9 = 7.828, Δ = +0.142 mag (0.124 from the V(0) baselines, 0.021
/// from the H-G surge) — the documented band below allows ±0.35 mag.
#[test]
fn neptune_v_magnitude_cross_validation() {
    let mut kernel = match cached_kernel() {
        Some(k) => k,
        None => {
            eprintln!("skipping: no de421/de440 kernel in starfield cache");
            return;
        }
    };
    let ts = Timescale::default();
    // Neptune opposition (2020-09-11): minimal phase angle, so the
    // Mallama-Hilton model's missing sub-1.9° phase term contributes the
    // least confound. Inside both DE421 and DE440 coverage.
    let t = ts.utc((2020, 9, 11));

    // (a) Mallama & Hilton 2018 via starfield, real observe() geometry.
    let earth = kernel.at("earth", &t).expect("earth in kernel");
    let neptune = earth
        .observe("neptune barycenter", &mut kernel, &t)
        .expect("neptune in kernel");
    let v_mh = planetary_magnitude(&neptune, &t).expect("Neptune is supported");
    assert!(v_mh.is_finite());

    // The same geometry the magnitude model used: observer->target vector
    // plus the observer's barycentric position (Sun ~ SSB, as in
    // magnitudelib).
    let obs_bary = neptune
        .observer_barycentric
        .as_ref()
        .expect("observe() sets observer_barycentric");
    let sun_to_observer = obs_bary.position;
    let sun_to_neptune = sun_to_observer + neptune.position;
    let r = sun_to_neptune.norm();
    let delta = neptune.position.norm();
    let alpha = phase_angle(r, delta, sun_to_observer.norm());

    // Sanity: opposition geometry (Δ ≈ r − 1, α well below 1.9° so the
    // M-H Neptune model applies no phase correction at all).
    assert!((29.0..31.0).contains(&r), "r = {r:.2} AU");
    assert!((delta - (r - 1.0)).abs() < 0.05, "Δ = {delta:.2} AU");
    assert!(alpha.to_degrees() < 1.9, "α = {:.3}°", alpha.to_degrees());

    // (b) p9-core generic photometry at the same r, Δ, α: Neptune's mass
    // through the (Neptune-anchored) mass-radius relation, V-band
    // geometric albedo 0.41, IAU H-G phase term with the default G.
    let radius_km = mass_radius_neptunian(NEPTUNE_MASS_EARTH) * EARTH_RADIUS_KM;
    let h = absolute_magnitude(radius_km, ALBEDO_NEPTUNE);
    let v_p9 = apparent_magnitude_with_phase(h, r, delta, alpha, DEFAULT_SLOPE_G);

    // Both must land at Neptune's familiar V ≈ 7.7-7.8 near opposition.
    assert!((7.5..8.2).contains(&v_mh), "V_MH = {v_mh:.3}");
    assert!((7.5..8.2).contains(&v_p9), "V_p9 = {v_p9:.3}");

    // Cross-validation band (see doc comment): measured Δ = +0.142 mag at
    // this epoch, dominated by the secular-brightening baseline (−7.00 vs
    // −6.88, 0.124 mag) plus the H-G opposition surge (0.021 mag at this
    // α) that Mallama-Hilton omits below 1.9°.
    let dv = v_p9 - v_mh;
    assert!(
        dv.abs() < 0.35,
        "p9-core vs Mallama-Hilton Neptune V: {v_p9:.3} vs {v_mh:.3} (Δ = {dv:+.3})"
    );
}

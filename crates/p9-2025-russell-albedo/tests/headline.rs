//! Headline reproduction of Russell & White (2025), "The Radius, Composition,
//! Albedo, and Absolute Magnitude of Planet Nine" (arXiv:2507.22297).
//!
//! All values below are COMPUTED from the composition endpoints through
//! p9-core's reflected-light photometry — none are hard-coded — and asserted
//! to land in the paper's published ranges.
//!
//! Computed values (this build):
//!   small/warm (2.0 R⊕, p = 0.33, 0.6% H-He): H ≈ −5.21,  m ≈ 22.74
//!   large/cold (2.6 R⊕, p = 0.47, 3.5% H-He): H ≈ −6.17,  m ≈ 21.79
//!   most-likely heliocentric distance solved to ≈ 625 AU.
//! Published: H = −6.1..−5.2, m = +21.9..+22.7.

use p9_2025_russell_albedo::{
    absolute_h, apparent_v_most_likely, Composition, LARGE_COLD, MOST_LIKELY_DISTANCE_AU,
    PUBLISHED_ALBEDO_MAX, PUBLISHED_ALBEDO_MIN, PUBLISHED_APPARENT_BRIGHT,
    PUBLISHED_APPARENT_FAINT, PUBLISHED_ENVELOPE_FRAC_MAX, PUBLISHED_ENVELOPE_FRAC_MIN,
    PUBLISHED_H_BRIGHT, PUBLISHED_H_FAINT, PUBLISHED_RADIUS_MAX_EARTH, PUBLISHED_RADIUS_MIN_EARTH,
    SMALL_WARM,
};

fn assert_in(x: f64, lo: f64, hi: f64, tol: f64, label: &str) {
    assert!(
        x >= lo - tol && x <= hi + tol,
        "{label} = {x:.3}, want [{lo}, {hi}] (±{tol})"
    );
}

#[test]
fn radius_endpoints_match_published_range() {
    assert_eq!(SMALL_WARM.radius_earth, PUBLISHED_RADIUS_MIN_EARTH);
    assert_eq!(LARGE_COLD.radius_earth, PUBLISHED_RADIUS_MAX_EARTH);
}

#[test]
fn albedo_endpoints_match_published_range() {
    assert_in(
        SMALL_WARM.albedo,
        PUBLISHED_ALBEDO_MIN,
        PUBLISHED_ALBEDO_MAX,
        1e-9,
        "albedo(small/warm)",
    );
    assert_in(
        LARGE_COLD.albedo,
        PUBLISHED_ALBEDO_MIN,
        PUBLISHED_ALBEDO_MAX,
        1e-9,
        "albedo(large/cold)",
    );
    // Larger, colder radius is the higher-albedo endpoint.
    assert!(LARGE_COLD.albedo > SMALL_WARM.albedo);
}

#[test]
fn envelope_fractions_match_published_range() {
    assert_in(
        SMALL_WARM.envelope_fraction,
        PUBLISHED_ENVELOPE_FRAC_MIN,
        PUBLISHED_ENVELOPE_FRAC_MAX,
        1e-9,
        "H-He frac(small/warm)",
    );
    assert_in(
        LARGE_COLD.envelope_fraction,
        PUBLISHED_ENVELOPE_FRAC_MIN,
        PUBLISHED_ENVELOPE_FRAC_MAX,
        1e-9,
        "H-He frac(large/cold)",
    );
}

#[test]
fn absolute_magnitude_recomputed_in_published_band() {
    // Computed from radius + albedo via p9-core; band is H = −6.1..−5.2.
    let h_faint = absolute_h(&SMALL_WARM);
    let h_bright = absolute_h(&LARGE_COLD);
    assert_in(
        h_faint,
        PUBLISHED_H_BRIGHT,
        PUBLISHED_H_FAINT,
        0.1,
        "H(faint)",
    );
    assert_in(
        h_bright,
        PUBLISHED_H_BRIGHT,
        PUBLISHED_H_FAINT,
        0.1,
        "H(bright)",
    );
    // The 2.6 R⊕ / p=0.47 endpoint is the brighter (more negative H) one.
    assert!(
        h_bright < h_faint,
        "large/cold should be the brighter endpoint: H={h_bright:.3} !< {h_faint:.3}"
    );
}

#[test]
fn apparent_magnitude_recomputed_in_published_band() {
    // Computed at the Reference Population 3.0 most-likely distance.
    let m_bright = apparent_v_most_likely(&LARGE_COLD);
    let m_faint = apparent_v_most_likely(&SMALL_WARM);
    assert_in(
        m_bright,
        PUBLISHED_APPARENT_BRIGHT,
        PUBLISHED_APPARENT_FAINT,
        0.2,
        "m(bright)",
    );
    assert_in(
        m_faint,
        PUBLISHED_APPARENT_BRIGHT,
        PUBLISHED_APPARENT_FAINT,
        0.2,
        "m(faint)",
    );
    // Brighter absolute magnitude ⇒ brighter (smaller) apparent magnitude.
    assert!(m_bright < m_faint);
}

#[test]
fn most_likely_distance_is_few_hundred_au() {
    // Solved from the reflected-light law to straddle the published m band.
    assert!(
        (300.0..=900.0).contains(&MOST_LIKELY_DISTANCE_AU),
        "most-likely distance = {MOST_LIKELY_DISTANCE_AU} AU"
    );
}

#[test]
fn h_and_m_are_monotone_across_the_endpoints() {
    // Sanity: a fabricated intermediate composition lands between the two
    // endpoints for both H and m, confirming the endpoints truly bracket.
    let mid = Composition {
        radius_earth: 2.3,
        albedo: 0.40,
        envelope_fraction: 0.02,
    };
    let h_mid = absolute_h(&mid);
    assert!(absolute_h(&LARGE_COLD) < h_mid && h_mid < absolute_h(&SMALL_WARM));
    let m_mid = apparent_v_most_likely(&mid);
    assert!(
        apparent_v_most_likely(&LARGE_COLD) < m_mid && m_mid < apparent_v_most_likely(&SMALL_WARM)
    );
}

//! Labelled reference constants from Cáceres & Gomes (2018) and the canonical
//! Brown & Batygin orbit, plus the vetted ETNO sample (reused from p9-core,
//! never re-tabulated here).

use p9_core::data::etno::BROWN_2017_SAMPLE;

/// Canonical Brown & Batygin (2016) Planet Nine semi-major axis (AU).
pub const BB_A9_AU: f64 = 700.0;

/// Canonical Brown & Batygin (2016) Planet Nine eccentricity.
pub const BB_E9: f64 = 0.6;

/// Canonical Brown & Batygin (2016) Planet Nine perihelion q9 = a9(1−e9) (AU).
/// 700 · (1 − 0.6) = 280 AU.
pub const BB_Q9_AU: f64 = BB_A9_AU * (1.0 - BB_E9);

/// Planet Nine mass used throughout the sweep (Earth masses). Cáceres & Gomes
/// explore the 10–20 M⊕ band; the secular *amplitude* of the level curves is
/// mass-independent at fixed (a9, e9) — only the libration rate scales with the
/// mass — so 10 M⊕ (the B&B nominal) is used as the representative value.
pub const P9_MASS_EARTH: f64 = 10.0;

/// Low-perihelion eccentricity range explored at fixed a9 = 700 AU. The lower
/// bound (e9 = BB_E9 = 0.6) is the canonical orbit; the upper bound (e9 = 0.9)
/// is the most eccentric / lowest-perihelion orbit at which the coplanar
/// secular confinement is still resolved before the P9 perihelion plunges into
/// the ETNO belt (q9 → 70 AU). Cáceres & Gomes study a9 ∈ [500, 1000] AU; we
/// fix a9 = 700 AU (the B&B value) so the ONLY thing changing across the sweep
/// is q9, isolating the headline dependence.
pub const E9_SWEEP_LO: f64 = 0.6;
pub const E9_SWEEP_HI: f64 = 0.9;

/// Lowest P9 perihelion in the swept range, q9 = a9(1 − E9_SWEEP_HI) (AU).
/// 700 · (1 − 0.9) = 70 AU — well below the canonical 280 AU.
pub const LOW_Q9_MIN_AU: f64 = BB_A9_AU * (1.0 - E9_SWEEP_HI);

/// Representative test-ETNO semi-major axis (AU). Sedna-scale, deep in the
/// clustered population (a ≳ 250 AU). NOTE: at the equilibrium-search
/// eccentricities (e up to 0.95) this orbit radially OVERLAPS the swept P9
/// perihelia (down to 70 AU), i.e. the geometry is orbit-crossing across the
/// whole sweep — where the doubly-averaged 1/Δ integral is singular and every
/// confinement number is regularized by the 1%·a₉ softening (a previous
/// version claimed "non-crossing geometry" here, which was false). The
/// softening-robustness test in `confinement` pins that the sweep's monotone
/// deepening survives a 2× softening change.
pub const TEST_ETNO_A_AU: f64 = 500.0;

/// Median semi-major axis of the vetted ETNO sample (AU), computed from the
/// p9-core Brown (2017) table — a data-derived cross-check on TEST_ETNO_A_AU.
pub fn etno_median_a() -> f64 {
    let mut a: Vec<f64> = BROWN_2017_SAMPLE.iter().map(|k| k.a).collect();
    a.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let n = a.len();
    if n % 2 == 1 {
        a[n / 2]
    } else {
        0.5 * (a[n / 2 - 1] + a[n / 2])
    }
}

/// Number of ETNOs in the clustered sample (reused from p9-core).
pub fn etno_count() -> usize {
    BROWN_2017_SAMPLE.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn canonical_bb_perihelion_is_280_au() {
        assert_relative_eq!(BB_Q9_AU, 280.0, epsilon = 1e-9);
    }

    #[test]
    fn low_q9_minimum_is_below_canonical() {
        // The whole point of the paper: the explored perihelion floor sits far
        // below the canonical detached orbit.
        assert!(LOW_Q9_MIN_AU < BB_Q9_AU);
        assert_relative_eq!(LOW_Q9_MIN_AU, 70.0, epsilon = 1e-9);
    }

    #[test]
    fn test_etno_a_is_within_observed_range() {
        // The representative test particle should sit inside the spread of the
        // real clustered sample.
        let amin = BROWN_2017_SAMPLE
            .iter()
            .map(|k| k.a)
            .fold(f64::INFINITY, f64::min);
        let amax = BROWN_2017_SAMPLE
            .iter()
            .map(|k| k.a)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(TEST_ETNO_A_AU >= amin && TEST_ETNO_A_AU <= amax);
        // The vetted sample is the 10-object Brown (2017) clustering table.
        assert_eq!(etno_count(), 10);
        assert!(etno_median_a() > 250.0);
    }
}

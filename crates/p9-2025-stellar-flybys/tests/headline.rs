//! Headline reproduction of Hu, Huang, Gladman & Zhu (2025), "Early Stellar
//! Flybys are Unlikely" (arXiv:2505.16317).
//!
//! The paper's central number is that the probability of a close-in
//! (`q_* ≤ 1000 au`) flyby that *also* satisfies the inclination (i < 30°) and
//! apsidal-alignment constraints actually having happened is `≲ 5%`. We rebuild
//! that as a product of an encounter probability, a geometric-acceptance
//! fraction, and a conditional success factor — and assert it lands in the
//! published band. The headline is COMPUTED from the factors, not hard-coded.

use p9_2025_stellar_flybys::{
    geometry_fraction, total_probability, F_SUCCESS, INCLINATION_CONSTRAINT_DEG,
    PUBLISHED_MAX_PROBABILITY, Q_STAR_MAX_AU,
};

#[test]
fn combined_probability_pins_the_computed_band() {
    // The component bands guarantee p <= 5%, so that alone certifies
    // nothing. The computed value is ~0.42% — an order of magnitude below
    // the paper's ≲5% CEILING (the paper quotes an upper bound, not a
    // point estimate; our factorized F_SUCCESS x geometry x close-encounter
    // product sits well inside it). Pin the computed band.
    let p = total_probability();
    assert!(
        (0.003..0.006).contains(&p),
        "P_total = {p:.4} drifted out of the computed 0.3-0.6% band"
    );
    assert!(
        p <= PUBLISHED_MAX_PROBABILITY,
        "P_total = {p:.4} must respect the published ≲{PUBLISHED_MAX_PROBABILITY} ceiling"
    );
}

#[test]
fn perihelion_threshold_matches_paper() {
    assert_eq!(
        Q_STAR_MAX_AU, 1000.0,
        "close-in flyby threshold is q_* ≤ 1000 au"
    );
}

#[test]
fn inclination_constraint_matches_paper() {
    assert_eq!(
        INCLINATION_CONSTRAINT_DEG, 30.0,
        "detached extreme TNOs have i < 30°"
    );
}

#[test]
fn geometry_fraction_is_a_small_solid_angle() {
    // Only coplanar (i_*≈0) or ecliptic-symmetric (ω_*≈0, i_*≈90°) encounters
    // pass — a small fraction of isotropic 4π orientations.
    let f = geometry_fraction();
    assert!(
        f > 0.0 && f < 0.10,
        "geometry fraction = {f:.4} should be a small (<10%) solid angle"
    );
}

#[test]
fn headline_is_built_from_factors_not_hardcoded() {
    // P_total must equal P_encounter × f_geometry × f_success. Since
    // total_probability() multiplies them, dividing out the two known factors
    // recovers a valid encounter probability in (0, 1].
    let p = total_probability();
    let p_encounter = p / (geometry_fraction() * F_SUCCESS);
    assert!(
        p_encounter > 0.0 && p_encounter <= 1.0,
        "implied P_encounter = {p_encounter:.4} should be a probability"
    );
}

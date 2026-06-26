//! Birth-cluster close-encounter rate for a Sun-mass star.
//!
//! A flyby that can detach sednoids must reach a small perihelion distance
//! `q_* ≤ Q_STAR_MAX_AU` (the paper's close-in threshold, 1000 au). We estimate
//! how often *any* such encounter occurs while the young Sun resides in its
//! birth cluster, using the standard encounter-rate formula
//!
//! ```text
//!   Γ = n_* · σ(q_*) · v_∞
//! ```
//!
//! with a gravitationally focused cross section
//!
//! ```text
//!   σ(q_*) = π q_*² · ( 1 + 2 G (M_⊙ + M_*) / (q_* v_∞²) ) .
//! ```
//!
//! The first term is the geometric "aim a disk of radius q_*" cross section;
//! the bracket is gravitational focusing, which dominates in cold open clusters
//! where the relative velocity dispersion is only ~1 km/s. The expected number
//! of close encounters over the cluster residence time `t` is `N = Γ t`, and the
//! probability of at least one is `1 − exp(−N)`.
//!
//! Cluster parameters are typical of the open clusters in which the Sun likely
//! formed (e.g. Adams 2010, ARA&A 48, 47; Pfalzner 2013). They are deliberately
//! documented constants so the whole bound is transparent and re-derivable.

use p9_core::constants::{GM_SUN, KMS_TO_AUDAY, PC_AU, YEAR_DAYS};
use std::f64::consts::PI;

/// Close-in perihelion threshold for a sednoid-detaching flyby (au).
/// Paper: encounters with `q_* ≤ 1000 au` are the ones that matter.
pub const Q_STAR_MAX_AU: f64 = 1000.0;

/// Typical open-cluster stellar number density (stars per pc³).
/// Solar-birth open clusters are inferred at ~100 stars/pc³
/// (Adams 2010; Pfalzner 2013). Conservative middle-of-the-road value.
pub const CLUSTER_DENSITY_PER_PC3: f64 = 100.0;

/// One-dimensional velocity dispersion of the birth cluster (km/s).
/// Open clusters are dynamically cold: σ_v ~ 1 km/s.
pub const CLUSTER_VELOCITY_DISPERSION_KMS: f64 = 1.0;

/// Cluster residence time of the young Sun (yr).
/// Open clusters disperse on ~1–10 Myr timescales; we take 10 Myr as a
/// generous upper bound so the encounter probability is *not* underestimated.
pub const CLUSTER_RESIDENCE_TIME_YR: f64 = 10.0e6;

/// Mass of the perturbing field star, in solar masses. The cluster mass
/// function is dominated by low-mass stars; ~0.5 M_⊙ is a representative
/// average perturber.
pub const PERTURBER_MASS_SOLAR: f64 = 0.5;

/// Gravitationally focused cross section (au²) for a closest approach within
/// `q_au`, given a relative velocity at infinity `v_inf_kms` (km/s) and a total
/// (Sun + perturber) mass `total_mass_solar` (M_⊙).
pub fn focused_cross_section_au2(q_au: f64, v_inf_kms: f64, total_mass_solar: f64) -> f64 {
    let v_inf = v_inf_kms * KMS_TO_AUDAY; // au/day
                                          // G(M_⊙+M_*) in au³/day² is GM_SUN scaled by the total mass in M_⊙.
    let gm_total = GM_SUN * total_mass_solar;
    let focusing = 1.0 + 2.0 * gm_total / (q_au * v_inf * v_inf);
    PI * q_au * q_au * focusing
}

/// Encounter rate Γ (per day) for closest approaches within `Q_STAR_MAX_AU`.
pub fn encounter_rate_per_day() -> f64 {
    // Number density: stars per pc³ → stars per au³.
    let n_per_au3 = CLUSTER_DENSITY_PER_PC3 / (PC_AU * PC_AU * PC_AU);
    let v_inf_auday = CLUSTER_VELOCITY_DISPERSION_KMS * KMS_TO_AUDAY;
    let total_mass = 1.0 + PERTURBER_MASS_SOLAR;
    let sigma =
        focused_cross_section_au2(Q_STAR_MAX_AU, CLUSTER_VELOCITY_DISPERSION_KMS, total_mass);
    n_per_au3 * sigma * v_inf_auday
}

/// Expected number of close (`q_* ≤ 1000 au`) encounters over the cluster
/// residence time, `N = Γ t`.
pub fn expected_close_encounters() -> f64 {
    let t_days = CLUSTER_RESIDENCE_TIME_YR * YEAR_DAYS;
    encounter_rate_per_day() * t_days
}

/// Probability of **at least one** close encounter during cluster residence,
/// `P = 1 − exp(−N)` (Poisson statistics).
pub fn p_close_encounter() -> f64 {
    1.0 - (-expected_close_encounters()).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn focusing_increases_cross_section() {
        let geom = PI * Q_STAR_MAX_AU * Q_STAR_MAX_AU;
        let focused = focused_cross_section_au2(Q_STAR_MAX_AU, 1.0, 1.5);
        assert!(
            focused > geom,
            "focused {focused:e} should exceed geometric {geom:e}"
        );
    }

    #[test]
    fn cold_cluster_is_focusing_dominated() {
        // At σ_v ~ 1 km/s and q_* = 1000 au, gravitational focusing dominates
        // the geometric term, so the focused/geometric ratio is well above 2.
        let geom = PI * Q_STAR_MAX_AU * Q_STAR_MAX_AU;
        let focused = focused_cross_section_au2(Q_STAR_MAX_AU, 1.0, 1.5);
        assert!(focused / geom > 2.0, "ratio = {:.2}", focused / geom);
    }

    #[test]
    fn close_encounter_is_plausible_in_birth_cluster() {
        // Close (q≤1000 au) encounters are not vanishingly rare in a dense open
        // cluster — the young Sun has an order-tens-of-percent chance of at
        // least one — but they are also far from guaranteed.
        let p = p_close_encounter();
        assert!((0.0..=1.0).contains(&p), "probability out of range: {p}");
        assert!(
            (0.1..0.6).contains(&p),
            "a 1000-au encounter probability should be order tens of percent, got {p:.3}"
        );
    }

    #[test]
    fn rate_is_positive() {
        assert!(encounter_rate_per_day() > 0.0);
        assert!(expected_close_encounters() > 0.0);
    }
}

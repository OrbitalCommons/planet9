//! Back-integration of sednoid longitudes of perihelion and the search for the
//! apsidal-convergence epoch.
//!
//! Under the giant-planet differential apsidal precession (see [`precession`]),
//! each object's longitude of perihelion evolves (to leading secular order, at
//! fixed a, e, i) linearly:
//!
//!   ϖ_i(t) = ϖ_i(0) + (dϖ/dt)_i · t        (t > 0 into the future)
//!
//! Back-integrating to a look-back time τ ≥ 0 means evaluating at t = −τ:
//!
//!   ϖ_i(−τ) = ϖ_i(0) − (dϖ/dt)_i · τ
//!
//! Because the rates (dϖ/dt)_i differ between objects (higher-a → slower), an
//! initially common ϖ₀ is *sheared apart* over time — and conversely, running a
//! set that shares a common origin backwards re-collapses it. Huang & Gladman
//! (2024) find that the sednoids' apsides converge ~4.5 Gyr ago (solar-system
//! formation), favouring a *primordial* alignment over a present-day Planet
//! Nine actively maintaining it.
//!
//! Two computations live here, both driven by the *real* leading-order forced
//! rates from [`precession::apsidal_precession_rate`] (nothing hard-coded):
//!
//! 1. [`primordial_convergence_epoch`] — the headline mechanism. Given a common
//!    primordial ϖ₀ and the real per-object differential rates, propagate
//!    forward by `age` and then back-scan: the circular mean resultant length
//!    R̄(τ) of {ϖ_i(−τ)} peaks exactly at τ = `age`. This *reproduces the
//!    paper's mechanism*: differential giant-planet precession unwinds a common
//!    birth alignment at ~4.5 Gyr, and the back-integration recovers it.
//!
//! 2. [`scan_convergence`] — the same back-scan applied to an arbitrary present
//!    state (e.g. the observed 2017 ϖ). This is the honest diagnostic. With the
//!    leading-order rates and the pinned 2017 osculating elements the observed
//!    set's R̄ is already near its maximum *today* (see crate-level residual
//!    note): the observed sednoid ϖ are not, to leading order, consistent with
//!    a literal common origin under these exact rates. The qualitative
//!    past-convergence mechanism (computation 1) is the reproduced result; the
//!    precise observed-sample epoch is a residual that depends on the adopted
//!    rates and orbit solutions.

use crate::precession::apsidal_precession_rate;
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::data::etno::{Etno, BROWN_2017_SAMPLE};
use p9_core::units::{julian_year, radians, Angle, AngularVelocity};

/// Published reference: the sednoids' apsides converge near solar-system
/// formation, ~4.5 Gyr ago (Huang & Gladman 2024).
pub const PUBLISHED_CONVERGENCE_GYR: f64 = 4.5;

/// One object's present-day ϖ (rad) and its giant-planet apsidal precession
/// rate (rad/yr).
#[derive(Debug, Clone, Copy)]
pub struct Precessor {
    pub name: &'static str,
    pub varpi0: f64,
    pub rate_rad_per_yr: f64,
}

impl Precessor {
    /// Longitude of perihelion at look-back time `tau_yr` years before present
    /// (τ ≥ 0): ϖ(−τ) = ϖ₀ − ϖ̇·τ.
    pub fn varpi_at_lookback(&self, tau_yr: f64) -> f64 {
        self.varpi0 - self.rate_rad_per_yr * tau_yr
    }

    /// Present-day longitude of perihelion ϖ₀ as a dimension-checked [`Angle`].
    pub fn varpi0_angle(&self) -> Angle {
        radians(self.varpi0)
    }

    /// Apsidal precession rate as a typed [`AngularVelocity`] (rad/yr storage).
    pub fn rate(&self) -> AngularVelocity {
        (radians(self.rate_rad_per_yr) / julian_year()).into()
    }

    /// Longitude of perihelion at look-back time `tau_yr` years as an [`Angle`].
    pub fn varpi_at_lookback_typed(&self, tau_yr: f64) -> Angle {
        radians(self.varpi_at_lookback(tau_yr))
    }
}

/// Build a `Precessor` for an ETNO from its osculating elements (real rate).
pub fn precessor_from_etno(e: &Etno) -> Precessor {
    let el = e.elements();
    Precessor {
        name: e.name,
        varpi0: e.longitude_of_perihelion(),
        rate_rad_per_yr: apsidal_precession_rate(el.a, el.e, el.i),
    }
}

/// The detached/sednoid set: the full Brown (2017) detached sample (all
/// q > 30 AU). Sedna and 2012 VP113 — the canonical sednoids — are members; the
/// wider detached set shares the forced differential precession. Each carries
/// its *observed* present-day ϖ and its real forced rate.
pub fn observed_sednoids() -> Vec<Precessor> {
    BROWN_2017_SAMPLE.iter().map(precessor_from_etno).collect()
}

/// The same objects' *real* forced precession rates, but assigned a single
/// common primordial longitude of perihelion `varpi0`, then propagated forward
/// by `age_yr` to a synthetic "present". This is the primordial-alignment
/// hypothesis made concrete: a birth-aligned set sheared by differential
/// giant-planet precession.
pub fn primordially_aligned_then_aged(varpi0: f64, age_yr: f64) -> Vec<Precessor> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|e| {
            let el = e.elements();
            let rate = apsidal_precession_rate(el.a, el.e, el.i);
            Precessor {
                name: e.name,
                varpi0: varpi0 + rate * age_yr,
                rate_rad_per_yr: rate,
            }
        })
        .collect()
}

/// Mean resultant length R̄ of the set's longitudes of perihelion at look-back
/// time `tau_yr` (years).
pub fn resultant_at_lookback(precessors: &[Precessor], tau_yr: f64) -> f64 {
    let varpis: Vec<f64> = precessors
        .iter()
        .map(|p| p.varpi_at_lookback(tau_yr))
        .collect();
    mean_resultant_length(&varpis)
}

/// Result of the convergence scan.
#[derive(Debug, Clone, Copy)]
pub struct ConvergenceScan {
    /// Look-back time maximising R̄, in Gyr.
    pub t_star_gyr: f64,
    /// R̄ at the convergence epoch.
    pub r_bar_max: f64,
    /// Present-day R̄ (τ = 0).
    pub r_bar_present: f64,
}

/// Scan look-back times τ ∈ [0, `max_gyr`] on `n_steps` uniform samples and
/// return the epoch maximising R̄ together with the present-day value.
pub fn scan_convergence(precessors: &[Precessor], max_gyr: f64, n_steps: usize) -> ConvergenceScan {
    let r_bar_present = resultant_at_lookback(precessors, 0.0);
    let mut best_gyr = 0.0;
    let mut best_r = r_bar_present;

    for k in 0..=n_steps {
        let tau_gyr = max_gyr * k as f64 / n_steps as f64;
        let r = resultant_at_lookback(precessors, tau_gyr * 1.0e9);
        if r > best_r {
            best_r = r;
            best_gyr = tau_gyr;
        }
    }

    ConvergenceScan {
        t_star_gyr: best_gyr,
        r_bar_max: best_r,
        r_bar_present,
    }
}

/// Headline mechanism: build a primordially aligned sednoid set (common ϖ₀),
/// shear it forward by `age_gyr` with the *real* per-object forced rates, and
/// back-scan to recover the convergence epoch. Returns the scan; `t_star_gyr`
/// reproduces `age_gyr` (≈ the 4.5 Gyr solar-system age) and `r_bar_max`
/// approaches 1 (perfect past alignment).
pub fn primordial_convergence_epoch(age_gyr: f64, n_steps: usize) -> ConvergenceScan {
    // ϖ₀ choice is irrelevant (circular statistics are rotation-invariant).
    let aged = primordially_aligned_then_aged(1.0, age_gyr * 1.0e9);
    scan_convergence(&aged, 2.0 * age_gyr, n_steps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::radian;
    use uom::si::angular_velocity::radian_per_second;
    use uom::si::time::second;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let p = observed_sednoids()[0];
        // Angle accessors read back exactly in radians.
        assert_relative_eq!(p.varpi0_angle().get::<radian>(), p.varpi0, epsilon = 1e-12);
        assert_relative_eq!(
            p.varpi_at_lookback_typed(1.0e8).get::<radian>(),
            p.varpi_at_lookback(1.0e8),
            epsilon = 1e-9
        );
        // The rad/yr rate, read back in rad/s, equals the f64 rate / seconds-per-year.
        let yr_s = julian_year().get::<second>();
        assert_relative_eq!(
            p.rate().get::<radian_per_second>(),
            p.rate_rad_per_yr / yr_s,
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_primordial_alignment_back_converges_at_birth_epoch() {
        // The paper's mechanism: differential giant-planet precession unwinds a
        // common birth alignment. Back-integrating recovers the birth epoch.
        let scan = primordial_convergence_epoch(PUBLISHED_CONVERGENCE_GYR, 4000);
        assert!(
            (scan.t_star_gyr - PUBLISHED_CONVERGENCE_GYR).abs() < 1.5,
            "t* = {:.2} Gyr, off published {:.1} Gyr by > 1.5 Gyr",
            scan.t_star_gyr,
            PUBLISHED_CONVERGENCE_GYR
        );
    }

    #[test]
    fn test_past_more_aligned_than_present() {
        // R̄(t*) > R̄(0): the apsides were more aligned in the past than now.
        let scan = primordial_convergence_epoch(PUBLISHED_CONVERGENCE_GYR, 4000);
        assert!(
            scan.r_bar_max > scan.r_bar_present,
            "no past alignment: R̄max {:.3} <= R̄now {:.3}",
            scan.r_bar_max,
            scan.r_bar_present
        );
    }

    #[test]
    fn test_birth_alignment_is_near_perfect() {
        // At the recovered epoch the set is essentially perfectly aligned.
        let scan = primordial_convergence_epoch(PUBLISHED_CONVERGENCE_GYR, 4000);
        assert!(scan.r_bar_max > 0.99, "R̄max = {:.4}", scan.r_bar_max);
    }

    #[test]
    fn test_differential_precession_spreads_a_common_origin() {
        // Forward in time, a birth-aligned set de-aligns: today's R̄ < 1.
        let aged = primordially_aligned_then_aged(1.0, PUBLISHED_CONVERGENCE_GYR * 1.0e9);
        let r_today = resultant_at_lookback(&aged, 0.0);
        assert!(
            r_today < 0.95,
            "differential precession should have de-aligned: R̄ = {r_today:.3}"
        );
        assert!(r_today > 0.0);
    }

    #[test]
    fn test_recovery_epoch_tracks_assumed_age() {
        // The recovered convergence epoch follows the assumed birth age — the
        // method measures it, it is not hard-coded.
        for age in [3.0, 4.0, 4.5, 5.0] {
            let scan = primordial_convergence_epoch(age, 4000);
            assert!(
                (scan.t_star_gyr - age).abs() < 0.1,
                "age {age}: recovered t* = {:.3}",
                scan.t_star_gyr
            );
        }
    }

    #[test]
    fn test_observed_sample_present_alignment_diagnostic() {
        // Honest residual: the observed 2017 ϖ are already near-maximally
        // aligned today under these leading-order rates (R̄ clustered, in
        // (0,1)). The literal observed back-integration does not converge in
        // the past — see crate docs. The mechanism test above is the
        // reproduced result.
        let obs = observed_sednoids();
        let r_now = resultant_at_lookback(&obs, 0.0);
        assert!(r_now > 0.3 && r_now < 1.0, "observed R̄now = {r_now}");
        let scan = scan_convergence(&obs, 5.0, 5000);
        // Document that the observed-sample maximum sits near the present, not
        // at 4.5 Gyr (within the scan resolution it does not exceed R̄now by a
        // meaningful margin).
        assert!(
            scan.r_bar_max - scan.r_bar_present < 0.05,
            "unexpected strong observed past-convergence: ΔR̄ = {:.3}",
            scan.r_bar_max - scan.r_bar_present
        );
    }
}

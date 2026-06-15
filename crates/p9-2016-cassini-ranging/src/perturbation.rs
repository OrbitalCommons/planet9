//! Differential acceleration of Planet Nine on Saturn vs the Sun, and its
//! quasi-static conversion to an Earth–Saturn range perturbation.
//!
//! # Physics
//!
//! In a heliocentric frame the relevant quantity is the acceleration P9
//! imparts to Saturn MINUS the acceleration it imparts to the Sun (the frame
//! origin). With `r_p9`, `r_sat` heliocentric positions,
//!
//! ```text
//! a_diff = GM_p9 · [ (r_p9 - r_sat)/|r_p9 - r_sat|³  -  r_p9/|r_p9|³ ].
//! ```
//!
//! The first term is the direct pull on Saturn; the second is the indirect
//! (frame) term — the pull on the Sun. Because P9 is far (`|r_p9| ≫ |r_sat|`)
//! this is the classic tidal field, scaling as `GM_p9 · r_sat / r_p9³`, so it
//! falls steeply with P9 heliocentric distance and varies strongly with ν for
//! the eccentric B&B orbit.
//!
//! The Cassini observable is the Earth–Saturn RANGE ρ, whose line of sight is
//! very nearly the Sun→Saturn direction (Earth's 1 AU offset is small against
//! Saturn's ~9.5 AU). We propagate Saturn along its REAL Cassini-epoch arc
//! ([`saturn_arc`], two-body from the J2000 DE440 state), accumulate the
//! P9-induced Saturn displacement by double time-integration, and project it on
//! the instantaneous line of sight to form the range signature `ρ(t)`
//! ([`range_signature`]). Its RMS about the mean is the PRE-fit residual
//! ([`prefit_amplitude`], the paper's blue curve); the RMS after removing the
//! low-order polynomial the INPOP orbit fit absorbs is the POST-fit residual
//! ([`range_perturbation_amplitude`], the paper's red curve). The favored
//! true-anomaly zone is the post-fit minimum among DETECTABLE positions
//! ([`favored_true_anomaly`]) — the paper's Fig 6 selection logic.

use nalgebra::Vector3;

use p9_core::constants::{AU_KM, GM_SUN, TWO_PI};
use p9_core::initial_conditions::planets::saturn_j2000;
use p9_core::integrator::kepler_step::kepler_drift;
use p9_core::types::P9Params;

use p9_core::types::position_at_true_anomaly as p9_position_at_true_anomaly;

/// Number of time samples along the Cassini Saturn arc.
const SATURN_PHASE_SAMPLES: usize = 64;

/// Cassini Saturn-orbit-insertion epoch (2004-07-01), days after J2000.
const CASSINI_START_DAYS: f64 = 1643.0;

/// Cassini end-of-mission epoch (2017-09-15), days after J2000. The ranging
/// data Fienga et al. fit run to ~2014.4; we use the full mission span as the
/// characteristic arc over which Saturn sweeps real true longitudes.
const CASSINI_END_DAYS: f64 = 6466.0;

/// Saturn's J2000 heliocentric position (AU), DE440 state from `p9_core`; the
/// starting point for the Cassini-arc propagation ([`saturn_arc`]).
pub fn saturn_position() -> Vector3<f64> {
    saturn_j2000().state.pos
}

/// Differential (tidal) acceleration P9 exerts on Saturn relative to the Sun,
/// in AU/day², for P9 at true anomaly `nu` (radians) and Saturn at
/// heliocentric position `r_sat`.
pub fn differential_acceleration_at(
    params: &P9Params,
    nu: f64,
    r_sat: Vector3<f64>,
) -> Vector3<f64> {
    let gm_p9 = params.gm();
    let r_p9 = p9_position_at_true_anomaly(params, nu).position;

    let d = r_p9 - r_sat; // P9 → Saturn separation vector (Sun-frame).
    let d3 = d.norm().powi(3);
    let r3 = r_p9.norm().powi(3);

    gm_p9 * (d / d3 - r_p9 / r3)
}

/// Differential (tidal) acceleration P9 exerts on Saturn relative to the Sun,
/// in AU/day², evaluated at Saturn's J2000 position. Convenience wrapper over
/// [`differential_acceleration_at`].
pub fn differential_acceleration(params: &P9Params, nu: f64) -> Vector3<f64> {
    differential_acceleration_at(params, nu, saturn_position())
}

/// Real Saturn heliocentric positions (AU) and sample times (days after J2000)
/// over the Cassini ranging arc. Saturn's J2000 DE440 state is propagated by
/// two-body Kepler motion ([`kepler_drift`]) from [`CASSINI_START_DAYS`] to
/// [`CASSINI_END_DAYS`], so the sampled geometry follows Saturn's TRUE orbital
/// longitudes during the mission (Saturn moves ~160° over the 13 yr span). This
/// is the geometry the Cassini line of sight actually constrained; using it —
/// rather than a synthetic arc — is what makes the range-residual curve depend
/// on P9's true anomaly through the real Sun–Saturn geometry, as the paper
/// emphasizes ("the perturbation … changes with the geometry of the problem,
/// and not only with the distance").
fn saturn_arc() -> (Vec<f64>, Vec<Vector3<f64>>) {
    let s0 = saturn_j2000().state;
    let n = SATURN_PHASE_SAMPLES;
    let mut times = Vec::with_capacity(n);
    let mut positions = Vec::with_capacity(n);
    for k in 0..n {
        let frac = (k as f64) / ((n - 1) as f64);
        let t = CASSINI_START_DAYS + frac * (CASSINI_END_DAYS - CASSINI_START_DAYS);
        let state = kepler_drift(&s0, t, GM_SUN);
        times.push(t - CASSINI_START_DAYS);
        positions.push(state.pos);
    }
    (times, positions)
}

/// The Earth–Saturn RANGE perturbation signature `ρ(t)` (AU) over the Cassini
/// arc, induced by P9 at true anomaly `nu`.
///
/// P9's own period (~18,500 yr) dwarfs the ~10 yr campaign, so its
/// differential pull on Saturn is treated as quasi-static in P9's position but
/// VARIES IN DIRECTION as SATURN sweeps its arc. Saturn's P9-induced
/// displacement accumulates by double time-integration of the (time-varying)
/// differential acceleration, `Δr(t) = ∬ a_diff(t') dt'`, and the range
/// signature is its projection on the instantaneous (Earth ≈ Sun)→Saturn line
/// of sight, `ρ(t) = Δr(t) · ŝ(t)`. Returns `(times_days, rho_au)`.
fn range_signature(params: &P9Params, nu: f64) -> (Vec<f64>, Vec<f64>) {
    let (times, positions) = saturn_arc();
    let n = positions.len();

    let mut vel = Vector3::zeros(); // accumulated velocity perturbation (AU/day)
    let mut disp = Vector3::zeros(); // accumulated displacement (AU)
    let mut rho = Vec::with_capacity(n);

    let mut t_prev = times[0];
    for (k, r_sat) in positions.iter().enumerate() {
        let dt = times[k] - t_prev;
        t_prev = times[k];
        let a_diff = differential_acceleration_at(params, nu, *r_sat);
        // Leapfrog-style accumulation of the displacement (velocity then drift).
        vel += a_diff * dt;
        disp += vel * dt;
        let los = r_sat.normalize();
        rho.push(disp.dot(&los));
    }
    (times, rho)
}

/// Remove the best-fit polynomial of degree `deg` from `y(t)` and return the
/// residual series. The INPOP global fit re-adjusts Saturn's six osculating
/// elements over the arc, which absorbs the low-order polynomial part of any
/// induced range signature (a constant bias, a linear range rate degenerate
/// with the mean motion, and the quadratic degenerate with the along-track
/// acceleration / semi-major-axis adjustment). Only the higher-order, geometry-
/// driven curvature survives in the POST-fit residuals. We remove up to the
/// quadratic (`deg = 2`).
fn detrend_poly(t: &[f64], y: &[f64], deg: usize) -> Vec<f64> {
    let n = t.len();
    // Normalize time to [-1, 1] for conditioning.
    let t_max = t[n - 1];
    let xs: Vec<f64> = t.iter().map(|&ti| 2.0 * ti / t_max - 1.0).collect();
    let m = deg + 1;

    // Normal equations for least-squares polynomial fit: (VᵀV) c = Vᵀy.
    let mut ata = vec![vec![0.0; m]; m];
    let mut aty = vec![0.0; m];
    for (xi, &yi) in xs.iter().zip(y.iter()) {
        let mut powers = vec![1.0; m];
        for j in 1..m {
            powers[j] = powers[j - 1] * xi;
        }
        for a in 0..m {
            aty[a] += powers[a] * yi;
            for b in 0..m {
                ata[a][b] += powers[a] * powers[b];
            }
        }
    }
    let coeffs = solve_linear(ata, aty);

    xs.iter()
        .zip(y.iter())
        .map(|(&xi, &yi)| {
            let mut fit = 0.0;
            let mut p = 1.0;
            for &c in &coeffs {
                fit += c * p;
                p *= xi;
            }
            yi - fit
        })
        .collect()
}

/// Tiny Gaussian-elimination solver for the small normal-equation systems in
/// [`detrend_poly`] (`m ≤ 3`); avoids pulling in a linear-algebra dependency
/// for a 3×3 solve.
fn solve_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot.
        let mut piv = col;
        for r in (col + 1)..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let d = a[col][col];
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = a[r][col] / d;
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    (0..n).map(|i| b[i] / a[i][i]).collect()
}

/// INPOP nominal post-fit RMS residual for the Cassini data, σ₀ ≈ 74.9 m
/// (Fienga et al. 2016, §4), expressed in km. P9 is only detectable when the
/// range signature it induces is above this precision floor; below it, the
/// position falls in the paper's "uncertainty zone" (Fig 6) where P9 is too
/// faint to be constrained either way.
pub const INPOP_RESIDUAL_FLOOR_KM: f64 = 0.0749;

/// *Pre-fit* Earth–Saturn range perturbation amplitude (km) for P9 at true
/// anomaly `nu`: the RMS, over the Cassini arc, of the raw range signature
/// `ρ(t)` ([`range_signature`]) about its mean. This is the paper's pre-fit
/// residual (Fig 2, blue): it measures how strongly P9 perturbs the Sun–Saturn
/// geometry over the campaign, independent of any orbit re-fit, and is the
/// detectability proxy used by [`favored_true_anomaly`].
pub fn prefit_amplitude(params: &P9Params, nu: f64) -> f64 {
    let (_t, rho) = range_signature(params, nu);
    let n = rho.len() as f64;
    let mean = rho.iter().sum::<f64>() / n;
    let var = rho.iter().map(|r| (r - mean) * (r - mean)).sum::<f64>() / n;
    var.sqrt() * AU_KM
}

/// *Post-fit* Earth–Saturn range residual amplitude (km) for P9 at true anomaly
/// `nu`.
///
/// Builds the range signature `ρ(t)` over the Cassini arc ([`range_signature`])
/// and removes the quadratic polynomial the INPOP orbit fit absorbs
/// ([`detrend_poly`] with `deg = 2`). The RMS of the surviving residual is the
/// post-fit range-residual amplitude. This is the analytic proxy for the
/// paper's post-fit RMS residual `σ(ν)` (Fig 2, red): it is LARGE where P9's
/// pull swings rapidly and un-absorbably across the arc — the perihelion-facing
/// zones the paper excludes — and small over the far half of the orbit.
pub fn range_perturbation_amplitude(params: &P9Params, nu: f64) -> f64 {
    let (times, rho) = range_signature(params, nu);
    let resid = detrend_poly(&times, &rho, 2);
    let n = resid.len() as f64;
    let rms_au = (resid.iter().map(|r| r * r).sum::<f64>() / n).sqrt(); // AU
    rms_au * AU_KM // km
}

/// The favored true anomaly (radians) for P9: the position that MINIMIZES the
/// post-fit range residual ([`range_perturbation_amplitude`]) among the
/// DETECTABLE positions (pre-fit amplitude above [`INPOP_RESIDUAL_FLOOR_KM`])
/// on the outbound, post-perihelion arc `ν ∈ (0, π)`.
///
/// This is the paper's selection logic (Fig 6). The global post-fit minimum
/// lies near aphelion, but there P9 is undetectable (the "uncertainty zone"),
/// so the data say nothing; the perihelion-facing arc is excluded (residuals
/// worsened). What remains is the gap between the two excluded lobes. P9's pull
/// on Saturn is symmetric about the apsidal line, so the detectable edge of
/// that gap appears on BOTH the post-perihelion (ν ≲ 130°) and pre-perihelion
/// (ν ≳ 230°) sides; the Cassini O−C residual breaks the tie in favor of the
/// post-perihelion side (the paper's `ν = 117.8°`, where the induced signature
/// REDUCES rather than worsens the residual). The observed-residual template is
/// not modelled here (see the crate docs); we encode that data-driven choice by
/// scanning the post-perihelion arc `ν ∈ (0, π)`, the half containing the
/// published favored window. `n` is the number of uniform ν samples scanned.
pub fn favored_true_anomaly(params: &P9Params, n: usize) -> f64 {
    let mut best_nu = std::f64::consts::FRAC_PI_2;
    let mut best_resid = f64::INFINITY;
    for k in 1..n {
        let nu = std::f64::consts::PI * (k as f64) / (n as f64); // (0, π)
        if prefit_amplitude(params, nu) < INPOP_RESIDUAL_FLOOR_KM {
            continue; // undetectable: uncertainty zone
        }
        let resid = range_perturbation_amplitude(params, nu);
        if resid < best_resid {
            best_resid = resid;
            best_nu = nu;
        }
    }
    best_nu
}

/// A sampled range-perturbation curve over true anomaly.
#[derive(Debug, Clone)]
pub struct RangePerturbationCurve {
    /// Sample true anomalies (radians), uniform in `[0, 2π)`.
    pub true_anomaly: Vec<f64>,
    /// Range-perturbation amplitude at each sample (km).
    pub amplitude_km: Vec<f64>,
    /// P9 heliocentric distance at each sample (AU).
    pub helio_distance_au: Vec<f64>,
}

impl RangePerturbationCurve {
    /// Index of the minimum-amplitude (best-fit) sample.
    pub fn argmin(&self) -> usize {
        let mut best = 0;
        for i in 1..self.amplitude_km.len() {
            if self.amplitude_km[i] < self.amplitude_km[best] {
                best = i;
            }
        }
        best
    }

    /// True anomaly (radians) that minimizes the range perturbation.
    pub fn argmin_true_anomaly(&self) -> f64 {
        self.true_anomaly[self.argmin()]
    }

    /// True anomaly (degrees) that minimizes the range perturbation.
    pub fn argmin_true_anomaly_deg(&self) -> f64 {
        self.argmin_true_anomaly().to_degrees()
    }

    /// Minimum amplitude (km).
    pub fn min_amplitude_km(&self) -> f64 {
        self.amplitude_km[self.argmin()]
    }

    /// Maximum amplitude (km).
    pub fn max_amplitude_km(&self) -> f64 {
        self.amplitude_km
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

/// Compute the range-perturbation curve over `n` uniform samples of true
/// anomaly in `[0, 2π)`.
pub fn range_perturbation_curve(params: &P9Params, n: usize) -> RangePerturbationCurve {
    assert!(n >= 2, "need at least two samples");
    let mut true_anomaly = Vec::with_capacity(n);
    let mut amplitude_km = Vec::with_capacity(n);
    let mut helio_distance_au = Vec::with_capacity(n);
    for k in 0..n {
        let nu = TWO_PI * (k as f64) / (n as f64);
        true_anomaly.push(nu);
        amplitude_km.push(range_perturbation_amplitude(params, nu));
        helio_distance_au.push(p9_position_at_true_anomaly(params, nu).distance);
    }
    RangePerturbationCurve {
        true_anomaly,
        amplitude_km,
        helio_distance_au,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::data::ephemeris_constraint::{
        brown_batygin_orbit, PREFERRED_TRUE_ANOMALY_DEG, TRUE_ANOMALY_TOLERANCE_DEG,
    };

    #[test]
    fn differential_acceleration_is_tidal_order() {
        // Order-of-magnitude check: |a_diff| ~ GM_p9 · r_sat / r_p9³ at the
        // far part of the orbit (indirect and direct terms nearly cancel).
        let p = brown_batygin_orbit();
        let nu = std::f64::consts::PI; // aphelion, r ≈ a(1+e)
        let a_diff = differential_acceleration(&p, nu).norm();
        let r_p9 = p9_position_at_true_anomaly(&p, nu).distance;
        let r_sat = saturn_position().norm();
        let tidal = p.gm() * r_sat / r_p9.powi(3);
        // Within a factor of a few of the pure-tidal estimate.
        let ratio = a_diff / tidal;
        assert!(
            (0.2..5.0).contains(&ratio),
            "a_diff/tidal = {ratio} (a_diff={a_diff:e}, tidal={tidal:e})"
        );
    }

    #[test]
    fn favored_zone_near_published_true_anomaly() {
        // The detectable post-fit minimum on the outbound arc is the favored
        // zone; it should land within tolerance of the published v = 117.8°.
        let p = brown_batygin_orbit();
        let nu_deg = favored_true_anomaly(&p, 1440).to_degrees();
        let diff = (nu_deg - PREFERRED_TRUE_ANOMALY_DEG).abs();
        assert!(
            diff <= TRUE_ANOMALY_TOLERANCE_DEG,
            "favored ν = {nu_deg:.1}° not within ±{}° of {}°",
            TRUE_ANOMALY_TOLERANCE_DEG,
            PREFERRED_TRUE_ANOMALY_DEG
        );
    }

    #[test]
    fn perihelion_arc_is_excluded_relative_to_favored_zone() {
        // The perihelion-facing arc (paper's forbidden zones) must give a much
        // larger post-fit residual than the favored zone.
        let p = brown_batygin_orbit();
        let nu_fav = favored_true_anomaly(&p, 1440);
        let resid_fav = range_perturbation_amplitude(&p, nu_fav);
        // Sample the perihelion-facing arc (paper's [-65°, 85°]).
        let resid_peri = range_perturbation_amplitude(&p, 20.0_f64.to_radians());
        assert!(
            resid_peri > 5.0 * resid_fav,
            "perihelion arc not excluded: peri={resid_peri:e} vs favored={resid_fav:e}"
        );
    }

    #[test]
    fn postfit_curve_has_two_excluded_lobes() {
        // Fienga Fig 2: the post-fit residual has two forbidden lobes on the
        // perihelion-facing arc, separated by a low gap on the far side.
        let p = brown_batygin_orbit();
        let curve = range_perturbation_curve(&p, 360);
        // The two excluded zones straddle perihelion (ν ≈ 0). Check that the
        // perihelion arc dominates the far-side gap.
        let peri = range_perturbation_amplitude(&p, 10.0_f64.to_radians());
        let gap = range_perturbation_amplitude(&p, 180.0_f64.to_radians());
        assert!(peri > 10.0 * gap, "peri={peri:e} gap={gap:e}");
        // The global maximum sits on the perihelion-facing arc, not the far side.
        let nu_max_deg = {
            let idx = (0..curve.amplitude_km.len())
                .max_by(|&a, &b| {
                    curve.amplitude_km[a]
                        .partial_cmp(&curve.amplitude_km[b])
                        .unwrap()
                })
                .unwrap();
            curve.true_anomaly[idx].to_degrees()
        };
        assert!(
            !(85.0..230.0).contains(&nu_max_deg),
            "post-fit maximum at ν = {nu_max_deg:.0}° should be on the perihelion arc"
        );
    }

    #[test]
    fn prefit_amplitude_minimizes_at_aphelion() {
        // The raw (pre-fit) perturbation is smallest where P9 is farthest:
        // near aphelion (ν = 180°). This is the paper's "uncertainty zone".
        let p = brown_batygin_orbit();
        let pre_apo = prefit_amplitude(&p, std::f64::consts::PI);
        let pre_peri = prefit_amplitude(&p, 20.0_f64.to_radians());
        assert!(
            pre_apo < pre_peri,
            "aphelion pre-fit {pre_apo:e} not below perihelion {pre_peri:e}"
        );
    }

    #[test]
    fn amplitude_scales_linearly_with_mass() {
        // Two masses at the same ν: amplitude ∝ GM_p9 ∝ mass.
        let mut p10 = brown_batygin_orbit();
        p10.mass_earth = 10.0;
        let mut p20 = brown_batygin_orbit();
        p20.mass_earth = 20.0;
        let nu = 20.0_f64.to_radians(); // on the high-signal perihelion arc
        let a10 = range_perturbation_amplitude(&p10, nu);
        let a20 = range_perturbation_amplitude(&p20, nu);
        assert!((a20 / a10 - 2.0).abs() < 1e-9, "ratio = {}", a20 / a10);
    }

    #[test]
    fn amplitude_falls_steeply_with_helio_distance() {
        // Two semi-major axes: a larger orbit puts P9 farther, killing the
        // tidal perturbation. Compare the peak (max) amplitude over the curve.
        let mut p_close = brown_batygin_orbit();
        p_close.a = 400.0;
        let mut p_far = brown_batygin_orbit();
        p_far.a = 800.0;
        let max_close = range_perturbation_curve(&p_close, 360).max_amplitude_km();
        let max_far = range_perturbation_curve(&p_far, 360).max_amplitude_km();
        // Doubling a roughly doubles distance → tidal ~1/r³ → factor several.
        let ratio = max_close / max_far;
        assert!(ratio > 4.0, "distance falloff too weak: ratio = {ratio:.2}");
    }

    #[test]
    fn preferred_zone_sits_near_600_au() {
        // At the published preferred ν = 117.8°, the B&B orbit places P9 near
        // the paper's quoted ~600 AU heliocentric distance.
        let p = brown_batygin_orbit();
        let nu = PREFERRED_TRUE_ANOMALY_DEG.to_radians();
        let r = p9_position_at_true_anomaly(&p, nu).distance;
        // a=700, e=0.6 → r(117.8°) ≈ 620 AU; allow a band around the quote.
        assert!(
            (450.0..800.0).contains(&r),
            "r(117.8°) = {r:.0} AU outside the expected ~600 AU zone"
        );
    }
}

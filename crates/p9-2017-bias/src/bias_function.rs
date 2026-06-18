//! Observational bias function computation.
//!
//! The bias function B(ϖ, Ω, i | a, e, H) gives the relative probability that
//! an object with orbital elements (a, e) and absolute magnitude H would have
//! been discovered with its perihelion at a given orientation on the sky.
//!
//! Key insight from Brown (2017): distant eccentric objects are only bright
//! enough to detect near perihelion, so the bias is dominated by survey
//! coverage near the *perihelion direction*. That coverage is strongly
//! longitude-dependent — surveys avoid the galactic plane, which the ecliptic
//! crosses at ecliptic longitudes λ ≈ 95° and 275° — which is what couples
//! the bias to the longitude of perihelion ϖ.
//!
//! Model components (each a documented assumption standing in for the full
//! MPC pointing/depth catalog used in the paper):
//!
//! 1. **Detectable arc fraction**: the fraction of the orbital period the
//!    object is brighter than the limiting magnitude, using the reflected
//!    sunlight law m = H + 5 log10(r·Δ) (flux ∝ r⁻⁴ at opposition; via
//!    `p9_core::analysis::photometry`). This replaces the previous hard
//!    r_max = 90 AU cutoff and actually uses the catalogued H magnitudes.
//! 2. **Galactic-plane longitude suppression**: a Gaussian suppression of
//!    discovery probability when the perihelion ecliptic longitude falls
//!    near the galactic-plane crossings at λ ≈ 95°/275°, with per-crossing
//!    width and floor. The crossing on the galactic-center side (λ ≈ 275°,
//!    toward Sagittarius) is far more heavily suppressed than the anticenter
//!    crossing (λ ≈ 95°): star density and confusion toward the bulge make
//!    slow-mover searches much shallower over a wider longitude range. This
//!    asymmetry is what gives the null a single dominant lobe in ϖ rather
//!    than two antipodal ones.
//! 3. **Latitude effect**: reduced (not zero) coverage when the perihelion
//!    ecliptic latitude is far off the ecliptic, plus a north/south asymmetry
//!    (the surveys that discovered this sample are predominantly
//!    northern-hemisphere — Palomar, Pan-STARRS, Subaru — so southern
//!    perihelion latitudes are down-weighted).

use std::f64::consts::PI;

use p9_core::analysis::photometry::{apparent_magnitude, opposition_delta};
use p9_core::constants::{DEG2RAD, TWO_PI};
use p9_core::units::{radians, Angle};

/// Ecliptic longitudes (radians) where the galactic plane crosses the
/// ecliptic: λ ≈ 95° (toward the galactic center side) and 275°.
pub const GALACTIC_CROSSING_LON_RAD: [f64; 2] = [95.0 * DEG2RAD, 275.0 * DEG2RAD];

/// Tunable parameters of the simplified bias model. All values are
/// assumptions documented here, not fits to the MPC catalog.
#[derive(Debug, Clone)]
pub struct BiasParams {
    /// Survey limiting magnitude (V). Assumption: the deep wide-field
    /// surveys that discovered the a > 230 AU sample reach m ≈ 24.5 (the
    /// faintest sample member, 2013 RF98, has m ≈ 24.2 at perihelion).
    pub limiting_mag: f64,
    /// 1σ widths (radians) of the Gaussian galactic-plane suppression in
    /// perihelion ecliptic longitude, one per crossing
    /// (`GALACTIC_CROSSING_LON_RAD` order: [anticenter λ ≈ 95°, center
    /// λ ≈ 275°]). Assumption: the dip toward the galactic center is wider
    /// (bulge star fields) than toward the anticenter.
    pub galactic_lon_sigmas: [f64; 2],
    /// Floors of the per-crossing longitude suppression (0 = total
    /// avoidance), same ordering. Assumption: some coverage survives at the
    /// anticenter crossing; essentially none toward the galactic center.
    pub galactic_lon_floors: [f64; 2],
    /// Ecliptic-latitude cutoff (radians) beyond which coverage is reduced.
    pub latitude_cutoff: f64,
    /// Penalty applied outside the latitude cutoff.
    pub latitude_penalty: f64,
    /// Penalty for perihelia at southern ecliptic latitudes (β < 0):
    /// the discovery surveys are predominantly northern-hemisphere.
    pub south_penalty: f64,
    /// Maximum inclination (radians) of the population for the sin(i) prior
    /// marginalization. The observed a > 230 AU sample has i ≲ 30°; the
    /// physical prior range is [0, i_max], not [0, π].
    pub i_max: f64,
}

impl BiasParams {
    /// Ecliptic-latitude cutoff as a dimension-checked [`Angle`].
    pub fn latitude_cutoff_typed(&self) -> Angle {
        radians(self.latitude_cutoff)
    }

    /// Maximum population inclination as a dimension-checked [`Angle`].
    pub fn i_max_typed(&self) -> Angle {
        radians(self.i_max)
    }

    /// Galactic-plane longitude suppression 1σ widths as dimension-checked
    /// [`Angle`]s (same ordering as [`BiasParams::galactic_lon_sigmas`]).
    pub fn galactic_lon_sigmas_typed(&self) -> [Angle; 2] {
        [
            radians(self.galactic_lon_sigmas[0]),
            radians(self.galactic_lon_sigmas[1]),
        ]
    }
}

impl Default for BiasParams {
    fn default() -> Self {
        Self {
            limiting_mag: 24.5,
            galactic_lon_sigmas: [15.0 * DEG2RAD, 40.0 * DEG2RAD],
            galactic_lon_floors: [0.5, 0.02],
            latitude_cutoff: 10.0 * DEG2RAD,
            latitude_penalty: 0.3,
            south_penalty: 0.5,
            i_max: 40.0 * DEG2RAD,
        }
    }
}

/// Fraction of the orbital period during which the object is brighter than
/// `limiting_mag`, given absolute magnitude `h_mag`.
///
/// The limiting heliocentric distance r_lim solves
/// H + 5 log10(r (r − 1)) = m_lim (opposition geometry, Δ = r − 1), and the
/// time fraction inside r_lim follows from Kepler's equation:
/// t(<r)/P = (E − e sin E)/π with cos E = (1 − r/a)/e.
pub fn detectable_fraction(a: f64, e: f64, h_mag: f64, limiting_mag: f64) -> f64 {
    let q = a * (1.0 - e);
    let aphelion = a * (1.0 + e);

    // Bisection for r_lim: m(r) is monotone increasing in r for r > 1.
    let m_at = |r: f64| apparent_magnitude(h_mag, r, opposition_delta(r));
    if m_at(q.max(1.5)) > limiting_mag {
        return 0.0; // too faint even at perihelion
    }
    if m_at(aphelion) <= limiting_mag {
        return 1.0; // detectable along the entire orbit
    }
    let (mut lo, mut hi) = (q.max(1.5), aphelion);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if m_at(mid) <= limiting_mag {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let r_lim = 0.5 * (lo + hi);

    // Kepler time fraction inside r_lim.
    let cos_ea = ((1.0 - r_lim / a) / e).clamp(-1.0, 1.0);
    let ea = cos_ea.acos();
    (ea - e * ea.sin()) / PI
}

/// Ecliptic longitude and latitude (radians) of the perihelion direction for
/// an orbit with argument of perihelion `omega`, node `omega_big` and
/// inclination `i`.
pub fn perihelion_ecliptic_lon_lat(omega: f64, omega_big: f64, i: f64) -> (f64, f64) {
    let beta = (omega.sin() * i.sin()).clamp(-1.0, 1.0).asin();
    let lon = omega_big + (omega.sin() * i.cos()).atan2(omega.cos());
    (lon.rem_euclid(TWO_PI), beta)
}

/// Angular part of the bias (∈ [0, 1]): galactic-plane longitude suppression
/// at the perihelion ecliptic longitude times the ecliptic-latitude coverage
/// penalty. Independent of (a, e, H) — see `perihelion_bias` for the full
/// bias including the brightness factor.
pub fn angular_bias(varpi: f64, omega: f64, i: f64, params: &BiasParams) -> f64 {
    let omega_big = varpi - omega;
    let (lon_peri, beta_peri) = perihelion_ecliptic_lon_lat(omega, omega_big, i);

    // Longitude suppression: Gaussian dips at the two galactic crossings,
    // each with its own width and floor (the galactic-center crossing is
    // much more heavily suppressed than the anticenter one).
    let mut lon_factor = 1.0_f64;
    for (k, &lon_gal) in GALACTIC_CROSSING_LON_RAD.iter().enumerate() {
        let mut d = (lon_peri - lon_gal).rem_euclid(TWO_PI);
        if d > PI {
            d -= TWO_PI;
        }
        let sigma = params.galactic_lon_sigmas[k];
        let floor = params.galactic_lon_floors[k];
        let z = d / sigma;
        let dip = floor + (1.0 - floor) * (1.0 - (-0.5 * z * z).exp());
        lon_factor = lon_factor.min(dip);
    }

    // Latitude effect: the *galactic* latitude is captured by the longitude
    // term; the latitude penalty models reduced wide-field coverage for
    // perihelia far off the ecliptic plane (high |β|), and the south penalty
    // the northern-hemisphere weighting of the discovery surveys.
    let mut lat_factor = if beta_peri.abs() > params.latitude_cutoff {
        params.latitude_penalty
    } else {
        1.0
    };
    if beta_peri < 0.0 {
        lat_factor *= params.south_penalty;
    }

    lon_factor * lat_factor
}

/// Simplified bias function B ∈ [0, 1].
///
/// `varpi` is the longitude of perihelion ϖ = ω + Ω and `omega` the argument
/// of perihelion; the node is recovered as Ω = ϖ − ω. The bias is the product
/// of the detectable arc fraction (brightness) and the angular factor
/// (galactic-plane longitude suppression + latitude penalty).
pub fn perihelion_bias(
    a: f64,
    e: f64,
    h_mag: f64,
    varpi: f64,
    omega: f64,
    i: f64,
    params: &BiasParams,
) -> f64 {
    let brightness = detectable_fraction(a, e, h_mag, params.limiting_mag);
    if brightness <= 0.0 {
        return 0.0;
    }
    brightness * angular_bias(varpi, omega, i, params)
}

/// Compute bias weights for a grid of (ϖ, Ω) values.
///
/// Returns a 2D array of bias values for each (ϖ, Ω) combination,
/// with inclination marginalized over a sin(i) prior restricted to the
/// physical range [0, i_max] of the observed population (not [0, π]).
pub fn bias_grid(
    a: f64,
    e: f64,
    h_mag: f64,
    n_varpi: usize,
    n_omega_big: usize,
    params: &BiasParams,
) -> Vec<Vec<f64>> {
    let d_varpi = TWO_PI / n_varpi as f64;
    let d_omega = TWO_PI / n_omega_big as f64;

    let n_i_samples = 10;

    let mut grid = vec![vec![0.0; n_omega_big]; n_varpi];

    for (j, row) in grid.iter_mut().enumerate() {
        let varpi = (j as f64 + 0.5) * d_varpi;
        for (k, cell) in row.iter_mut().enumerate() {
            let omega_big = (k as f64 + 0.5) * d_omega;
            let omega = varpi - omega_big;

            // Marginalize over inclination with a sin(i) prior on [0, i_max].
            let mut total_weight = 0.0;
            let mut total_sin = 0.0;
            for l in 0..n_i_samples {
                let i = (l as f64 + 0.5) / n_i_samples as f64 * params.i_max;
                let sin_i = i.sin();
                let w = perihelion_bias(a, e, h_mag, varpi, omega, i, params);
                total_weight += w * sin_i;
                total_sin += sin_i;
            }

            *cell = total_weight / total_sin;
        }
    }

    grid
}

/// Maximum of the raw bias over (ϖ, Ω) for a given (a, e, H, i): the true
/// envelope needed for un-clamped rejection sampling.
pub fn bias_max(a: f64, e: f64, h_mag: f64, i: f64, params: &BiasParams) -> f64 {
    let n = 72;
    let mut max = 0.0_f64;
    for j in 0..n {
        let varpi = TWO_PI * j as f64 / n as f64;
        for k in 0..n {
            let omega = TWO_PI * k as f64 / n as f64;
            max = max.max(perihelion_bias(a, e, h_mag, varpi, omega, i, params));
        }
    }
    // Small safety margin for the finite grid resolution; the bias varies on
    // the galactic_lon_sigma ≈ 15° scale, well resolved by 5° sampling.
    max * 1.05
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn typed_bias_param_angles_match_f64() {
        let p = BiasParams::default();
        assert_relative_eq!(
            (p.latitude_cutoff_typed() / radians(1.0)).value,
            p.latitude_cutoff,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (p.i_max_typed() / radians(1.0)).value,
            p.i_max,
            max_relative = 1e-12
        );
        let sig = p.galactic_lon_sigmas_typed();
        assert_relative_eq!(
            (sig[1] / radians(1.0)).value,
            p.galactic_lon_sigmas[1],
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_detectable_fraction_limits() {
        // Sedna (H = 1.6) is detectable near perihelion (q ≈ 76 AU) at
        // m_lim = 24, but not over the whole a = 506 AU orbit.
        let f = detectable_fraction(506.0, 0.85, 1.6, 24.0);
        assert!(f > 0.0 && f < 1.0, "f = {f}");

        // A faint object (H = 8.7) with a distant perihelion is undetectable.
        let f_faint = detectable_fraction(1000.0, 0.85, 8.7, 24.0);
        assert_eq!(f_faint, 0.0);

        // A bright object on a small bright orbit is always detectable.
        let f_bright = detectable_fraction(45.0, 0.1, -2.0, 24.0);
        assert_eq!(f_bright, 1.0);
    }

    #[test]
    fn test_perihelion_bias_depends_on_varpi() {
        // The defining fix: at fixed (a, e, H, ω, i), moving the perihelion
        // longitude onto a galactic crossing must reduce the bias.
        let params = BiasParams::default();
        let (a, e, h, i) = (500.0, 0.85, 4.0, 0.1);
        // ω = 0 puts the perihelion at the node, so λ_peri ≈ Ω = ϖ.
        let b_clear = perihelion_bias(a, e, h, 0.0, 0.0, i, &params);
        let b_center = perihelion_bias(a, e, h, 275.0 * DEG2RAD, 0.0, i, &params);
        let b_anticenter = perihelion_bias(a, e, h, 95.0 * DEG2RAD, 0.0, i, &params);
        assert!(b_clear > 0.0);
        assert!(
            b_center < 0.1 * b_clear,
            "bias at the galactic-center crossing ({b_center:.3}) should be heavily suppressed vs clear sky ({b_clear:.3})"
        );
        // The anticenter crossing is suppressed too, but less deeply.
        assert!(b_anticenter < 0.7 * b_clear && b_anticenter > b_center);
    }

    #[test]
    fn test_perihelion_bias_too_faint() {
        // q = 150 AU with H = 6: far beyond m = 24 at perihelion.
        let b = perihelion_bias(1000.0, 0.85, 6.0, 0.0, 0.0, 0.1, &BiasParams::default());
        assert_eq!(b, 0.0);
    }

    #[test]
    fn test_perihelion_lon_lat() {
        // ω = 90°, i = 30°: perihelion at maximum latitude β = 30°.
        let (_, beta) = perihelion_ecliptic_lon_lat(PI / 2.0, 0.0, 30.0 * DEG2RAD);
        assert!((beta - 30.0 * DEG2RAD).abs() < 1e-12);
        // ω = 0: perihelion on the ecliptic at the node longitude.
        let (lon, beta) = perihelion_ecliptic_lon_lat(0.0, 1.0, 0.5);
        assert!(beta.abs() < 1e-12);
        assert!((lon - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_bias_grid_shape_and_nonuniform() {
        let params = BiasParams::default();
        let grid = bias_grid(300.0, 0.85, 5.0, 36, 36, &params);
        assert_eq!(grid.len(), 36);
        assert_eq!(grid[0].len(), 36);

        let mut min = f64::INFINITY;
        let mut max = 0.0_f64;
        for row in &grid {
            for &val in row {
                assert!(val >= 0.0);
                min = min.min(val);
                max = max.max(val);
            }
        }
        // The longitudinal galactic suppression must make the grid
        // non-uniform in ϖ (the previous model was flat in ϖ).
        assert!(max > 0.0);
        assert!(
            min < 0.8 * max,
            "grid should vary: min {min:.3}, max {max:.3}"
        );
    }

    #[test]
    fn test_bias_max_bounds_bias() {
        let params = BiasParams::default();
        let (a, e, h, i) = (350.0, 0.86, 6.5, 0.3);
        let max = bias_max(a, e, h, i, &params);
        for j in 0..50 {
            let varpi = TWO_PI * (j as f64 + 0.31) / 50.0;
            let omega = TWO_PI * (j as f64 + 0.77) / 31.0;
            let b = perihelion_bias(a, e, h, varpi, omega, i, &params);
            assert!(b <= max, "bias {b} exceeds claimed max {max}");
        }
    }
}

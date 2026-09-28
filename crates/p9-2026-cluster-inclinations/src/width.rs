//! Inclination width of a high-q population, measured the way the paper does
//! (Section 3.1): select on (a, q, i), rotate to the sample's mean orbital
//! pole, and fit the width `w` of
//!
//! ```text
//!   f_t(i) = A · sin i · exp(−i² / 2w²)
//! ```
//!
//! by maximum likelihood over the selected inclinations, with the likelihood
//! normalised over the selection window `i < i_max`.

use nalgebra::{Matrix3, Unit, Vector3};
use p9_core::analysis::poles::{mean_pole_direction, pole_vector};
use p9_core::constants::{DEG2RAD, RAD2DEG};
use p9_core::types::OrbitalElements;

/// Orbital cuts selecting the "high-q TNO" subset of a population.
#[derive(Debug, Clone, Copy)]
pub struct SelectionCuts {
    pub a_min: f64,
    pub a_max: f64,
    pub q_min: f64,
    pub q_max: f64,
    /// Inclination ceiling (radians), applied in the ecliptic frame before
    /// the mean pole is computed.
    pub i_max: f64,
}

impl SelectionCuts {
    /// The paper's inclination-width selection: a ∈ (200, 2000) AU,
    /// q ∈ (40, 80) AU, i < 40°.
    pub fn width_sample() -> Self {
        Self {
            a_min: crate::sample::A_MIN_AU,
            a_max: crate::sample::A_MAX_AU,
            q_min: crate::sample::Q_MIN_AU,
            q_max: crate::sample::Q_MAX_AU,
            i_max: crate::sample::I_MAX_DEG * DEG2RAD,
        }
    }

    /// The whole distant population regardless of perihelion: a ∈ (100,
    /// 5000) AU, i < 40°. Diagnoses whether a change in the windowed width
    /// is dynamical heating/cooling or the flow of orbits through the q
    /// window.
    pub fn whole_population() -> Self {
        Self {
            a_min: 100.0,
            a_max: 5000.0,
            q_min: 0.0,
            q_max: f64::INFINITY,
            i_max: crate::sample::I_MAX_DEG * DEG2RAD,
        }
    }

    /// The paper's perihelion-clustering selection (Section 3.2):
    /// a ≥ 250 AU, q ∈ (40, 100) AU, i ≤ 40°.
    pub fn clustering_sample() -> Self {
        Self {
            a_min: 250.0,
            a_max: f64::INFINITY,
            q_min: 40.0,
            q_max: 100.0,
            i_max: 40.0 * DEG2RAD,
        }
    }

    pub fn passes(&self, e: &OrbitalElements) -> bool {
        let q = e.a * (1.0 - e.e);
        e.e < 1.0
            && e.a > self.a_min
            && e.a < self.a_max
            && q > self.q_min
            && q < self.q_max
            && e.i < self.i_max
    }
}

/// Elements passing the cuts.
pub fn select(elements: &[OrbitalElements], cuts: &SelectionCuts) -> Vec<OrbitalElements> {
    elements
        .iter()
        .copied()
        .filter(|e| cuts.passes(e))
        .collect()
}

/// Mean orbital pole (unit vector, ecliptic frame) of a set of orbits.
pub fn mean_pole(elements: &[OrbitalElements]) -> Vector3<f64> {
    let poles: Vec<[f64; 3]> = elements
        .iter()
        .map(|e| pole_vector(e.i, e.omega_big))
        .collect();
    let p = mean_pole_direction(&poles);
    Vector3::new(p[0], p[1], p[2])
}

/// Rotation taking `pole` onto +z (Rodrigues rotation about `pole × ẑ`).
/// Positions and pole vectors rotated by this matrix are expressed in the
/// frame whose reference plane is the population's mean plane.
pub fn rotation_to_pole(pole: &Vector3<f64>) -> Matrix3<f64> {
    let z = Vector3::z();
    let axis = pole.cross(&z);
    let sin_t = axis.norm();
    let cos_t = pole.dot(&z).clamp(-1.0, 1.0);
    if sin_t < 1e-12 {
        return if cos_t > 0.0 {
            Matrix3::identity()
        } else {
            Matrix3::from_diagonal(&Vector3::new(1.0, -1.0, -1.0))
        };
    }
    let k = Unit::new_normalize(axis);
    let kx = Matrix3::new(0.0, -k.z, k.y, k.z, 0.0, -k.x, -k.y, k.x, 0.0);
    Matrix3::identity() + sin_t * kx + (1.0 - cos_t) * kx * kx
}

/// Inclinations (radians) of each orbit relative to `pole`.
pub fn relative_inclinations(elements: &[OrbitalElements], pole: &Vector3<f64>) -> Vec<f64> {
    elements
        .iter()
        .map(|e| {
            let p = pole_vector(e.i, e.omega_big);
            let p = Vector3::new(p[0], p[1], p[2]);
            p.dot(pole).clamp(-1.0, 1.0).acos()
        })
        .collect()
}

/// Unnormalised intrinsic density sin i · exp(−i²/2w²).
pub fn intrinsic_density(i: f64, w: f64) -> f64 {
    i.sin() * (-0.5 * (i / w).powi(2)).exp()
}

/// Normalisation ∫₀^{i_max} sin i · exp(−i²/2w²) di (Simpson, 400 panels).
pub fn normalisation(w: f64, i_max: f64) -> f64 {
    let n = 400;
    let h = i_max / n as f64;
    let mut s = intrinsic_density(0.0, w) + intrinsic_density(i_max, w);
    for k in 1..n {
        let x = k as f64 * h;
        s += if k % 2 == 1 { 4.0 } else { 2.0 } * intrinsic_density(x, w);
    }
    s * h / 3.0
}

/// Log-likelihood of width `w` for inclinations `incs` (radians) confined to
/// `i < i_max`.
pub fn log_likelihood(incs: &[f64], w: f64, i_max: f64) -> f64 {
    let ln_z = normalisation(w, i_max).ln();
    incs.iter()
        .map(|&i| i.sin().max(1e-300).ln() - 0.5 * (i / w).powi(2) - ln_z)
        .sum()
}

/// Maximum-likelihood width `w` (degrees) on a 0.1° grid over 1°–80°.
pub fn fit_width(incs: &[f64], i_max: f64) -> f64 {
    assert!(!incs.is_empty(), "cannot fit a width to an empty sample");
    let mut best = (f64::NEG_INFINITY, 1.0);
    let mut k = 10;
    while k <= 800 {
        let w_deg = k as f64 / 10.0;
        let ll = log_likelihood(incs, w_deg * DEG2RAD, i_max);
        if ll > best.0 {
            best = (ll, w_deg);
        }
        k += 1;
    }
    best.1
}

/// The paper's width measurement on a population: select, rotate to the
/// mean pole, fit. Returns `(w_deg, n_selected)`; `w_deg` is NaN when fewer
/// than three orbits pass the cuts.
pub fn population_width(elements: &[OrbitalElements], cuts: &SelectionCuts) -> (f64, usize) {
    let sel = select(elements, cuts);
    if sel.len() < 3 {
        return (f64::NAN, sel.len());
    }
    let pole = mean_pole(&sel);
    let incs = relative_inclinations(&sel, &pole);
    let i_max = cuts.i_max;
    let incs: Vec<f64> = incs.into_iter().filter(|&i| i < i_max).collect();
    (fit_width(&incs, i_max), sel.len())
}

/// Convenience: width in degrees of inclinations given in degrees, with the
/// paper's 40° window.
pub fn width_of_degrees(incs_deg: &[f64]) -> f64 {
    let incs: Vec<f64> = incs_deg.iter().map(|d| d * DEG2RAD).collect();
    fit_width(&incs, crate::sample::I_MAX_DEG * DEG2RAD)
}

/// Draw one inclination from f_t(i) on [0, i_max] by rejection (radians).
pub fn draw_inclination<R: rand::Rng>(w: f64, i_max: f64, rng: &mut R) -> f64 {
    // Envelope: the density peaks at or below i = w; bound it by its maximum.
    let mut peak: f64 = 0.0;
    for k in 0..=200 {
        peak = peak.max(intrinsic_density(i_max * k as f64 / 200.0, w));
    }
    loop {
        let i = rng.gen_range(0.0..i_max);
        if rng.gen_range(0.0..peak) < intrinsic_density(i, w) {
            return i;
        }
    }
}

/// Degrees helper for readability in tests and reports.
pub fn deg(x: f64) -> f64 {
    x * RAD2DEG
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn rotation_sends_pole_to_z() {
        let pole = Vector3::new(0.3, -0.4, 0.866).normalize();
        let r = rotation_to_pole(&pole);
        let z = r * pole;
        assert!((z - Vector3::z()).norm() < 1e-12);
        // Proper rotation.
        assert!((r.determinant() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn fit_recovers_generating_width() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let i_max = 40.0 * DEG2RAD;
        for &w_true in &[8.0, 15.0, 25.0] {
            let incs: Vec<f64> = (0..4000)
                .map(|_| draw_inclination(w_true * DEG2RAD, i_max, &mut rng))
                .collect();
            let w = fit_width(&incs, i_max);
            assert!((w - w_true).abs() < 1.5, "w_true = {w_true}, fit = {w}");
        }
    }

    #[test]
    fn relative_inclination_of_pole_itself_is_zero() {
        let e = OrbitalElements {
            a: 300.0,
            e: 0.8,
            i: 0.3,
            omega_big: 1.0,
            omega: 0.0,
            mean_anomaly: 0.0,
        };
        let pole = mean_pole(&[e]);
        assert!(relative_inclinations(&[e], &pole)[0] < 1e-12);
    }
}

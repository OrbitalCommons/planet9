//! Orbital-angle clustering statistics of the DES extreme-TNO subset, under a
//! flat null and under the DES selection function.
//!
//! For each angle (ϖ, ω, Ω) we compute the mean resultant length R̄ and the
//! Rayleigh p-value with `p9_core::analysis::circular`, then a Monte-Carlo
//! p-value comparing the observed R̄ to two nulls:
//!
//!   * `NullModel::Uniform` — isotropic angles, no selection. The "clustering
//!     vs isotropy" baseline.
//!   * `NullModel::DesSelection` — isotropic underlying angles folded through a
//!     documented DES footprint acceptance (perihelion direction must fall in
//!     the ~5000 deg² southern field). This is the paper's actual null.
//!
//! Reproducing Bernardinelli et al. (2021): under the DES selection null, the
//! extreme-TNO orbital angles — in particular ϖ, the Planet 9 alignment angle
//! — are consistent with azimuthal isotropy. A small footprint cannot fully
//! erase the very high Ω concentration of this n = 8 sample (the residual
//! flagged by the 2020 isotropy paper and its reproduction crate); we report Ω
//! honestly rather than forcing it to isotropy.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use p9_core::analysis::circular::{mean_resultant_length, rayleigh_p_value};
use p9_core::constants::{DEG2RAD, TWO_PI};

use crate::catalog::{angle_values, Angle, DesTno};
use p9_core::analysis::surveys::des_footprint_contains;

/// Ecliptic (λ, β) of the perihelion direction for an orbit (i, ω, Ω), radians.
pub fn perihelion_direction(i: f64, omega: f64, node: f64) -> (f64, f64) {
    let beta = (i.sin() * omega.sin()).asin();
    let lambda = (node + (omega.sin() * i.cos()).atan2(omega.cos())).rem_euclid(TWO_PI);
    (lambda, beta)
}

/// Equatorial (RA, Dec) in radians from ecliptic (λ, β), obliquity 23.439°.
pub fn ecliptic_to_equatorial(lambda: f64, beta: f64) -> (f64, f64) {
    let eps = 23.439_291 * DEG2RAD;
    let dec = (beta.sin() * eps.cos() + beta.cos() * eps.sin() * lambda.sin()).asin();
    let ra = (lambda.sin() * eps.cos() - beta.tan() * eps.sin()).atan2(lambda.cos());
    (ra.rem_euclid(TWO_PI), dec)
}

/// Null model for the clustering Monte-Carlo p-values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NullModel {
    /// Isotropic angles, no selection (flat baseline).
    Uniform,
    /// Isotropic angles folded through the DES footprint selection (paper null).
    DesSelection,
}

/// Draw one isotropic (ω, Ω) pair for an inclination i under the chosen null.
fn draw_angles(rng: &mut StdRng, i: f64, null: NullModel) -> (f64, f64) {
    // Hard membership in the shared DES box-union footprint
    // (p9-core::analysis::surveys), attempt-capped like the isotropy crate;
    // the earlier Gaussian-band × 0.25-RA soft model with a 0.02 floor
    // accepted an effective ~10–12 kdeg² — 2–3× the real footprint — and
    // flattened this null toward uniform.
    let mut attempts = 0usize;
    loop {
        let omega = rng.gen::<f64>() * TWO_PI;
        let node = rng.gen::<f64>() * TWO_PI;
        match null {
            NullModel::Uniform => return (omega, node),
            NullModel::DesSelection => {
                attempts += 1;
                if attempts > 100_000 {
                    return (omega, node);
                }
                let (lambda, beta) = perihelion_direction(i, omega, node);
                let (ra, dec) = ecliptic_to_equatorial(lambda, beta);
                if des_footprint_contains(ra.to_degrees(), dec.to_degrees()) {
                    return (omega, node);
                }
            }
        }
    }
}

/// Clustering summary for one angle of a sample.
#[derive(Debug, Clone, Copy)]
pub struct ClusteringResult {
    pub angle: Angle,
    /// Mean resultant length R̄ of the observed angles.
    pub r_bar: f64,
    /// Analytic Rayleigh p (uniformity), small-n corrected.
    pub rayleigh_p: f64,
    /// Monte-Carlo p: fraction of null draws with R̄ ≥ observed.
    pub mc_p: f64,
    pub null: NullModel,
}

impl ClusteringResult {
    /// Consistent with isotropy at the 5% level under this null.
    pub fn consistent_with_isotropy(&self) -> bool {
        self.mc_p > 0.05
    }
}

/// Compute R̄, Rayleigh p, and a selection-aware Monte-Carlo p for one angle.
pub fn clustering_for_angle(
    sample: &[DesTno],
    angle: Angle,
    null: NullModel,
    seed: u64,
    iters: usize,
) -> ClusteringResult {
    let observed = angle_values(sample, angle);
    let observed_r = mean_resultant_length(&observed);
    let rayleigh_p = rayleigh_p_value(&observed);

    let mut rng = StdRng::seed_from_u64(seed);
    let mut ge = 0usize;
    let mut buf = vec![0.0; sample.len()];
    for _ in 0..iters {
        for (k, o) in sample.iter().enumerate() {
            let (omega, node) = draw_angles(&mut rng, o.i_deg * DEG2RAD, null);
            buf[k] = match angle {
                Angle::Node => node,
                Angle::ArgPeri => omega,
                Angle::Varpi => (omega + node).rem_euclid(TWO_PI),
            };
        }
        if mean_resultant_length(&buf) >= observed_r {
            ge += 1;
        }
    }
    let mc_p = (ge as f64 + 1.0) / (iters as f64 + 1.0);

    ClusteringResult {
        angle,
        r_bar: observed_r,
        rayleigh_p,
        mc_p,
        null,
    }
}

/// Run the clustering battery for all three angles under a null.
pub fn clustering_battery(
    sample: &[DesTno],
    null: NullModel,
    seed: u64,
    iters: usize,
) -> Vec<ClusteringResult> {
    Angle::all()
        .iter()
        .enumerate()
        .map(|(k, &angle)| clustering_for_angle(sample, angle, null, seed + k as u64, iters))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::catalog::DES_EXTREME_TNOS;

    const ITERS: usize = 20_000;

    #[test]
    fn varpi_consistent_with_isotropy_under_des_selection() {
        // The headline: the Planet-9 alignment angle ϖ is consistent with
        // azimuthal isotropy once the DES footprint selection is folded in.
        let r = clustering_for_angle(
            &DES_EXTREME_TNOS,
            Angle::Varpi,
            NullModel::DesSelection,
            7,
            ITERS,
        );
        assert!(
            r.consistent_with_isotropy(),
            "varpi mc_p = {:.3} (R-bar = {:.3}) should exceed 0.05 under DES selection",
            r.mc_p,
            r.r_bar
        );
    }

    #[test]
    fn rayleigh_p_is_a_valid_probability() {
        for &angle in &Angle::all() {
            let r = clustering_for_angle(&DES_EXTREME_TNOS, angle, NullModel::Uniform, 1, 2000);
            assert!(
                r.rayleigh_p > 0.0 && r.rayleigh_p <= 1.0,
                "{}: rayleigh p = {}",
                angle.label(),
                r.rayleigh_p
            );
            assert!(r.r_bar >= 0.0 && r.r_bar <= 1.0);
        }
    }

    #[test]
    fn selection_relaxes_clustering_significance() {
        // Folding in the DES footprint raises the MC p-value (orbits look less
        // clustered relative to the selection-shaped null than to a flat null)
        // for ϖ — the angle the paper reports as consistent with isotropy.
        let flat = clustering_for_angle(
            &DES_EXTREME_TNOS,
            Angle::Varpi,
            NullModel::Uniform,
            11,
            ITERS,
        );
        let sel = clustering_for_angle(
            &DES_EXTREME_TNOS,
            Angle::Varpi,
            NullModel::DesSelection,
            11,
            ITERS,
        );
        assert!(
            sel.mc_p >= flat.mc_p,
            "selection p {:.3} should be >= flat p {:.3}",
            sel.mc_p,
            flat.mc_p
        );
    }

    #[test]
    fn node_is_the_labelled_residual() {
        // Documented residual (cf. the 2020 isotropy crate): the n = 8 sample's
        // Omega concentration is extreme; the static footprint cannot fully
        // explain it away. We assert the high observed R-bar rather than
        // forcing Omega to isotropy.
        let r = clustering_for_angle(
            &DES_EXTREME_TNOS,
            Angle::Node,
            NullModel::DesSelection,
            3,
            ITERS,
        );
        assert!(
            r.r_bar > 0.5,
            "Omega R-bar = {:.3} expected to be strongly concentrated",
            r.r_bar
        );
    }
}

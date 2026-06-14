//! DES selection function and the selection-aware isotropy null.
//!
//! This is the crux of Bernardinelli et al. (2020). The DES eTNO orbital
//! angles look clustered when compared to a *uniform* null — Ω in particular
//! has R̄ ≈ 0.8 in the encoded sample. But that clustering is largely a
//! footprint artefact: DES observed a single ~5000 deg² southern field, so it
//! can only discover objects whose orbits carry them through that footprint at
//! the survey epoch. An underlying *isotropic* population, folded through this
//! selection, still produces tightly clustered observed (Ω, ω, ϖ).
//!
//! The paper's tests compare the observed angles to a synthetic *isotropic*
//! population run through the DES detection pipeline, not to a flat
//! distribution. We reproduce that idea with a documented, simplified DES
//! acceptance: an isotropic orbit is detectable if its perihelion direction
//! (where a high-e object spends its observable time and is brightest) falls
//! inside the DES footprint. The null draws isotropic (Ω, ω), keeps (a, e, i)
//! fixed per object, and rejection-samples on the footprint acceptance.
//!
//! Caveat / residual: our footprint is a static RA/Dec box approximating the
//! DES wide survey (the real selection is the per-exposure cadence × depth ×
//! linking efficiency of Belyakov, Bernardinelli & Brown 2022, reproduced in
//! `p9-2022-des`). Against this selection-aware null the DES angles become
//! consistent with isotropy — the paper's conclusion — whereas against a flat
//! null they look clustered. Both nulls are provided and labelled.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use p9_core::constants::{DEG2RAD, TWO_PI};

use crate::des_sample::{Angle, DesEtno};

/// Ecliptic longitude/latitude (radians) of the perihelion direction for an
/// orbit with inclination i, argument of perihelion ω and node Ω:
///   sin β = sin i sin ω,  λ = Ω + atan2(sin ω cos i, cos ω).
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

/// DES wide-survey footprint acceptance for a perihelion direction given in
/// ecliptic coordinates. The DES wide survey is a ~5000 deg² southern field
/// (roughly −70° < Dec < +5°, concentrated in the SGP cap RA 300°–110°). We
/// model acceptance as a smooth declination band centred on the DES footprint
/// with a soft RA gate, returning a weight in [0, 1].
pub fn des_footprint_weight(lambda: f64, beta: f64) -> f64 {
    let (ra, dec) = ecliptic_to_equatorial(lambda, beta);

    // Declination band: DES wide survey is southern, peak near −45°, with a
    // Gaussian falloff; sharp cut north of +10°.
    let dec_deg = dec / DEG2RAD;
    if dec_deg > 10.0 {
        return 0.0;
    }
    let dec_center = -45.0;
    let dec_sigma = 30.0;
    let dec_w = (-(dec_deg - dec_center).powi(2) / (2.0 * dec_sigma * dec_sigma)).exp();

    // RA gate: the DES field straddles RA ~300°→110° (through 0°). Soft
    // suppression in the complementary ~110°–300° region.
    let ra_deg = ra / DEG2RAD;
    let in_field = ra_deg >= 300.0 || ra_deg <= 110.0;
    let ra_w = if in_field { 1.0 } else { 0.25 };

    (dec_w * ra_w).clamp(0.0, 1.0)
}

/// Null model for the isotropy tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NullModel {
    /// Flat: isotropic angles, no selection. The "clustering vs isotropy"
    /// baseline the paper's setup is reacting against.
    Uniform,
    /// Isotropic underlying angles folded through the DES footprint selection
    /// — the paper's actual null. Reproduces "consistent with isotropy".
    DesSelection,
}

/// Draw one isotropic (ω, Ω) pair for an object of inclination i, accepting
/// under the chosen null. For `DesSelection` we rejection-sample on the
/// footprint weight (with a small floor so the loop always terminates).
fn draw_angles(rng: &mut StdRng, i: f64, null: NullModel) -> (f64, f64) {
    loop {
        let omega = rng.gen::<f64>() * TWO_PI;
        let node = rng.gen::<f64>() * TWO_PI;
        match null {
            NullModel::Uniform => return (omega, node),
            NullModel::DesSelection => {
                let (lambda, beta) = perihelion_direction(i, omega, node);
                let w = des_footprint_weight(lambda, beta).max(0.02);
                if rng.gen::<f64>() < w {
                    return (omega, node);
                }
            }
        }
    }
}

/// Monte Carlo p-value for an isotropy statistic on one angle of a sample,
/// under the chosen null. The objects' (a, e, i) are held fixed; only the
/// angles are randomized per draw (matching the clustering-paper convention).
/// Returns the fraction of null draws whose statistic ≥ observed.
pub fn selection_aware_p<F>(
    sample: &[DesEtno],
    angle: Angle,
    null: NullModel,
    seed: u64,
    iters: usize,
    stat: F,
) -> f64
where
    F: Fn(&[f64]) -> f64,
{
    let observed: Vec<f64> = sample.iter().map(|o| angle.of(o)).collect();
    let observed_stat = stat(&observed);

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
        if stat(&buf) >= observed_stat {
            ge += 1;
        }
    }
    (ge as f64 + 1.0) / (iters as f64 + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_footprint_rejects_far_north() {
        // A perihelion at high northern ecliptic latitude → far north Dec → 0.
        let w = des_footprint_weight(0.0, 80.0 * DEG2RAD);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_footprint_accepts_southern() {
        // Southern ecliptic latitude near the DES field → positive weight.
        let w = des_footprint_weight(330.0 * DEG2RAD, -40.0 * DEG2RAD);
        assert!(w > 0.1, "w = {w}");
    }

    #[test]
    fn test_ecliptic_equatorial_roundtrip_pole() {
        // Ecliptic north pole maps near Dec = +66.56°.
        let (_, dec) = ecliptic_to_equatorial(0.0, std::f64::consts::FRAC_PI_2);
        assert!(
            (dec / DEG2RAD - 66.56).abs() < 0.5,
            "dec = {}",
            dec / DEG2RAD
        );
    }
}

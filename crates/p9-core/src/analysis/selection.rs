//! The labeled family of longitude-of-perihelion selection stand-ins.
//!
//! Whether the observed ETNO ϖ clustering is real or a survey artifact is
//! decided by which selection model the null is folded through. The
//! workspace previously scattered three mutually inconsistent stand-ins
//! across crates (a cosine lobe peaked AT the cluster, Gaussian dips at the
//! galactic crossings, and a floored two-crossing suppression); this module
//! is the single home for that family, each member labeled with its origin,
//! and the verdict sensitivity is pinned by a test below: the same data give
//! p ≪ 0.05 under the plane-avoidance members and p ≫ 0.05 under the
//! cluster-aligned lobe. The choice of stand-in IS the verdict, which is why
//! every consumer must say which member it uses.

use crate::constants::{DEG2RAD, TWO_PI};
use std::f64::consts::PI;

/// Ecliptic longitudes where the galactic plane crosses the ecliptic
/// (anticenter λ ≈ 95°, center λ ≈ 275°) — where wide surveys lose depth to
/// stellar crowding.
pub const GALACTIC_CROSSING_LON_DEG: [f64; 2] = [95.0, 275.0];

/// One member of the ϖ-selection stand-in family. All weights are relative
/// detection probabilities in [0, 1]-ish scale (unnormalized), evaluated at
/// a longitude of perihelion (radians).
#[derive(Debug, Clone, Copy)]
pub enum VarpiSelection {
    /// Two-harmonic cosine lobe peaked at `phi1` (radians) — the
    /// p9-2021-napier-critique composite stand-in for the OSSOS/DES/S&T
    /// pointing histories. With the lobe placed on the observed cluster
    /// direction, isotropic draws routinely reproduce the observed R̄: the
    /// selection-artifact mechanism.
    CosineLobes {
        a1: f64,
        phi1: f64,
        a2: f64,
        phi2: f64,
    },
    /// Gaussian suppression dips of equal `depth` and width `sigma`
    /// (radians) at the two galactic crossings — the p9-2021-orbit
    /// survey-bias null.
    CrossingDips { depth: f64, sigma: f64 },
    /// Floored, per-crossing suppression (min over crossings of
    /// `floor + (1 − floor)(1 − exp(−Δ²/2σ²))`) — the p9-2017-bias
    /// longitude factor, with the galactic-center crossing much more
    /// heavily suppressed than the anticenter one. Order matches
    /// [`GALACTIC_CROSSING_LON_DEG`].
    FlooredCrossingDips { sigmas: [f64; 2], floors: [f64; 2] },
}

/// Signed smallest angular difference a − b on the circle (radians, in
/// (−π, π]).
fn wrap_diff(a: f64, b: f64) -> f64 {
    let mut d = (a - b).rem_euclid(TWO_PI);
    if d > PI {
        d -= TWO_PI;
    }
    d
}

impl VarpiSelection {
    /// Relative detection weight at longitude of perihelion `varpi`
    /// (radians). Non-negative.
    pub fn weight(&self, varpi: f64) -> f64 {
        match *self {
            VarpiSelection::CosineLobes { a1, phi1, a2, phi2 } => {
                (1.0 + a1 * (varpi - phi1).cos() + a2 * (2.0 * (varpi - phi2)).cos()).max(0.0)
            }
            VarpiSelection::CrossingDips { depth, sigma } => {
                let mut w = 1.0;
                for lon_deg in GALACTIC_CROSSING_LON_DEG {
                    let d = wrap_diff(varpi, lon_deg * DEG2RAD);
                    w -= depth * (-d * d / (2.0 * sigma * sigma)).exp();
                }
                w.max(0.0)
            }
            VarpiSelection::FlooredCrossingDips { sigmas, floors } => {
                let mut w = 1.0_f64;
                for (k, lon_deg) in GALACTIC_CROSSING_LON_DEG.iter().enumerate() {
                    let z = wrap_diff(varpi, lon_deg * DEG2RAD) / sigmas[k];
                    let dip = floors[k] + (1.0 - floors[k]) * (1.0 - (-0.5 * z * z).exp());
                    w = w.min(dip);
                }
                w
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::circular::mean_resultant_length;
    use crate::data::etno::longitudes_of_perihelion;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    /// Draw one ϖ from the density ∝ w(ϖ) by rejection.
    fn draw_varpi<R: Rng>(sel: &VarpiSelection, w_max: f64, rng: &mut R) -> f64 {
        loop {
            let v = rng.gen::<f64>() * TWO_PI;
            if rng.gen::<f64>() * w_max < sel.weight(v) {
                return v;
            }
        }
    }

    /// Selection-aware MC p-value: fraction of null resamples (n objects
    /// drawn from the selection density) whose R̄ meets the observed one.
    fn selection_p(sel: &VarpiSelection, varpi_obs: &[f64], seed: u64, iters: usize) -> f64 {
        let r_obs = mean_resultant_length(varpi_obs);
        let w_max = (0..3600)
            .map(|k| sel.weight(k as f64 * TWO_PI / 3600.0))
            .fold(0.0_f64, f64::max);
        let mut rng = StdRng::seed_from_u64(seed);
        let mut hits = 0usize;
        for _ in 0..iters {
            let draw: Vec<f64> = (0..varpi_obs.len())
                .map(|_| draw_varpi(sel, w_max, &mut rng))
                .collect();
            if mean_resultant_length(&draw) >= r_obs {
                hits += 1;
            }
        }
        hits as f64 / iters as f64
    }

    #[test]
    fn weights_are_nonnegative_and_dip_where_advertised() {
        let members = [
            VarpiSelection::CosineLobes {
                a1: 0.90,
                phi1: 52.0 * DEG2RAD,
                a2: 0.09,
                phi2: 52.0 * DEG2RAD,
            },
            VarpiSelection::CrossingDips {
                depth: 0.85,
                sigma: 20.0 * DEG2RAD,
            },
            VarpiSelection::FlooredCrossingDips {
                sigmas: [15.0 * DEG2RAD, 40.0 * DEG2RAD],
                floors: [0.5, 0.02],
            },
        ];
        for sel in &members {
            for k in 0..720 {
                assert!(sel.weight(k as f64 * TWO_PI / 720.0) >= 0.0);
            }
        }
        // The lobe peaks at its phase; the dip members bottom out at a
        // crossing.
        assert!(members[0].weight(52.0 * DEG2RAD) > members[0].weight(232.0 * DEG2RAD));
        assert!(members[1].weight(95.0 * DEG2RAD) < 0.3);
        assert!(members[2].weight(275.0 * DEG2RAD) < 0.05);
    }

    #[test]
    fn verdict_depends_on_the_selection_stand_in() {
        // THE point of housing the family in one place: on the same vetted
        // ETNO ϖ sample, the plane-avoidance members leave the clustering
        // significant, while the cluster-aligned lobe dissolves it. Any
        // consumer quoting a selection-debiased p must therefore name its
        // member; the numbers below pin the spread.
        let varpi = longitudes_of_perihelion();

        let lobe = VarpiSelection::CosineLobes {
            a1: 0.90,
            phi1: 52.0 * DEG2RAD,
            a2: 0.09,
            phi2: 52.0 * DEG2RAD,
        };
        let dips = VarpiSelection::CrossingDips {
            depth: 0.85,
            sigma: 20.0 * DEG2RAD,
        };
        let floored = VarpiSelection::FlooredCrossingDips {
            sigmas: [15.0 * DEG2RAD, 40.0 * DEG2RAD],
            floors: [0.5, 0.02],
        };

        let p_lobe = selection_p(&lobe, &varpi, 11, 4000);
        let p_dips = selection_p(&dips, &varpi, 13, 4000);
        let p_floored = selection_p(&floored, &varpi, 17, 4000);

        assert!(p_lobe > 0.10, "cluster-aligned lobe: p = {p_lobe}");
        assert!(p_dips < 0.05, "crossing dips: p = {p_dips}");
        assert!(p_floored < 0.05, "floored dips: p = {p_floored}");
    }
}

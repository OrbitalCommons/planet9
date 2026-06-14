//! Intrinsic (isotropic) ETNO population generator.
//!
//! Shankman et al. (2017) test whether the apparent clustering in longitude
//! of perihelion ϖ of the large-a TNO sample can be produced by detection
//! bias acting on an *intrinsically isotropic* parent population. The null
//! population therefore draws ϖ uniformly on [0, 2π) — every orbital
//! orientation equally likely — so its intrinsic mean resultant length R̄ is
//! ~0 (Rayleigh-consistent with uniformity) before any survey selection.
//!
//! The orbital-element distributions follow the paper's detached/scattering
//! large-a regime:
//!
//! - **a** uniform in [150, 800] AU (the "large semimajor axis" range of the
//!   title; the striking biases set in for a ≳ 150 AU where objects are only
//!   detectable near perihelion).
//! - **q** (perihelion) uniform in a detached range [30, 80] AU, with
//!   e = 1 − q/a; objects with q below Neptune are scattered too quickly to
//!   matter for the apsidally-confined sample.
//! - **i** half-normal with σ ≈ 16° (the dynamically excited but not
//!   isotropic inclination distribution of the scattering disk).
//! - **ϖ = Ω + ω** drawn uniformly (isotropic), with ω uniform and the node
//!   Ω = ϖ − ω; the mean anomaly is uniform so objects sit at a random orbital
//!   phase.
//!
//! Each generated object carries an absolute magnitude H drawn from an
//! exponential (the standard log-uniform TNO size distribution truncated to
//! the bright end actually reachable by wide-field surveys).

use p9_core::constants::{DEG2RAD, TWO_PI};
use p9_core::types::OrbitalElements;
use rand::Rng;
use rand_distr::{Distribution, Exp};

/// One synthetic ETNO: orbital elements plus absolute magnitude.
#[derive(Debug, Clone, Copy)]
pub struct SyntheticEtno {
    pub elements: OrbitalElements,
    /// Absolute magnitude H.
    pub h_mag: f64,
}

impl SyntheticEtno {
    /// Longitude of perihelion ϖ = Ω + ω, wrapped to [0, 2π).
    pub fn longitude_of_perihelion(&self) -> f64 {
        (self.elements.omega_big + self.elements.omega).rem_euclid(TWO_PI)
    }
}

/// Parameters of the intrinsic population.
#[derive(Debug, Clone)]
pub struct PopulationParams {
    /// Semi-major-axis range (AU).
    pub a_range: (f64, f64),
    /// Perihelion-distance range (AU) of the detached population.
    pub q_range: (f64, f64),
    /// 1σ of the half-normal inclination distribution (radians).
    pub i_sigma: f64,
    /// Mean of the exponential H distribution above `h_min`.
    pub h_scale: f64,
    /// Bright-end cutoff of the H distribution.
    pub h_min: f64,
}

impl Default for PopulationParams {
    fn default() -> Self {
        Self {
            a_range: (150.0, 800.0),
            q_range: (30.0, 80.0),
            i_sigma: 16.0 * DEG2RAD,
            h_scale: 1.2,
            h_min: 4.0,
        }
    }
}

/// Generate `n` synthetic ETNOs with the given RNG. Longitude of perihelion
/// is uniform (isotropic), so the intrinsic ϖ sample is statistically
/// indistinguishable from a uniform circular distribution.
pub fn generate_population<R: Rng>(
    n: usize,
    params: &PopulationParams,
    rng: &mut R,
) -> Vec<SyntheticEtno> {
    let exp = Exp::new(1.0 / params.h_scale).expect("positive scale");
    let mut out = Vec::with_capacity(n);
    for _ in 0..n {
        let a = rng.gen_range(params.a_range.0..params.a_range.1);
        // Perihelion confined to the detached range, but never above aphelion.
        let q_hi = params.q_range.1.min(a);
        let q = rng.gen_range(params.q_range.0..q_hi);
        let e = (1.0 - q / a).clamp(0.0, 0.999);

        // Half-normal inclination: |N(0, σ)|.
        let i = (rng.sample::<f64, _>(rand_distr::StandardNormal) * params.i_sigma).abs();

        // ISOTROPIC longitude of perihelion: ϖ uniform, ω uniform, Ω = ϖ − ω.
        let varpi = rng.gen::<f64>() * TWO_PI;
        let omega = rng.gen::<f64>() * TWO_PI;
        let omega_big = (varpi - omega).rem_euclid(TWO_PI);

        let mean_anomaly = rng.gen::<f64>() * TWO_PI;

        let h_mag = params.h_min + exp.sample(rng);

        out.push(SyntheticEtno {
            elements: OrbitalElements {
                a,
                e,
                i,
                omega_big,
                omega,
                mean_anomaly,
            },
            h_mag,
        });
    }
    out
}

/// Longitudes of perihelion of a population (radians).
pub fn longitudes_of_perihelion(pop: &[SyntheticEtno]) -> Vec<f64> {
    pop.iter().map(|o| o.longitude_of_perihelion()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::analysis::circular::{mean_resultant_length, rayleigh_p_value};
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_population_within_ranges() {
        let mut rng = StdRng::seed_from_u64(7);
        let p = PopulationParams::default();
        let pop = generate_population(5000, &p, &mut rng);
        for o in &pop {
            assert!(o.elements.a >= p.a_range.0 && o.elements.a <= p.a_range.1);
            let q = o.elements.perihelion();
            assert!(
                q >= p.q_range.0 - 1e-9 && q <= p.q_range.1 + 1e-9,
                "q = {q}"
            );
            assert!(o.elements.i >= 0.0);
            assert!(o.h_mag >= p.h_min);
        }
    }

    #[test]
    fn test_intrinsic_varpi_is_isotropic() {
        // The whole point of the null: uniform ϖ → R̄ ≈ 0. For a truly
        // isotropic sample R̄ ≈ 1/√n (here ≈ 0.007 at n = 20 000), so a tight
        // R̄ < 0.05 bound holds for every seed.
        //
        // The Rayleigh p-value is, by construction, uniform on [0, 1] under the
        // null — so ~5% of seeds WILL fall below 0.05 purely by chance; that is
        // not evidence of clustering. We therefore assert the median behaviour
        // across seeds (most seeds non-significant) rather than demanding every
        // seed clear 0.05, which would be a statistically incorrect test.
        let mut n_significant = 0;
        for seed in [1u64, 2, 3, 42, 100, 7, 13, 99] {
            let mut rng = StdRng::seed_from_u64(seed);
            let pop = generate_population(20_000, &PopulationParams::default(), &mut rng);
            let varpis = longitudes_of_perihelion(&pop);
            let r_bar = mean_resultant_length(&varpis);
            let p = rayleigh_p_value(&varpis);
            assert!(r_bar < 0.05, "seed {seed}: intrinsic R̄ = {r_bar:.4} not ~0");
            // No seed should look strongly clustered.
            assert!(
                p > 1e-3,
                "seed {seed}: intrinsic Rayleigh p = {p:.4} far too small"
            );
            if p < 0.05 {
                n_significant += 1;
            }
        }
        // Most seeds non-significant (a handful near the 5% tail is expected).
        assert!(
            n_significant <= 2,
            "{n_significant}/8 isotropic seeds flagged significant — too many"
        );
    }
}

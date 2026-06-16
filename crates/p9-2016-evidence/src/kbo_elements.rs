//! Clone-stability screening and the headline joint clustering statistic for
//! the 6 dynamically stable KBOs from Batygin & Brown (2016).
//!
//! The KBO table itself (`KboRecord`, `stable_kbos`, `longitude_of_perihelion`)
//! lives in [`p9_core::data::stable_kbos`]. Clone generation uses the documented
//! assumed dispersions below.

use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::constants::*;
use p9_core::data::stable_kbos::{longitude_of_perihelion, stable_kbos, KboRecord};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets::neptune_only_j2000;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{cartesian_to_elements, OrbitalElements, SimConfig};
use p9_core::units::{days, Time};

/// Assumed 1σ dispersions for clone generation. These are *assumptions*
/// standing in for the per-object orbit-fit covariances (not transcribed
/// here), sized to typical multi-opposition ETNO fit uncertainties.
pub const CLONE_SIGMA_A_FRAC: f64 = 0.005;
pub const CLONE_SIGMA_E: f64 = 0.002;
pub const CLONE_SIGMA_I_DEG: f64 = 0.01;
pub const CLONE_SIGMA_ANGLE_DEG: f64 = 0.05;

/// Generate `n_clones` of a KBO, perturbed within the assumed dispersions
/// above.
pub fn generate_clones(
    kbo: &KboRecord,
    n_clones: usize,
    rng: &mut impl rand::Rng,
) -> Vec<OrbitalElements> {
    use rand_distr::{Distribution, Normal};

    let a_norm = Normal::new(0.0, CLONE_SIGMA_A_FRAC * kbo.elements.a).unwrap();
    let e_norm = Normal::new(0.0, CLONE_SIGMA_E).unwrap();
    let i_norm = Normal::new(0.0, CLONE_SIGMA_I_DEG * DEG2RAD).unwrap();
    let angle_norm = Normal::new(0.0, CLONE_SIGMA_ANGLE_DEG * DEG2RAD).unwrap();

    (0..n_clones)
        .map(|_| {
            let mut elem = kbo.elements;
            elem.a += a_norm.sample(rng);
            elem.e = (elem.e + e_norm.sample(rng)).clamp(0.0, 0.9999);
            elem.i += i_norm.sample(rng);
            elem.omega += angle_norm.sample(rng);
            elem.omega_big += angle_norm.sample(rng);
            elem.mean_anomaly += angle_norm.sample(rng);
            elem
        })
        .collect()
}

/// Result of clone-stability screening for one object.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CloneStabilityResult {
    pub name: String,
    pub n_clones: usize,
    /// Clones that stayed bound, were never removed, and kept |Δa|/a below
    /// `STABILITY_MAX_DA_FRAC`.
    pub n_stable: usize,
    /// Largest |Δa|/a over the surviving clones.
    pub max_da_frac: f64,
    /// RNG seed used (recorded for reproducibility).
    pub seed: u64,
    /// Integration horizon (days).
    pub t_total: f64,
}

impl CloneStabilityResult {
    /// Integration horizon as a typed [`Time`] (days).
    pub fn integration_time(&self) -> Time {
        days(self.t_total)
    }
}

/// Stability screen: a clone counts as stable when its semi-major axis
/// changes by less than this fraction over the horizon (semi-major-axis
/// diffusion is the instability channel for this sample).
pub const STABILITY_MAX_DA_FRAC: f64 = 0.10;

/// Screen one KBO's clones for dynamical stability by direct integration:
/// Neptune integrated directly, Jupiter+Saturn+Uranus absorbed into the
/// J2-averaged solar quadrupole (`ExtraForce::J2Jsu`), no Planet Nine.
///
/// Reduced-scale by default (call with t_total ≪ 4 Gyr); the 4-Gyr
/// full-horizon screen lives behind an `#[ignore]`d test. `dt` must resolve
/// Neptune (≤ P_N/20 ≈ 3000 d).
pub fn clone_stability_screen(
    kbo: &KboRecord,
    n_clones: usize,
    t_total: f64,
    dt: f64,
    seed: u64,
) -> CloneStabilityResult {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let clones = generate_clones(kbo, n_clones, &mut rng);
    let a0: Vec<f64> = clones.iter().map(|c| c.a).collect();
    let mut particles: Vec<_> = clones.iter().map(|c| c.to_state_vector(GM_SUN)).collect();
    let mut active = vec![true; particles.len()];

    let mut bodies = neptune_only_j2000();
    let integrator = WhmIntegrator::with_extra_forces(vec![ExtraForce::J2Jsu]);

    let config = SimConfig {
        dt,
        t_start: 0.0,
        t_end: t_total,
        removal_inner_au: 5.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: t_total,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    let n_steps = (t_total / dt).ceil() as usize;
    for _ in 0..n_steps {
        integrator.step(&mut bodies, &mut particles, &mut active, dt, &config);
    }

    let mut n_stable = 0;
    let mut max_da_frac: f64 = 0.0;
    for (k, particle) in particles.iter().enumerate() {
        if !active[k] {
            continue;
        }
        let elem = cartesian_to_elements(particle, GM_SUN);
        if elem.e >= 1.0 || elem.a <= 0.0 {
            continue;
        }
        let da_frac = ((elem.a - a0[k]) / a0[k]).abs();
        max_da_frac = max_da_frac.max(da_frac);
        if da_frac < STABILITY_MAX_DA_FRAC {
            n_stable += 1;
        }
    }

    CloneStabilityResult {
        name: kbo.name.to_string(),
        n_clones,
        n_stable,
        max_da_frac,
        seed,
        t_total,
    }
}

/// Orbital pole unit vector in ecliptic coordinates.
fn pole_vector(elem: &OrbitalElements) -> [f64; 3] {
    let (si, ci) = (elem.i.sin(), elem.i.cos());
    let (so, co) = (elem.omega_big.sin(), elem.omega_big.cos());
    [si * so, -si * co, ci]
}

/// Mean resultant length of a set of pole unit vectors: |Σ p̂|/n.
fn pole_resultant_length(poles: &[[f64; 3]]) -> f64 {
    let mut s = [0.0; 3];
    for p in poles {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt() / poles.len() as f64
}

/// Observed clustering statistics from the 6-KBO sample, computed with
/// circular (not arithmetic) means: (mean ϖ, mean Ω, mean ω) in radians.
pub fn observed_clustering_stats() -> (f64, f64, f64) {
    let kbos = stable_kbos();
    let varpis: Vec<f64> = kbos
        .iter()
        .map(|k| longitude_of_perihelion(&k.elements))
        .collect();
    let omegas: Vec<f64> = kbos.iter().map(|k| k.elements.omega).collect();
    let nodes: Vec<f64> = kbos.iter().map(|k| k.elements.omega_big).collect();

    (
        circular_mean(&varpis).expect("clustered sample has a mean direction"),
        circular_mean(&nodes).expect("clustered sample has a mean direction"),
        circular_mean(&omegas).expect("clustered sample has a mean direction"),
    )
}

/// Result of the joint ϖ + orbital-pole clustering Monte Carlo — the
/// paper's headline statistic (quoted as P ≈ 0.007%).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct JointClusteringResult {
    /// Observed mean resultant length of the 6 longitudes of perihelion.
    pub r_bar_varpi_obs: f64,
    /// Observed mean resultant length of the 6 orbital pole vectors.
    pub r_pole_obs: f64,
    /// P(R̄_ϖ ≥ observed) under uniform ϖ.
    pub p_varpi: f64,
    /// P(R_pole ≥ observed) under uniform Ω (inclinations held at the
    /// observed values, as in the paper's null).
    pub p_pole: f64,
    /// Joint false-alarm probability: P(both exceedances simultaneously).
    pub p_joint: f64,
    pub n_trials: usize,
    /// RNG seed used (recorded for reproducibility).
    pub seed: u64,
}

/// Monte Carlo estimate of the joint significance of the ϖ and orbital-pole
/// clustering of the 6 stable KBOs.
///
/// Null hypothesis: each object's ϖ and Ω are independent and uniform on
/// the circle, with the observed inclinations held fixed (so the pole null
/// randomizes the node, not the tilt). Statistics: the circular mean
/// resultant length R̄ of the 6 ϖ's, and |Σ p̂|/6 over the 6 pole unit
/// vectors.
pub fn joint_clustering_significance(n_trials: usize, seed: u64) -> JointClusteringResult {
    use rand::Rng;
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let kbos = stable_kbos();
    let n = kbos.len();
    let varpis_obs: Vec<f64> = kbos
        .iter()
        .map(|k| longitude_of_perihelion(&k.elements))
        .collect();
    let poles_obs: Vec<[f64; 3]> = kbos.iter().map(|k| pole_vector(&k.elements)).collect();
    let inclinations: Vec<f64> = kbos.iter().map(|k| k.elements.i).collect();

    let r_bar_varpi_obs = mean_resultant_length(&varpis_obs);
    let r_pole_obs = pole_resultant_length(&poles_obs);

    let mut n_varpi = 0usize;
    let mut n_pole = 0usize;
    let mut n_joint = 0usize;

    let mut sim_varpis = vec![0.0; n];
    let mut sim_poles = vec![[0.0; 3]; n];

    for _ in 0..n_trials {
        for k in 0..n {
            sim_varpis[k] = rng.gen_range(0.0..TWO_PI);
            let node = rng.gen_range(0.0..TWO_PI);
            let (si, ci) = (inclinations[k].sin(), inclinations[k].cos());
            sim_poles[k] = [si * node.sin(), -si * node.cos(), ci];
        }

        let varpi_hit = mean_resultant_length(&sim_varpis) >= r_bar_varpi_obs;
        let pole_hit = pole_resultant_length(&sim_poles) >= r_pole_obs;
        if varpi_hit {
            n_varpi += 1;
        }
        if pole_hit {
            n_pole += 1;
        }
        if varpi_hit && pole_hit {
            n_joint += 1;
        }
    }

    let nt = n_trials as f64;
    JointClusteringResult {
        r_bar_varpi_obs,
        r_pole_obs,
        p_varpi: n_varpi as f64 / nt,
        p_pole: n_pole as f64 / nt,
        p_joint: n_joint as f64 / nt,
        n_trials,
        seed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_observed_stats_are_circular() {
        // Circular means must land in the observed clusters: ϖ near 71°,
        // Ω near 113° (paper values), where an arithmetic mean of wrapped
        // angles can be arbitrarily far off.
        let (varpi_mean, node_mean, _omega_mean) = observed_clustering_stats();
        let varpi_deg = varpi_mean * RAD2DEG;
        let node_deg = node_mean * RAD2DEG;
        assert!(
            (varpi_deg - 71.0).abs() < 20.0,
            "mean ϖ = {varpi_deg:.1}° (paper: 71 ± 16°)"
        );
        assert!(
            (node_deg - 113.0).abs() < 20.0,
            "mean Ω = {node_deg:.1}° (paper: 113 ± 13°)"
        );
    }

    #[test]
    fn test_clone_generation() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let kbos = stable_kbos();
        let clones = generate_clones(&kbos[0], 6, &mut rng);

        assert_eq!(clones.len(), 6);
        for clone in &clones {
            assert!(
                (clone.a - kbos[0].elements.a).abs()
                    < 10.0 * CLONE_SIGMA_A_FRAC * kbos[0].elements.a
            );
            assert!(clone.e > 0.0 && clone.e < 1.0);
        }
    }

    #[test]
    fn test_joint_clustering_smoke() {
        // Fast seeded MC: the joint probability is small and bounded above
        // by each marginal.
        let result = joint_clustering_significance(50_000, 7);
        assert!(
            result.r_bar_varpi_obs > 0.7,
            "R̄_ϖ = {}",
            result.r_bar_varpi_obs
        );
        assert!(result.r_pole_obs > 0.9, "R_pole = {}", result.r_pole_obs);
        assert!(result.p_varpi < 0.05, "p_ϖ = {}", result.p_varpi);
        assert!(result.p_pole < 0.10, "p_pole = {}", result.p_pole);
        assert!(result.p_joint <= result.p_varpi.min(result.p_pole) + 1e-12);
        assert!(result.p_joint < 5e-3, "p_joint = {}", result.p_joint);
        assert_eq!(result.seed, 7);
    }

    #[test]
    #[ignore = "high-N Monte Carlo (~10^7 trials); validates the paper's ~0.007% headline"]
    fn test_joint_clustering_matches_paper_order_of_magnitude() {
        let result = joint_clustering_significance(10_000_000, 7);
        // Paper: P ≈ 0.007% = 7e-5. Assert the same order of magnitude.
        assert!(
            result.p_joint > 5e-6 && result.p_joint < 1e-3,
            "p_joint = {:.2e} (paper ~7e-5)",
            result.p_joint
        );
    }

    #[test]
    fn test_integration_time_typed_matches_f64() {
        use approx::assert_relative_eq;
        let kbos = stable_kbos();
        let result = clone_stability_screen(&kbos[0], 2, 1e5 * YEAR_DAYS, 3000.0, 11);
        assert_relative_eq!(
            (result.integration_time() / days(1.0)).value,
            result.t_total,
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_clone_stability_screen_reduced() {
        // Reduced-scale screen: 1e5 yr, Neptune direct + J2-averaged JSU.
        // The sample objects all have q > 35 AU and must screen stable.
        let kbos = stable_kbos();
        let result = clone_stability_screen(&kbos[0], 3, 1e5 * YEAR_DAYS, 3000.0, 11);
        assert_eq!(result.n_clones, 3);
        assert_eq!(
            result.n_stable, 3,
            "Sedna clones must be stable over 1e5 yr (max Δa/a = {:.3})",
            result.max_da_frac
        );
        assert_eq!(result.seed, 11);
    }

    #[test]
    #[ignore = "full 4-Gyr clone-stability screen (hours)"]
    fn test_clone_stability_screen_full() {
        for kbo in &stable_kbos() {
            let result = clone_stability_screen(kbo, 8, 4.0 * GYR_DAYS, 3000.0, 11);
            assert!(
                result.n_stable * 2 > result.n_clones,
                "{}: only {}/{} clones stable",
                kbo.name,
                result.n_stable,
                result.n_clones
            );
        }
    }
}

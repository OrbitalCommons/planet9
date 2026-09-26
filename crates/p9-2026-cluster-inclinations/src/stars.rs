//! Passing field stars in the impulse approximation (paper Section 2.2;
//! Heisler & Tremaine 1986).
//!
//! A star of mass `M_*` passing at speed `v_*` on a straight line gives an
//! object at perpendicular offset `b` from that line the velocity kick
//! `Δv = 2 G M_* / (v_* b)` toward the line. The Sun receives its own kick,
//! so the heliocentric change of a particle is
//!
//! ```text
//!   Δv = (2 G M_* / v_*) · ( b_p / |b_p|² − b_s / |b_s|² ) ,
//! ```
//!
//! with `b_s` the impact vector from the Sun to the star's path and
//! `b_p = b_s − r_⊥` the same for the particle (`r_⊥` its heliocentric
//! position projected perpendicular to the path). Encounters are a Poisson
//! process with rate `Γ = n_* π b_max² v_*`.

use nalgebra::Vector3;
use p9_core::constants::{GM_SUN, KMS_TO_AUDAY, PC_AU, YEAR_DAYS};
use rand::Rng;

/// Passing-star population model.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct PassingStarModel {
    /// Stellar number density (pc⁻³).
    pub density_per_pc3: f64,
    /// Encounter speed (km/s). The paper sets every passing star to
    /// √2 × 1 km/s, the RMS speed of isotropic motion in the birth cluster.
    pub speed_kms: f64,
    /// Largest impact parameter simulated (pc).
    pub b_max_pc: f64,
    /// Perturber mass (solar masses).
    pub mass_solar: f64,
}

impl PassingStarModel {
    /// The paper's prescription: cluster-like √2 km/s encounters at the
    /// solar-neighbourhood density of 0.1 pc⁻³, mean-IMF 0.5 M☉ stars, out
    /// to 1 pc.
    pub fn paper() -> Self {
        Self {
            density_per_pc3: 0.1,
            speed_kms: std::f64::consts::SQRT_2,
            b_max_pc: 1.0,
            mass_solar: 0.5,
        }
    }

    /// Encounter rate Γ (per day) inside `b_max`.
    pub fn rate_per_day(&self) -> f64 {
        let n_au3 = self.density_per_pc3 / PC_AU.powi(3);
        let b_max = self.b_max_pc * PC_AU;
        n_au3 * std::f64::consts::PI * b_max * b_max * self.speed_kms * KMS_TO_AUDAY
    }

    /// Encounters expected per Myr.
    pub fn rate_per_myr(&self) -> f64 {
        self.rate_per_day() * 1.0e6 * YEAR_DAYS
    }

    /// Draw one encounter: isotropic direction of motion, impact parameter
    /// with density ∝ b on (0, b_max], impact vector uniform in azimuth
    /// around the path.
    pub fn sample_encounter<R: Rng>(&self, rng: &mut R) -> Encounter {
        let v_hat = random_unit(rng);
        let b = self.b_max_pc * PC_AU * rng.gen::<f64>().sqrt();
        // Perpendicular basis.
        let helper = if v_hat.x.abs() < 0.9 {
            Vector3::x()
        } else {
            Vector3::y()
        };
        let e1 = v_hat.cross(&helper).normalize();
        let e2 = v_hat.cross(&e1);
        let phi = rng.gen_range(0.0..std::f64::consts::TAU);
        Encounter {
            b_sun: b * (phi.cos() * e1 + phi.sin() * e2),
            v_hat,
            gm_over_v: 2.0 * GM_SUN * self.mass_solar / (self.speed_kms * KMS_TO_AUDAY),
        }
    }

    /// Poisson schedule of encounter epochs (days) on [0, t_end).
    pub fn schedule<R: Rng>(&self, t_end_days: f64, rng: &mut R) -> Vec<(f64, Encounter)> {
        let rate = self.rate_per_day();
        let mut t = 0.0;
        let mut out = Vec::new();
        loop {
            t += -rng.gen::<f64>().max(1e-300).ln() / rate;
            if t >= t_end_days {
                break out;
            }
            out.push((t, self.sample_encounter(rng)));
        }
    }
}

/// One straight-line stellar passage.
#[derive(Debug, Clone, Copy)]
pub struct Encounter {
    /// Impact vector from the Sun to the star's closest approach (AU).
    pub b_sun: Vector3<f64>,
    /// Unit direction of the star's motion.
    pub v_hat: Vector3<f64>,
    /// 2 G M_* / v_* (AU²/day).
    pub gm_over_v: f64,
}

impl Encounter {
    /// Heliocentric velocity kick (AU/day) on a particle at heliocentric
    /// position `pos` (AU).
    pub fn impulse(&self, pos: &Vector3<f64>) -> Vector3<f64> {
        let r_perp = pos - pos.dot(&self.v_hat) * self.v_hat;
        let b_p = self.b_sun - r_perp;
        let kick = |b: Vector3<f64>| b / b.norm_squared().max(1e-6);
        self.gm_over_v * (kick(b_p) - kick(self.b_sun))
    }
}

fn random_unit<R: Rng>(rng: &mut R) -> Vector3<f64> {
    let z: f64 = rng.gen_range(-1.0..1.0);
    let phi = rng.gen_range(0.0..std::f64::consts::TAU);
    let s = (1.0 - z * z).sqrt();
    Vector3::new(s * phi.cos(), s * phi.sin(), z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn rate_is_a_fraction_per_myr_at_field_density() {
        // n π b² v with n = 0.1 pc⁻³, b = 1 pc, v = 1.41 km/s ≈ 1.45 pc/Myr.
        let r = PassingStarModel::paper().rate_per_myr();
        assert!((0.4..0.5).contains(&r), "rate = {r} / Myr");
    }

    #[test]
    fn impulse_is_tidal_and_grows_with_distance() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let enc = PassingStarModel::paper().sample_encounter(&mut rng);
        let near = enc.impulse(&Vector3::new(100.0, 0.0, 0.0)).norm();
        let far = enc.impulse(&Vector3::new(10_000.0, 0.0, 0.0)).norm();
        assert!(far > near);
        // A particle at the Sun feels no heliocentric kick.
        assert!(enc.impulse(&Vector3::zeros()).norm() < 1e-30);
    }

    #[test]
    fn impulse_matches_tidal_limit() {
        // For r ≪ b the kick is the tidal expansion, |Δv| ~ (2GM/v) r / b².
        let enc = Encounter {
            b_sun: Vector3::new(1.0e5, 0.0, 0.0),
            v_hat: Vector3::z(),
            gm_over_v: 1.0,
        };
        let r = 100.0;
        let dv = enc.impulse(&Vector3::new(r, 0.0, 0.0)).norm();
        let tidal = r / 1.0e10;
        assert!((dv - tidal).abs() / tidal < 5e-3, "{dv} vs {tidal}");
    }
}

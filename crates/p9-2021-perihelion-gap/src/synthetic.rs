//! Seeded synthetic perihelion populations.
//!
//! Oldroyd & Trujillo (2021) argue the gap is unlikely to come from a single
//! continuous distribution and instead reflects two separate populations. We
//! make that concrete: draw (a) a seeded two-component population — a low-q
//! scattering component plus a high-q detached component with a sparse bridge
//! between them — which reproduces a perihelion gap; and (b) a single
//! continuous control population over the same q range, which does not. The
//! dip statistic (`crate::distribution`) is then computed on both.

use rand::SeedableRng;
use rand_distr::{Distribution, Normal, Uniform};

use p9_core::units::{au, Length};

/// Parameters of the seeded two-component (scattering + detached) population.
#[derive(Debug, Clone, Copy)]
pub struct TwoPopulationParams {
    /// Number of low-q scattering objects.
    pub n_scattering: usize,
    /// Mean perihelion of the scattering component (AU).
    pub scattering_q_mean: f64,
    /// Spread (std dev) of the scattering component (AU).
    pub scattering_q_sigma: f64,
    /// Number of high-q detached objects.
    pub n_detached: usize,
    /// Mean perihelion of the detached component (AU).
    pub detached_q_mean: f64,
    /// Spread (std dev) of the detached component (AU).
    pub detached_q_sigma: f64,
    /// Number of "bridge" objects spread thinly across the gap — the small
    /// real population the paper predicts (∼20% of the IOC density).
    pub n_bridge: usize,
    /// Bridge window (AU): [lo, hi).
    pub bridge_lo: f64,
    pub bridge_hi: f64,
}

impl Default for TwoPopulationParams {
    /// Defaults that place the scattering peak below the gap and the detached
    /// peak above it, matching the q ≈ 50–65 AU gap.
    fn default() -> Self {
        Self {
            n_scattering: 120,
            scattering_q_mean: 43.0,
            scattering_q_sigma: 4.0,
            n_detached: 60,
            detached_q_mean: 78.0,
            detached_q_sigma: 8.0,
            n_bridge: 6,
            bridge_lo: 50.0,
            bridge_hi: 65.0,
        }
    }
}

impl TwoPopulationParams {
    // ---- typed boundary accessors (dimension-checked views of the f64 fields) ----

    /// Mean perihelion of the scattering component as a [`Length`].
    pub fn scattering_mean_perihelion(&self) -> Length {
        au(self.scattering_q_mean)
    }
    /// Perihelion spread of the scattering component as a [`Length`].
    pub fn scattering_perihelion_sigma(&self) -> Length {
        au(self.scattering_q_sigma)
    }
    /// Mean perihelion of the detached component as a [`Length`].
    pub fn detached_mean_perihelion(&self) -> Length {
        au(self.detached_q_mean)
    }
    /// Perihelion spread of the detached component as a [`Length`].
    pub fn detached_perihelion_sigma(&self) -> Length {
        au(self.detached_q_sigma)
    }
    /// Lower edge of the bridge perihelion window as a [`Length`].
    pub fn bridge_lower(&self) -> Length {
        au(self.bridge_lo)
    }
    /// Upper edge of the bridge perihelion window as a [`Length`].
    pub fn bridge_upper(&self) -> Length {
        au(self.bridge_hi)
    }
}

/// Draw a seeded two-component perihelion population (AU). Perihelia are
/// clamped to be positive; the two Gaussian components flank the gap and the
/// bridge populates it sparsely.
pub fn two_population_perihelia(seed: u64, p: &TwoPopulationParams) -> Vec<f64> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let scat = Normal::new(p.scattering_q_mean, p.scattering_q_sigma).unwrap();
    let det = Normal::new(p.detached_q_mean, p.detached_q_sigma).unwrap();
    let bridge = Uniform::new(p.bridge_lo, p.bridge_hi);

    let mut out = Vec::with_capacity(p.n_scattering + p.n_detached + p.n_bridge);
    for _ in 0..p.n_scattering {
        out.push(scat.sample(&mut rng).max(30.0));
    }
    for _ in 0..p.n_detached {
        out.push(det.sample(&mut rng).max(30.0));
    }
    for _ in 0..p.n_bridge {
        out.push(bridge.sample(&mut rng));
    }
    out
}

/// Draw a seeded single continuous control population (AU): one broad uniform
/// distribution over [`lo`, `hi`) with the same total count. This is the
/// "single continuous distribution" null the paper rules out — it has no gap.
pub fn single_population_perihelia(seed: u64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let u = Uniform::new(lo, hi);
    (0..n).map(|_| u.sample(&mut rng)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::{dip_statistic, locate_gap_center, Histogram};
    use crate::published::{GAP_Q_CENTER_AU, GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
    use approx::assert_relative_eq;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_param_accessors_match_f64() {
        let p = TwoPopulationParams::default();
        assert_relative_eq!(
            p.scattering_mean_perihelion().get::<astronomical_unit>(),
            p.scattering_q_mean,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.scattering_perihelion_sigma().get::<astronomical_unit>(),
            p.scattering_q_sigma,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.detached_mean_perihelion().get::<astronomical_unit>(),
            p.detached_q_mean,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.detached_perihelion_sigma().get::<astronomical_unit>(),
            p.detached_q_sigma,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.bridge_lower().get::<astronomical_unit>(),
            p.bridge_lo,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.bridge_upper().get::<astronomical_unit>(),
            p.bridge_hi,
            epsilon = 1e-9
        );
    }

    #[test]
    fn two_population_is_reproducible() {
        let p = TwoPopulationParams::default();
        let a = two_population_perihelia(7, &p);
        let b = two_population_perihelia(7, &p);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x, y);
        }
    }

    #[test]
    fn two_population_shows_a_gap() {
        let p = TwoPopulationParams::default();
        let q = two_population_perihelia(2021, &p);
        let d = dip_statistic(&q, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, 12.0);
        // The seeded two-component draw has a clear deficit in the gap.
        assert!(
            d.dip_ratio < 0.6,
            "expected a deficit; dip_ratio = {}",
            d.dip_ratio
        );
    }

    #[test]
    fn single_continuous_population_has_no_gap() {
        let q = single_population_perihelia(2021, 200, 35.0, 95.0);
        let d = dip_statistic(&q, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, 12.0);
        // A flat distribution should have a dip ratio near 1 (no gap).
        assert!(
            (d.dip_ratio - 1.0).abs() < 0.4,
            "uniform control should not show a deficit; dip_ratio = {}",
            d.dip_ratio
        );
    }

    #[test]
    fn located_gap_matches_published_center() {
        let p = TwoPopulationParams::default();
        let q = two_population_perihelia(99, &p);
        // 13 bins of 5 AU over [35, 100): bin centres land on ...,52.5,57.5,62.5,...
        let h = Histogram::new(&q, 35.0, 100.0, 13);
        let center = locate_gap_center(&h).expect("gap present");
        // The located minimum-density bin sits within the published 50-65 AU
        // window, near the centre.
        assert!(
            center >= GAP_Q_LOW_AU && center <= GAP_Q_HIGH_AU,
            "gap center {center} outside published window"
        );
        assert!(
            (center - GAP_Q_CENTER_AU).abs() < 8.0,
            "gap center {center} far from published {GAP_Q_CENTER_AU}"
        );
    }

    #[test]
    fn bridge_density_is_well_below_detached_density() {
        // The paper predicts gap objects are only ~20% as numerous as the
        // nearby IOC population; our bridge is built sparse to match.
        let p = TwoPopulationParams::default();
        let q = two_population_perihelia(5, &p);
        let d = dip_statistic(&q, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, 15.0);
        assert!(d.gap_density < 0.5 * d.high_flank_density);
    }
}

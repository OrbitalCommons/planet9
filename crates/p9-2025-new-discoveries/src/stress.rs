//! How the two new objects degrade the ETNO apsidal-clustering signal.
//!
//! For the longitude of perihelion ϖ = ω + Ω we compute, reusing
//! `p9_core::analysis::circular`:
//!   - the mean resultant length R̄ (clustering strength), and
//!   - the small-n-corrected Rayleigh p-value (significance vs isotropy)
//! for the baseline `BROWN_2017_SAMPLE` and for that sample with each new
//! object (and both) appended. We also report each object's angular offset of
//! ϖ from the baseline circular mean.
//!
//! The Planet Nine "stress" is that adding an object the planet should have
//! aligned *lowers* R̄ and *raises* the Rayleigh p — i.e. it pushes the
//! clustering toward non-significance rather than reinforcing it.

use p9_core::analysis::circular::{circular_mean, mean_resultant_length, rayleigh_p_value};
use p9_core::constants::{RAD2DEG, TWO_PI};
use p9_core::data::etno::{Etno, BROWN_2017_SAMPLE};

use crate::discoveries::NewDiscovery;

/// Clustering statistics for one sample of longitudes of perihelion.
#[derive(Debug, Clone, Copy)]
pub struct ClusteringStat {
    /// Number of objects in the sample.
    pub n: usize,
    /// Mean resultant length R̄ ∈ [0, 1].
    pub r_bar: f64,
    /// Small-n-corrected Rayleigh p-value (vs circular uniformity).
    pub rayleigh_p: f64,
    /// Circular mean of ϖ (degrees, [0, 360)).
    pub mean_varpi_deg: f64,
}

/// Compute R̄, the Rayleigh p-value and the circular mean ϖ for a slice of
/// `Etno` entries.
pub fn clustering_stat(sample: &[Etno]) -> ClusteringStat {
    let varpis: Vec<f64> = sample.iter().map(|e| e.longitude_of_perihelion()).collect();
    ClusteringStat {
        n: sample.len(),
        r_bar: mean_resultant_length(&varpis),
        rayleigh_p: rayleigh_p_value(&varpis),
        mean_varpi_deg: circular_mean(&varpis)
            .map(|m| m * RAD2DEG)
            .unwrap_or(f64::NAN),
    }
}

/// Angular offset (degrees, in [0, 180]) of an object's ϖ from a reference
/// longitude of perihelion (radians).
pub fn varpi_offset_deg(object: &Etno, reference_varpi: f64) -> f64 {
    let mut d = (object.longitude_of_perihelion() - reference_varpi).rem_euclid(TWO_PI);
    if d > TWO_PI / 2.0 {
        d = TWO_PI - d;
    }
    d * RAD2DEG
}

/// Full stress analysis: baseline, +OF201, +Ammonite, +both, plus each
/// object's ϖ offset from the baseline circular mean.
#[derive(Debug, Clone, Copy)]
pub struct StressResult {
    /// Baseline `BROWN_2017_SAMPLE` clustering.
    pub baseline: ClusteringStat,
    /// Baseline + 2017 OF201.
    pub with_of201: ClusteringStat,
    /// Baseline + Ammonite.
    pub with_ammonite: ClusteringStat,
    /// Baseline + both new objects.
    pub with_both: ClusteringStat,
    /// 2017 OF201 ϖ offset from the baseline circular mean (deg).
    pub of201_offset_deg: f64,
    /// Ammonite ϖ offset from the baseline circular mean (deg).
    pub ammonite_offset_deg: f64,
}

/// Run the stress analysis on the vetted baseline sample and the two new
/// objects from [`crate::discoveries`].
pub fn run_stress(of201: &NewDiscovery, ammonite: &NewDiscovery) -> StressResult {
    let baseline_vec: Vec<Etno> = BROWN_2017_SAMPLE.to_vec();

    let mut with_of201 = baseline_vec.clone();
    with_of201.push(of201.etno);

    let mut with_ammonite = baseline_vec.clone();
    with_ammonite.push(ammonite.etno);

    let mut with_both = baseline_vec.clone();
    with_both.push(of201.etno);
    with_both.push(ammonite.etno);

    let baseline = clustering_stat(&baseline_vec);
    let baseline_mean_varpi = baseline.mean_varpi_deg / RAD2DEG;

    StressResult {
        baseline,
        with_of201: clustering_stat(&with_of201),
        with_ammonite: clustering_stat(&with_ammonite),
        with_both: clustering_stat(&with_both),
        of201_offset_deg: varpi_offset_deg(&of201.etno, baseline_mean_varpi),
        ammonite_offset_deg: varpi_offset_deg(&ammonite.etno, baseline_mean_varpi),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::discoveries::{ammonite, of201};

    fn result() -> StressResult {
        run_stress(&of201(), &ammonite())
    }

    #[test]
    fn baseline_is_clustered() {
        let r = result();
        // The defining baseline feature: clustered ϖ, significant Rayleigh.
        assert!(r.baseline.r_bar > 0.5, "R̄ = {:.3}", r.baseline.r_bar);
        assert!(
            r.baseline.rayleigh_p < 0.05,
            "p = {:.4}",
            r.baseline.rayleigh_p
        );
        assert_eq!(r.baseline.n, 10);
    }

    #[test]
    fn both_new_objects_far_from_cluster() {
        let r = result();
        // Each new object's ϖ offset from the baseline circular mean
        // exceeds ~60° — they are NOT apsidally aligned with the cluster.
        assert!(
            r.of201_offset_deg > 60.0,
            "2017 OF201 offset = {:.1}°",
            r.of201_offset_deg
        );
        assert!(
            r.ammonite_offset_deg > 60.0,
            "Ammonite offset = {:.1}°",
            r.ammonite_offset_deg
        );
    }

    #[test]
    fn adding_objects_reduces_clustering() {
        let r = result();
        // Each addition lowers R̄ relative to baseline.
        assert!(r.with_of201.r_bar < r.baseline.r_bar);
        assert!(r.with_ammonite.r_bar < r.baseline.r_bar);
        assert!(r.with_both.r_bar < r.baseline.r_bar);
        // And the strongest degradation is from adding both.
        assert!(r.with_both.r_bar < r.with_of201.r_bar);
        assert!(r.with_both.r_bar < r.with_ammonite.r_bar);
    }

    #[test]
    fn adding_objects_raises_rayleigh_p() {
        let r = result();
        // Adding either object weakens the significance; adding both pushes
        // the clustering across the conventional p = 0.05 threshold into
        // non-significance.
        assert!(r.with_of201.rayleigh_p > r.baseline.rayleigh_p);
        assert!(r.with_ammonite.rayleigh_p > r.baseline.rayleigh_p);
        assert!(r.with_both.rayleigh_p > r.baseline.rayleigh_p);
        assert!(
            r.with_both.rayleigh_p > 0.05,
            "p(baseline+both) = {:.4} should lose significance",
            r.with_both.rayleigh_p
        );
    }

    #[test]
    fn residuals_documented() {
        // Honest record of the computed numbers (JPL elements, epoch
        // 2461200.5). These are the values the assertions above key on:
        //   baseline      : R̄ ≈ 0.649, p ≈ 0.011, mean ϖ ≈ 52°
        //   +2017 OF201   : R̄ ≈ 0.572, p ≈ 0.024  (offset ≈ 106°)
        //   +Ammonite     : R̄ ≈ 0.522, p ≈ 0.047  (offset ≈ 141°)
        //   +both         : R̄ ≈ 0.472, p ≈ 0.067  (no longer significant)
        let r = result();
        assert!((r.baseline.r_bar - 0.649).abs() < 0.02, "{:?}", r.baseline);
        assert!((r.baseline.mean_varpi_deg - 52.3).abs() < 2.0);
        assert!(
            (r.with_both.r_bar - 0.472).abs() < 0.02,
            "{:?}",
            r.with_both
        );
        assert!((r.with_both.rayleigh_p - 0.067).abs() < 0.02);
        assert!((r.of201_offset_deg - 105.7).abs() < 2.0);
        assert!((r.ammonite_offset_deg - 141.5).abs() < 2.0);
    }
}

//! The serialized data contract between the Rust survey computation and the
//! Python plotter.
//!
//! This is the *only* interface the plotter depends on: the `bin/survey`
//! binary writes a [`SurveyDataset`] as JSON, and `scripts/plot_survey.py`
//! reads it back through a mirror of these structs (typed dataclasses with
//! field-by-field validation). Keep the two sides in lock-step — bump
//! [`SCHEMA_VERSION`] on any breaking change so the Python loader rejects a
//! stale file instead of silently mis-plotting.
//!
//! Units are encoded in field names (`_deg`, `_au`, `_deg2`, `_mag`). All
//! angles are degrees, distances AU, areas square degrees, magnitudes V (or
//! the telescope's stated band, treated as a reflected-light proxy).

use p9_core::analysis::surveys::Footprint;
use p9_core::units::{au, degrees, earth_masses, Angle, Length, Mass};
use serde::{Deserialize, Serialize};

/// Schema version. The Python loader asserts the file matches this exactly.
pub const SCHEMA_VERSION: u32 = 2;

/// One paper's Planet Nine orbit solution and the (documented, assumed)
/// 1σ spreads used to turn it into a sky-position probability cloud.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitSolution {
    /// Short label, e.g. "2021 Brown & Batygin (MCMC)".
    pub name: String,
    /// Bibliographic citation.
    pub citation: String,
    /// arXiv identifier (empty if none).
    pub arxiv: String,
    /// Perturber mass (Earth masses).
    pub mass_earth: f64,
    /// Assumed geometric albedo (V band) for the brightness estimate.
    pub albedo: f64,
    /// Semi-major axis (AU) and its assumed 1σ.
    pub a_au: f64,
    pub a_sigma_au: f64,
    /// Eccentricity and its assumed 1σ.
    pub e: f64,
    pub e_sigma: f64,
    /// Inclination (deg) and its assumed 1σ.
    pub i_deg: f64,
    pub i_sigma_deg: f64,
    /// Argument of perihelion (deg) and its assumed 1σ.
    pub omega_deg: f64,
    pub omega_sigma_deg: f64,
    /// Longitude of ascending node (deg) and its assumed 1σ.
    pub omega_big_deg: f64,
    pub omega_big_sigma_deg: f64,
    /// Provenance note: where each number comes from, what is assumed.
    pub note: String,
}

impl OrbitSolution {
    // ---- typed boundary accessors (dimension-checked views of the wire
    // fields; the serialized schema is unchanged) ----

    /// Perturber mass as a dimension-checked [`Mass`] (Earth-mass storage).
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }
    /// Semi-major axis as a dimension-checked [`Length`] (AU storage).
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }
    /// Inclination as a dimension-checked [`Angle`] (degree storage).
    pub fn inclination(&self) -> Angle {
        degrees(self.i_deg)
    }
    /// Argument of perihelion as a dimension-checked [`Angle`] (degree storage).
    pub fn argument_of_perihelion(&self) -> Angle {
        degrees(self.omega_deg)
    }
    /// Longitude of ascending node as a dimension-checked [`Angle`] (degree
    /// storage).
    pub fn longitude_of_node(&self) -> Angle {
        degrees(self.omega_big_deg)
    }
}

/// The fixed RA/Dec grid the probability maps live on. Cell `(i_ra, i_dec)`
/// is stored at flat index `i_dec * n_ra + i_ra` (row-major, Dec outer).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkyGrid {
    pub ra_min_deg: f64,
    pub ra_max_deg: f64,
    pub n_ra: usize,
    pub dec_min_deg: f64,
    pub dec_max_deg: f64,
    pub n_dec: usize,
}

impl SkyGrid {
    pub fn dra(&self) -> f64 {
        (self.ra_max_deg - self.ra_min_deg) / self.n_ra as f64
    }
    pub fn ddec(&self) -> f64 {
        (self.dec_max_deg - self.dec_min_deg) / self.n_dec as f64
    }
    pub fn len(&self) -> usize {
        self.n_ra * self.n_dec
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Flat index for an (RA, Dec) sample, or `None` when the sample falls
    /// outside the grid. (A previous version clamped out-of-band samples
    /// into the edge rows — with the ±60° Dec grid and the sampled i/Ω
    /// dispersions, high-declination draws DO occur, and the clamp piled a
    /// few percent of spurious probability onto the ±58.5° rows, distorting
    /// area68/95 and the edge-row targets.)
    pub fn index(&self, ra_deg: f64, dec_deg: f64) -> Option<usize> {
        let fx = (ra_deg - self.ra_min_deg) / self.dra();
        let fy = (dec_deg - self.dec_min_deg) / self.ddec();
        if fx < 0.0 || fy < 0.0 {
            return None;
        }
        let (ix, iy) = (fx as usize, fy as usize);
        if ix >= self.n_ra || iy >= self.n_dec {
            return None;
        }
        Some(iy * self.n_ra + ix)
    }
    /// Dec at the center of row `i_dec` (for cos-weighted solid angle).
    pub fn dec_center(&self, i_dec: usize) -> f64 {
        self.dec_min_deg + (i_dec as f64 + 0.5) * self.ddec()
    }
}

/// Computed sky/brightness summary for one study.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StudyResult {
    pub solution: OrbitSolution,
    /// Position probability per grid cell, normalized to sum = 1
    /// (length `grid.len()`, same indexing as [`SkyGrid::index`]).
    pub prob: Vec<f64>,
    /// Most-probable cell center.
    pub peak_ra_deg: f64,
    pub peak_dec_deg: f64,
    /// Smallest cos(Dec)-weighted area enclosing 68% / 95% of the probability.
    pub area68_deg2: f64,
    pub area95_deg2: f64,
    /// Current heliocentric distance distribution (dwell-weighted), AU.
    pub dist_p16_au: f64,
    pub dist_median_au: f64,
    pub dist_p84_au: f64,
    /// Apparent V (reflected light) distribution.
    pub v_p16: f64,
    pub v_median: f64,
    pub v_p84: f64,
}

impl StudyResult {
    /// Most-probable-cell RA as a dimension-checked [`Angle`] (degree storage).
    pub fn peak_ra(&self) -> Angle {
        degrees(self.peak_ra_deg)
    }
    /// Most-probable-cell Dec as a dimension-checked [`Angle`] (degree storage).
    pub fn peak_dec(&self) -> Angle {
        degrees(self.peak_dec_deg)
    }
    /// 16th-percentile heliocentric distance as a [`Length`] (AU storage).
    pub fn dist_p16(&self) -> Length {
        au(self.dist_p16_au)
    }
    /// Median heliocentric distance as a [`Length`] (AU storage).
    pub fn dist_median(&self) -> Length {
        au(self.dist_median_au)
    }
    /// 84th-percentile heliocentric distance as a [`Length`] (AU storage).
    pub fn dist_p84(&self) -> Length {
        au(self.dist_p84_au)
    }
}

/// A telescope/survey: a limiting magnitude in some band plus a footprint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Telescope {
    pub name: String,
    pub band: String,
    pub limiting_mag: f64,
    pub footprint: Footprint,
    pub space_based: bool,
    pub note: String,
}

/// One telescope's detection odds against one study.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StudyDetection {
    pub study_name: String,
    /// P(in footprint AND brighter than depth), including coverage fraction.
    pub detection_prob: f64,
    /// P(in footprint), including coverage fraction (geometry only).
    pub prob_in_footprint: f64,
    /// P(brighter than the limiting magnitude), ignoring footprint.
    pub prob_bright_enough: f64,
}

/// A telescope plus its detection odds against every study.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TelescopeResult {
    pub telescope: Telescope,
    pub per_study: Vec<StudyDetection>,
}

/// Sky overlays (precomputed in Rust from p9-core frame matrices so the
/// plotter does no astronomy): polylines of [RA_deg, Dec_deg].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Overlays {
    pub galactic_plane: Vec<[f64; 2]>,
    pub ecliptic: Vec<[f64; 2]>,
}

/// Search-narrowing constraints for the plotter to overlay: the Rubin cede
/// line and the Cassini/Iorio ephemeris favored-ν sky arc.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Constraints {
    /// Rubin/LSST northern limit (deg); JBT cedes everything at/below it.
    pub rubin_dec_max_deg: f64,
    /// Favored true-anomaly interval (deg) from the Cassini ephemeris fit.
    pub favored_nu_lo_deg: f64,
    pub favored_nu_hi_deg: f64,
    /// Favored-ν zone mapped onto the sky (RA, Dec deg) for the 2021 orbit.
    pub favored_arc: Vec<[f64; 2]>,
}

impl Constraints {
    /// Rubin/LSST northern declination limit as a dimension-checked [`Angle`]
    /// (degree storage).
    pub fn rubin_dec_max(&self) -> Angle {
        degrees(self.rubin_dec_max_deg)
    }
    /// Lower edge of the favored true-anomaly interval as an [`Angle`].
    pub fn favored_nu_lo(&self) -> Angle {
        degrees(self.favored_nu_lo_deg)
    }
    /// Upper edge of the favored true-anomaly interval as an [`Angle`].
    pub fn favored_nu_hi(&self) -> Angle {
        degrees(self.favored_nu_hi_deg)
    }
}

/// The complete dataset written to JSON and read by the plotter.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurveyDataset {
    pub schema_version: u32,
    pub generated_by: String,
    /// Monte Carlo samples per study (provenance for reproducibility).
    pub samples_per_study: usize,
    pub rng_seed: u64,
    pub grid: SkyGrid,
    pub studies: Vec<StudyResult>,
    pub telescopes: Vec<TelescopeResult>,
    pub overlays: Overlays,
    pub constraints: Constraints,
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::degree;
    use uom::si::length::astronomical_unit;
    use uom::si::mass::kilogram;

    fn sample_solution() -> OrbitSolution {
        OrbitSolution {
            name: "test".into(),
            citation: String::new(),
            arxiv: String::new(),
            mass_earth: 6.0,
            albedo: 0.3,
            a_au: 460.0,
            a_sigma_au: 50.0,
            e: 0.2,
            e_sigma: 0.05,
            i_deg: 16.0,
            i_sigma_deg: 5.0,
            omega_deg: 150.0,
            omega_sigma_deg: 10.0,
            omega_big_deg: 100.0,
            omega_big_sigma_deg: 10.0,
            note: String::new(),
        }
    }

    #[test]
    fn orbit_solution_typed_accessors_match_fields() {
        let s = sample_solution();
        assert_relative_eq!(
            s.mass().get::<kilogram>(),
            earth_masses(s.mass_earth).get::<kilogram>(),
            max_relative = 1e-12
        );
        assert_relative_eq!(
            s.semi_major_axis().get::<astronomical_unit>(),
            s.a_au,
            epsilon = 1e-12
        );
        assert_relative_eq!(s.inclination().get::<degree>(), s.i_deg, epsilon = 1e-12);
        assert_relative_eq!(
            s.argument_of_perihelion().get::<degree>(),
            s.omega_deg,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            s.longitude_of_node().get::<degree>(),
            s.omega_big_deg,
            epsilon = 1e-12
        );
    }

    #[test]
    fn study_result_typed_accessors_match_fields() {
        let r = StudyResult {
            solution: sample_solution(),
            prob: vec![1.0],
            peak_ra_deg: 85.0,
            peak_dec_deg: 20.0,
            area68_deg2: 100.0,
            area95_deg2: 300.0,
            dist_p16_au: 400.0,
            dist_median_au: 500.0,
            dist_p84_au: 650.0,
            v_p16: 19.0,
            v_median: 21.0,
            v_p84: 23.0,
        };
        assert_relative_eq!(r.peak_ra().get::<degree>(), r.peak_ra_deg, epsilon = 1e-12);
        assert_relative_eq!(
            r.peak_dec().get::<degree>(),
            r.peak_dec_deg,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            r.dist_median().get::<astronomical_unit>(),
            r.dist_median_au,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            r.dist_p16().get::<astronomical_unit>(),
            r.dist_p16_au,
            epsilon = 1e-12
        );
    }

    #[test]
    fn constraints_typed_accessors_match_fields() {
        let c = Constraints {
            rubin_dec_max_deg: 12.0,
            favored_nu_lo_deg: 30.0,
            favored_nu_hi_deg: 60.0,
            favored_arc: vec![],
        };
        assert_relative_eq!(
            c.rubin_dec_max().get::<degree>(),
            c.rubin_dec_max_deg,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            c.favored_nu_lo().get::<degree>(),
            c.favored_nu_lo_deg,
            epsilon = 1e-12
        );
    }

    #[test]
    fn out_of_grid_samples_are_dropped_not_clamped() {
        let grid = crate::default_grid();
        assert!(grid.index(100.0, 0.0).is_some());
        assert!(
            grid.index(100.0, 65.0).is_none(),
            "above +60 must not clamp"
        );
        assert!(grid.index(100.0, -61.0).is_none());
        assert!(grid.index(-1.0, 0.0).is_none());
    }
}

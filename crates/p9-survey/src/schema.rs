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
    /// Flat index for an (RA, Dec) sample, clamped into the grid.
    pub fn index(&self, ra_deg: f64, dec_deg: f64) -> usize {
        let ix = (((ra_deg - self.ra_min_deg) / self.dra()) as isize)
            .clamp(0, self.n_ra as isize - 1) as usize;
        let iy = (((dec_deg - self.dec_min_deg) / self.ddec()) as isize)
            .clamp(0, self.n_dec as isize - 1) as usize;
        iy * self.n_ra + ix
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

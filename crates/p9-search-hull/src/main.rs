//! Where to point next: cross the Planet Nine current-position posterior against
//! the combined coverage hull of every wide optical search, and report the
//! residual — the high-probability sky that has NOT yet been searched deeply
//! enough — restricted to what a small all-sky space telescope (the JBT 0.5 m /
//! SPENCER concept, V ≈ 24.5) could actually detect.
//!
//! Two products, both written into one JSON for `scripts/plot_search_hull.py`:
//!
//! 1. **Sky priority map.** Per direction: the ensemble posterior probability,
//!    the predicted distance and apparent V, the deepest survey that has
//!    reached it (the coverage hull), the probability it would already have been
//!    caught, and the residual un-searched probability. Plus a ranked list of
//!    space-telescope pointings.
//! 2. **Parameter-space reach hull.** The (distance × apparent V) plane with the
//!    posterior cloud and the survey-depth cuts that bound what has been seen.
//!
//! Everything is computed from the workspace's own models (p9-core photometry +
//! coordinates, p9-survey orbit solutions), not hand-drawn.

mod hull;
mod posterior;
mod surveys;

use p9_core::constants::DEG2RAD;
use serde::Serialize;

use p9_survey::{default_grid, schema::SkyGrid, studies, telescope};

/// MC draws per orbit solution.
const PER_STUDY: usize = 60_000;
/// RNG seed (reproducible).
const SEED: u64 = 20_260_626;
/// Posterior (distance, V) cloud points to keep for the parameter-space panel.
const CLOUD_KEEP: usize = 4_000;
/// How many ranked pointings to emit.
const N_TARGETS: usize = 40;

#[derive(Serialize)]
struct GridSpec {
    ra_min_deg: f64,
    ra_max_deg: f64,
    n_ra: usize,
    dec_min_deg: f64,
    dec_max_deg: f64,
    n_dec: usize,
    ra_centers_deg: Vec<f64>,
    dec_centers_deg: Vec<f64>,
}

#[derive(Serialize)]
struct SurveyOut {
    name: &'static str,
    depth: f64,
    dec_min_deg: f64,
    dec_max_deg: f64,
    galactic_lat_min_deg: f64,
    coverage_fraction: f64,
    reference: &'static str,
}

#[derive(Serialize)]
struct TelescopeOut {
    name: String,
    limiting_mag: f64,
    dec_min_deg: f64,
    dec_max_deg: f64,
    galactic_lat_min_deg: f64,
    coverage_fraction: f64,
}

/// Flat per-cell arrays, row-major (idx = i_dec * n_ra + i_ra).
#[derive(Serialize)]
struct Cells {
    gal_b_deg: Vec<f64>,
    prob: Vec<f64>,
    dist_mean_au: Vec<Option<f64>>,
    v_mean: Vec<Option<f64>>,
    best_depth: Vec<Option<f64>>,
    best_survey: Vec<Option<&'static str>>,
    exclusion_prob: Vec<f64>,
    residual_prob: Vec<f64>,
    reach_dist_au: Vec<Option<f64>>,
    space_reachable: Vec<bool>,
}

#[derive(Serialize)]
struct Target {
    rank: usize,
    ra_deg: f64,
    dec_deg: f64,
    gal_b_deg: f64,
    prob: f64,
    residual_prob: f64,
    v_mean: f64,
    dist_mean_au: f64,
    best_depth: Option<f64>,
    best_survey: Option<&'static str>,
    /// How many magnitudes brighter than the space-telescope limit P9 sits here.
    space_margin_mag: f64,
}

#[derive(Serialize)]
struct CloudPoint {
    dist_au: f64,
    v_mag: f64,
}

#[derive(Serialize)]
struct AllSkyDepth {
    name: &'static str,
    depth: f64,
    reach_au: Option<f64>,
}

#[derive(Serialize)]
struct Hull {
    fiducial_mass_earth: f64,
    distance_au: Vec<f64>,
    v_curve: Vec<f64>,
    allsky_depths: Vec<AllSkyDepth>,
    space_depth: f64,
    space_reach_au: Option<f64>,
    cloud: Vec<CloudPoint>,
}

#[derive(Serialize)]
struct Summary {
    /// Posterior probability already excluded by the optical coverage hull.
    excluded_prob: f64,
    /// Posterior probability that remains un-searched.
    residual_prob: f64,
    /// Residual probability our space telescope could actually detect.
    residual_space_reachable_prob: f64,
    /// Sky area (deg²) searched deeper than the given V thresholds.
    area_searched_deg2: Vec<(f64, f64)>,
    sky_area_in_grid_deg2: f64,
}

#[derive(Serialize)]
struct Dataset {
    generated_by: &'static str,
    samples_per_study: usize,
    rng_seed: u64,
    grid: GridSpec,
    space_telescope: TelescopeOut,
    surveys: Vec<SurveyOut>,
    cells: Cells,
    targets: Vec<Target>,
    hull: Hull,
    summary: Summary,
}

fn ra_center(grid: &SkyGrid, i_ra: usize) -> f64 {
    grid.ra_min_deg + (i_ra as f64 + 0.5) * grid.dra()
}

fn main() {
    let grid = default_grid();
    let studies = studies::catalog();
    let surveys = surveys::real_surveys();

    // The space telescope we are pointing.
    let scope = telescope::catalog()
        .into_iter()
        .find(|t| t.space_based)
        .expect("a space-based telescope in the catalog");
    let scope_depth = scope.limiting_mag;
    let scope_fp = scope.footprint.clone();

    // Ensemble current-position posterior with per-cell predicted distance / V.
    let post = posterior::ensemble(&studies, &grid, PER_STUDY, SEED, CLOUD_KEEP);

    let n = grid.len();
    let mut gal_b = vec![0.0_f64; n];
    let mut dist_mean = vec![None; n];
    let mut v_mean = vec![None; n];
    let mut best_depth = vec![None; n];
    let mut best_survey: Vec<Option<&'static str>> = vec![None; n];
    let mut excl = vec![0.0_f64; n];
    let mut residual = vec![0.0_f64; n];
    let mut reach_dist = vec![None; n];
    let mut space_reachable = vec![false; n];

    let cell_base = grid.dra() * grid.ddec(); // deg² before cos(Dec) weighting
    let mut sky_area = 0.0_f64;
    let thresholds = [20.5_f64, 21.5, 23.0];
    let mut area_deeper = [0.0_f64; 3];

    let mut excluded_prob = 0.0;
    let mut residual_prob = 0.0;
    let mut residual_reach = 0.0;

    for iy in 0..grid.n_dec {
        let dec = grid.dec_center(iy);
        let solid = cell_base * (dec * DEG2RAD).cos();
        for ix in 0..grid.n_ra {
            let idx = iy * grid.n_ra + ix;
            sky_area += solid;

            // Galactic latitude of this direction (via an ecliptic-frame proxy
            // position on the unit sphere at the cell center).
            let ra = ra_center(&grid, ix);
            gal_b[idx] = galactic_lat_at(ra, dec);

            // Coverage hull at this direction.
            let cov = surveys::cell_coverage(&surveys, dec, gal_b[idx]);
            best_depth[idx] = cov.best_depth;
            best_survey[idx] = cov.best_survey;
            if let Some(d) = cov.best_depth {
                reach_dist[idx] = hull::reach_distance(hull::FIDUCIAL_MASS_EARTH, d);
                if d >= thresholds[0] {
                    area_deeper[0] += solid;
                }
                if d >= thresholds[1] {
                    area_deeper[1] += solid;
                }
                if d >= thresholds[2] {
                    area_deeper[2] += solid;
                }
            }

            let p = post.prob[idx];
            if p <= 0.0 || !post.v_mean[idx].is_finite() {
                continue;
            }
            let v = post.v_mean[idx];
            let r = post.dist_mean[idx];
            v_mean[idx] = Some(v);
            dist_mean[idx] = Some(r);

            let pe = surveys::exclusion_probability(&surveys, dec, gal_b[idx], v);
            excl[idx] = pe;
            residual[idx] = p * (1.0 - pe);

            excluded_prob += p * pe;
            residual_prob += residual[idx];

            let reachable = scope_fp.accepts(dec, gal_b[idx]) && v <= scope_depth;
            space_reachable[idx] = reachable;
            if reachable {
                residual_reach += residual[idx];
            }
        }
    }

    // Rank space-reachable cells by residual probability.
    let mut ranked: Vec<usize> = (0..n)
        .filter(|&k| space_reachable[k] && residual[k] > 0.0)
        .collect();
    ranked.sort_by(|&a, &b| residual[b].partial_cmp(&residual[a]).unwrap());
    let targets: Vec<Target> = ranked
        .iter()
        .take(N_TARGETS)
        .enumerate()
        .map(|(rank, &k)| {
            let iy = k / grid.n_ra;
            let ix = k % grid.n_ra;
            Target {
                rank: rank + 1,
                ra_deg: ra_center(&grid, ix),
                dec_deg: grid.dec_center(iy),
                gal_b_deg: gal_b[k],
                prob: post.prob[k],
                residual_prob: residual[k],
                v_mean: v_mean[k].unwrap(),
                dist_mean_au: dist_mean[k].unwrap(),
                best_depth: best_depth[k],
                best_survey: best_survey[k],
                space_margin_mag: scope_depth - v_mean[k].unwrap(),
            }
        })
        .collect();

    // Parameter-space hull.
    let (distance_au, v_curve) = hull::magnitude_curve(120);
    let allsky_depths = vec![
        allsky("CRTS", &surveys),
        allsky("ZTF", &surveys),
        allsky("PS1 3pi", &surveys),
    ];
    let hull = Hull {
        fiducial_mass_earth: hull::FIDUCIAL_MASS_EARTH,
        distance_au,
        v_curve,
        allsky_depths,
        space_depth: scope_depth,
        space_reach_au: hull::reach_distance(hull::FIDUCIAL_MASS_EARTH, scope_depth),
        cloud: post
            .cloud
            .iter()
            .map(|s| CloudPoint {
                dist_au: s.dist_au,
                v_mag: s.v_mag,
            })
            .collect(),
    };

    let dataset = Dataset {
        generated_by: "p9-search-hull",
        samples_per_study: PER_STUDY,
        rng_seed: SEED,
        grid: GridSpec {
            ra_min_deg: grid.ra_min_deg,
            ra_max_deg: grid.ra_max_deg,
            n_ra: grid.n_ra,
            dec_min_deg: grid.dec_min_deg,
            dec_max_deg: grid.dec_max_deg,
            n_dec: grid.n_dec,
            ra_centers_deg: (0..grid.n_ra).map(|i| ra_center(&grid, i)).collect(),
            dec_centers_deg: (0..grid.n_dec).map(|i| grid.dec_center(i)).collect(),
        },
        space_telescope: TelescopeOut {
            name: scope.name.clone(),
            limiting_mag: scope_depth,
            dec_min_deg: scope_fp.dec_min_deg,
            dec_max_deg: scope_fp.dec_max_deg,
            galactic_lat_min_deg: scope_fp.galactic_lat_min_deg,
            coverage_fraction: scope_fp.coverage_fraction,
        },
        surveys: surveys
            .iter()
            .map(|s| SurveyOut {
                name: s.name,
                depth: s.depth,
                dec_min_deg: s.footprint.dec_min_deg,
                dec_max_deg: s.footprint.dec_max_deg,
                galactic_lat_min_deg: s.footprint.galactic_lat_min_deg,
                coverage_fraction: s.footprint.coverage_fraction,
                reference: s.reference,
            })
            .collect(),
        cells: Cells {
            gal_b_deg: gal_b,
            prob: post.prob.clone(),
            dist_mean_au: dist_mean,
            v_mean,
            best_depth,
            best_survey,
            exclusion_prob: excl,
            residual_prob: residual,
            reach_dist_au: reach_dist,
            space_reachable,
        },
        targets,
        hull,
        summary: Summary {
            excluded_prob,
            residual_prob,
            residual_space_reachable_prob: residual_reach,
            area_searched_deg2: thresholds
                .iter()
                .zip(area_deeper.iter())
                .map(|(&t, &a)| (t, a))
                .collect(),
            sky_area_in_grid_deg2: sky_area,
        },
    };

    let json = serde_json::to_string_pretty(&dataset).expect("serialize");
    std::fs::write("figures/search_hull.json", json).expect("write figures/search_hull.json");
    eprintln!(
        "wrote figures/search_hull.json  (excluded {:.0}%, residual {:.0}%, space-reachable residual {:.0}%)",
        100.0 * excluded_prob,
        100.0 * residual_prob,
        100.0 * residual_reach,
    );
}

/// Galactic latitude (deg) of an equatorial direction (ra, dec), via a unit
/// ecliptic-frame vector — reuses the posterior module's frame chain.
fn galactic_lat_at(ra_deg: f64, dec_deg: f64) -> f64 {
    use nalgebra::Vector3;
    use p9_core::coords::sky::equatorial_to_ecliptic;
    let (lon, lat) = equatorial_to_ecliptic(ra_deg.to_radians(), dec_deg.to_radians());
    let v = Vector3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
    posterior::galactic_latitude_deg(&v)
}

fn allsky(name: &'static str, surveys: &[surveys::RealSurvey]) -> AllSkyDepth {
    let depth = surveys.iter().find(|s| s.name == name).unwrap().depth;
    AllSkyDepth {
        name,
        depth,
        reach_au: hull::reach_distance(hull::FIDUCIAL_MASS_EARTH, depth),
    }
}

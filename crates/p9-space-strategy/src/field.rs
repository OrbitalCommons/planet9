//! The value of every tile: how much probability is there, how much of it
//! nobody else will collect, and what a given exposure would capture.

use serde::Serialize;

use crate::crowding::completeness;
use crate::instrument::SpaceTelescope;
use crate::prior::Draw;
use crate::tiles::Tiles;

/// Stars this much fainter than the detection limit still cost area.
pub const CROWDING_MARGIN_MAG: f64 = 0.5;

/// How much of Rubin's eventual harvest the plan concedes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum RubinPolicy {
    /// Plan as if Rubin did not exist.
    Ignore,
    /// Concede the wide-fast-deep footprint (Dec ≤ +12°).
    Wfd,
    /// Concede wide-fast-deep and the North Ecliptic Spur.
    WfdPlusNes,
}

impl RubinPolicy {
    pub fn label(self) -> &'static str {
        match self {
            RubinPolicy::Ignore => "Race Rubin everywhere",
            RubinPolicy::Wfd => "Concede Rubin's main survey",
            RubinPolicy::WfdPlusNes => "Concede all of Rubin's footprint",
        }
    }

    fn p_rubin(self, d: &Draw) -> f64 {
        match self {
            RubinPolicy::Ignore => 0.0,
            RubinPolicy::Wfd => d.p_rubin_wfd,
            RubinPolicy::WfdPlusNes => d.p_rubin_nes,
        }
    }
}

/// What the plan is asked to maximise.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Policy {
    pub rubin: RubinPolicy,
    /// Weight draws by the Cassini-ranging orbital-phase prior.
    pub cassini_phase: bool,
}

impl Policy {
    pub const BASELINE: Policy = Policy {
        rubin: RubinPolicy::WfdPlusNes,
        cassini_phase: false,
    };
}

/// One tile's share of the prior.
#[derive(Debug, Clone, Serialize)]
pub struct Tile {
    pub index: usize,
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub ra_width_deg: f64,
    pub dec_height_deg: f64,
    pub area_deg2: f64,
    pub gal_l_deg: f64,
    pub gal_b_deg: f64,
    pub ecl_lon_deg: f64,
    pub ecl_lat_deg: f64,
    /// Prior probability that Planet Nine is in this tile.
    pub prior: f64,
    /// ... and has not already been found from the ground.
    pub residual: f64,
    /// ... and will not be found by Rubin either (per the policy).
    pub unique: f64,
    /// Fraction of the tile the space telescope keeps after crowding at its
    /// reference depth (tiers are scored at their own depth).
    pub crowd: f64,
    /// Unique weight and magnitude of every draw in the tile.
    #[serde(skip)]
    pub draws: Vec<(f64, f64)>,
}

impl Tile {
    /// Fraction of the tile that survives crowding in an exposure reaching
    /// `depth`: stars within [`CROWDING_MARGIN_MAG`] below the detection
    /// limit still leave residuals.
    pub fn crowding(&self, scope: &SpaceTelescope, depth: f64) -> f64 {
        completeness(
            depth + CROWDING_MARGIN_MAG,
            scope.mask_radius_arcsec,
            self.gal_l_deg,
            self.gal_b_deg,
        )
    }

    /// Probability captured by surveying this tile to `depth` per visit.
    pub fn captured(&self, scope: &SpaceTelescope, depth: f64) -> f64 {
        self.crowding(scope, depth)
            * self
                .draws
                .iter()
                .map(|&(w, v)| w * scope.link_probability(v, depth))
                .sum::<f64>()
    }

    /// Probability density of the unique residual (per deg²).
    pub fn unique_density(&self) -> f64 {
        self.unique / self.area_deg2
    }
}

/// The scored sky for one prior and policy.
#[derive(Debug, Clone, Serialize)]
pub struct Sky {
    pub policy: Policy,
    pub tiles: Vec<Tile>,
    pub prior_total: f64,
    pub found_by_ground: f64,
    pub residual_total: f64,
    pub conceded_to_rubin: f64,
    pub unique_total: f64,
}

/// Accumulate scored draws onto the tiling.
pub fn score_sky(draws: &[Draw], tiling: &Tiles, scope: &SpaceTelescope, policy: Policy) -> Sky {
    let weight_of = |d: &Draw| {
        if policy.cassini_phase {
            d.nu_weight
        } else {
            1.0
        }
    };
    let norm: f64 = draws.iter().map(weight_of).sum();

    let mut prior = vec![0.0; tiling.len()];
    let mut residual = vec![0.0; tiling.len()];
    let mut unique = vec![0.0; tiling.len()];
    let mut members: Vec<Vec<(f64, f64)>> = vec![Vec::new(); tiling.len()];
    for d in draws {
        let w = weight_of(d) / norm;
        let k = tiling.index(d.ra_deg, d.dec_deg);
        let r = w * (1.0 - d.p_ground);
        let u = r * (1.0 - policy.rubin.p_rubin(d));
        prior[k] += w;
        residual[k] += r;
        unique[k] += u;
        if u > 0.0 {
            members[k].push((u, d.v_mag));
        }
    }

    let tiles: Vec<Tile> = (0..tiling.len())
        .map(|k| {
            use p9_core::coords::sky::{equatorial_to_ecliptic_deg, equatorial_to_galactic};
            let (ra, dec) = tiling.centre(k);
            let (w, h) = tiling.extent(k);
            let (l, b) = equatorial_to_galactic(ra.to_radians(), dec.to_radians());
            let (lon, lat) = equatorial_to_ecliptic_deg(ra, dec);
            Tile {
                index: k,
                ra_deg: ra,
                dec_deg: dec,
                ra_width_deg: w,
                dec_height_deg: h,
                area_deg2: tiling.area_deg2(k),
                gal_l_deg: l.to_degrees(),
                gal_b_deg: b.to_degrees(),
                ecl_lon_deg: lon.rem_euclid(360.0),
                ecl_lat_deg: lat,
                prior: prior[k],
                residual: residual[k],
                unique: unique[k],
                crowd: completeness(
                    scope.reference_depth + CROWDING_MARGIN_MAG,
                    scope.mask_radius_arcsec,
                    l.to_degrees(),
                    b.to_degrees(),
                ),
                draws: std::mem::take(&mut members[k]),
            }
        })
        .collect();

    let prior_total: f64 = tiles.iter().map(|t| t.prior).sum();
    let residual_total: f64 = tiles.iter().map(|t| t.residual).sum();
    let unique_total: f64 = tiles.iter().map(|t| t.unique).sum();
    Sky {
        policy,
        tiles,
        prior_total,
        found_by_ground: prior_total - residual_total,
        residual_total,
        conceded_to_rubin: residual_total - unique_total,
        unique_total,
    }
}

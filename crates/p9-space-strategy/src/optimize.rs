//! Spend a telescope-time budget where it buys the most probability.
//!
//! Every tile can be left alone or surveyed at one of the exposure tiers;
//! each choice has a cost (wall-clock hours) and a gain (probability
//! captured). This is a multiple-choice knapsack. Its greedy solution on the
//! upper concave hull of each tile's (cost, gain) options is optimal at every
//! budget that falls on a hull vertex and within one increment of optimal
//! elsewhere: take increments in order of gain per hour until the budget is
//! spent.

use serde::Serialize;

use crate::field::{Sky, Tile};
use crate::instrument::SpaceTelescope;

/// Per-visit integrations considered (seconds).
pub const TIERS_S: [f64; 7] = [60.0, 120.0, 300.0, 600.0, 1200.0, 2400.0, 4800.0];

/// One step up a tile's hull.
#[derive(Debug, Clone, Copy)]
struct Increment {
    tile: usize,
    tier: usize,
    hours: f64,
    gain: f64,
}

impl Increment {
    fn rate(&self) -> f64 {
        self.gain / self.hours
    }
}

/// Probability of finding Planet Nine in `tile` at per-visit `depth`.
pub fn planet_nine_gain(tile: &Tile, scope: &SpaceTelescope, depth: f64) -> f64 {
    tile.captured(scope, depth)
}

fn options(
    tile: &Tile,
    scope: &SpaceTelescope,
    gain: &impl Fn(&Tile, &SpaceTelescope, f64) -> f64,
) -> Vec<(f64, f64)> {
    TIERS_S
        .iter()
        .map(|&t| {
            (
                scope.wall_hours(tile.area_deg2, t),
                gain(tile, scope, scope.depth(t, tile.ecl_lat_deg)),
            )
        })
        .collect()
}

fn hull_increments(tile_pos: usize, opts: &[(f64, f64)]) -> Vec<Increment> {
    let mut out = Vec::new();
    let (mut cost, mut gain) = (0.0, 0.0);
    let mut from = 0;
    loop {
        let best = (from..opts.len())
            .filter(|&k| opts[k].1 > gain + 1e-15 && opts[k].0 > cost)
            .max_by(|&a, &b| {
                let ra = (opts[a].1 - gain) / (opts[a].0 - cost);
                let rb = (opts[b].1 - gain) / (opts[b].0 - cost);
                ra.partial_cmp(&rb).unwrap()
            });
        let Some(k) = best else { break };
        out.push(Increment {
            tile: tile_pos,
            tier: k,
            hours: opts[k].0 - cost,
            gain: opts[k].1 - gain,
        });
        (cost, gain) = opts[k];
        from = k + 1;
    }
    out
}

/// A tile in the plan.
#[derive(Debug, Clone, Serialize)]
pub struct Allocation {
    /// Position in `Sky::tiles`.
    pub tile: usize,
    /// Index into [`TIERS_S`].
    pub tier: usize,
    pub integration_s: f64,
    pub depth: f64,
    pub hours: f64,
    pub captured: f64,
}

/// A plan for one budget.
#[derive(Debug, Clone, Serialize)]
pub struct Plan {
    pub budget_hours: f64,
    pub hours: f64,
    pub captured: f64,
    pub area_deg2: f64,
    pub n_fields: f64,
    pub allocations: Vec<Allocation>,
}

/// The efficient frontier and the means to cut a plan from it.
#[derive(Debug, Clone)]
pub struct Frontier {
    increments: Vec<Increment>,
    /// Cumulative (hours, captured) after each increment.
    pub curve: Vec<(f64, f64)>,
}

impl Frontier {
    /// The frontier for finding Planet Nine.
    pub fn build(sky: &Sky, scope: &SpaceTelescope) -> Self {
        Self::build_for(sky, scope, planet_nine_gain)
    }

    /// The frontier for any per-tile objective `gain(tile, scope, depth)`.
    pub fn build_for(
        sky: &Sky,
        scope: &SpaceTelescope,
        gain: impl Fn(&Tile, &SpaceTelescope, f64) -> f64,
    ) -> Self {
        let mut increments: Vec<Increment> = sky
            .tiles
            .iter()
            .enumerate()
            .flat_map(|(pos, t)| hull_increments(pos, &options(t, scope, &gain)))
            .collect();
        increments.sort_by(|a, b| b.rate().partial_cmp(&a.rate()).unwrap());
        let mut curve = Vec::with_capacity(increments.len() + 1);
        let (mut h, mut g) = (0.0, 0.0);
        curve.push((h, g));
        for inc in &increments {
            h += inc.hours;
            g += inc.gain;
            curve.push((h, g));
        }
        Self { increments, curve }
    }

    /// Probability captured at `hours` (interpolated along the frontier).
    pub fn captured_at(&self, hours: f64) -> f64 {
        let k = self.curve.partition_point(|p| p.0 <= hours);
        if k == 0 {
            return 0.0;
        }
        if k == self.curve.len() {
            return self.curve[k - 1].1;
        }
        let (a, b) = (self.curve[k - 1], self.curve[k]);
        a.1 + (b.1 - a.1) * (hours - a.0) / (b.0 - a.0)
    }

    /// Hours at which the marginal return falls to `rate` (probability per
    /// hour): the natural stopping point for a given patience.
    pub fn hours_at_rate(&self, rate: f64) -> f64 {
        let k = self.increments.partition_point(|i| i.rate() >= rate);
        self.curve[k].0
    }

    /// The plan that spends at most `budget_hours` on finding Planet Nine.
    pub fn plan(&self, sky: &Sky, scope: &SpaceTelescope, budget_hours: f64) -> Plan {
        self.plan_for(sky, scope, budget_hours, planet_nine_gain)
    }

    /// The plan that spends at most `budget_hours`, with `captured` reported
    /// in units of the objective `gain`.
    pub fn plan_for(
        &self,
        sky: &Sky,
        scope: &SpaceTelescope,
        budget_hours: f64,
        gain: impl Fn(&Tile, &SpaceTelescope, f64) -> f64,
    ) -> Plan {
        let mut tier_of: std::collections::BTreeMap<usize, usize> = Default::default();
        let mut spent = 0.0;
        for inc in &self.increments {
            if spent + inc.hours > budget_hours {
                break;
            }
            spent += inc.hours;
            tier_of.insert(inc.tile, inc.tier);
        }
        let allocations: Vec<Allocation> = tier_of
            .into_iter()
            .map(|(pos, tier)| {
                let t = &sky.tiles[pos];
                let depth = scope.depth(TIERS_S[tier], t.ecl_lat_deg);
                Allocation {
                    tile: pos,
                    tier,
                    integration_s: TIERS_S[tier],
                    depth,
                    hours: scope.wall_hours(t.area_deg2, TIERS_S[tier]),
                    captured: gain(t, scope, depth),
                }
            })
            .collect();
        let area: f64 = allocations
            .iter()
            .map(|a| sky.tiles[a.tile].area_deg2)
            .sum();
        Plan {
            budget_hours,
            hours: allocations.iter().map(|a| a.hours).sum(),
            captured: allocations.iter().map(|a| a.captured).sum(),
            area_deg2: area,
            n_fields: scope.fields_for(area),
            allocations,
        }
    }
}

/// Probability `plan` (cut from `planned`) would capture if the truth were
/// the prior behind `actual`. Tiles are matched by their sky index.
pub fn evaluate(plan: &Plan, planned: &Sky, actual: &Sky, scope: &SpaceTelescope) -> f64 {
    let by_index: std::collections::HashMap<usize, &Tile> =
        actual.tiles.iter().map(|t| (t.index, t)).collect();
    plan.allocations
        .iter()
        .filter_map(|a| {
            by_index
                .get(&planned.tiles[a.tile].index)
                .map(|t| t.captured(scope, a.depth))
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hull_is_concave_and_skips_dominated_tiers() {
        // Tier 1 is dominated (worse rate than jumping straight to tier 2).
        let opts = [(1.0, 0.12), (2.0, 0.11), (3.0, 0.30), (6.0, 0.36)];
        let inc = hull_increments(0, &opts);
        let tiers: Vec<usize> = inc.iter().map(|i| i.tier).collect();
        assert_eq!(tiers, vec![0, 2, 3]);
        for w in inc.windows(2) {
            assert!(w[0].rate() >= w[1].rate());
        }
        let total: f64 = inc.iter().map(|i| i.gain).sum();
        assert!((total - 0.36).abs() < 1e-12);
    }
}

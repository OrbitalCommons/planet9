//! Assemble the strategy: score the sky, cut the plans, test them against
//! other priors, and write the dataset, the tables and the figures.

use std::collections::BTreeMap;
use std::fmt::Write;

use serde::Serialize;

use crate::crowding::completeness;
use crate::field::{score_sky, Policy, RubinPolicy, Sky, Tile};
use crate::instrument::{sky_rate_arcsec_per_hr, SpaceTelescope};
use crate::optimize::{evaluate, planet_nine_gain, Frontier, Plan, TIERS_S};
use crate::prior::{sample, Draw, PriorKind, RUBIN_SEEING};
use crate::season::{opposition_month, window_half_width_days, MONTHS};
use crate::svg::{self, Axis, Svg};
use crate::tiles::Tiles;
use crate::tno;
use crate::zones::{in_rubin_footprint, zone_of, Zone};

/// Draws per prior.
pub const N_DRAWS: usize = 600_000;
/// Seed (the date this strategy was first cut).
pub const SEED: u64 = 20_260_927;
/// Marginal return at which the reference campaign stops: one percent of
/// the prior per thousand hours.
pub const STOP_RATE_PER_HOUR: f64 = 1.0e-5;
/// Budgets tabulated alongside the reference campaign (wall-clock hours).
pub const BUDGETS_H: [f64; 6] = [250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0];

#[derive(Serialize)]
pub struct PolicyRow {
    pub policy: &'static str,
    pub cassini_phase: bool,
    pub found_by_ground: f64,
    pub residual: f64,
    pub conceded_to_rubin: f64,
    pub unique: f64,
    /// Probability captured at each of [`BUDGETS_H`].
    pub captured: Vec<f64>,
}

#[derive(Serialize)]
pub struct RobustnessRow {
    pub prior: &'static str,
    pub unique: f64,
    /// Captured by the reference plan (cut for Brown & Batygin 2021).
    pub reference_plan: f64,
    /// Captured by the plan cut for this prior at the same hours.
    pub own_plan: f64,
}

#[derive(Serialize)]
pub struct ZoneRow {
    pub zone: &'static str,
    pub why: &'static str,
    pub area_deg2: f64,
    pub fields: f64,
    pub hours: f64,
    pub captured: f64,
    pub median_integration_s: f64,
    pub median_depth: f64,
    pub ra_range_h: (f64, f64),
    pub dec_range_deg: (f64, f64),
    pub opposition_months: String,
    /// Median distance and brightness of the Planet Nine draws captured.
    pub median_dist_au: f64,
    pub median_v: f64,
    /// Hours needed between linked visits at that distance, at opposition.
    pub min_baseline_hr: f64,
    pub tnos: f64,
    pub tnos_new: f64,
}

#[derive(Serialize)]
pub struct TileRow {
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub ra_width_deg: f64,
    pub dec_height_deg: f64,
    pub gal_l_deg: f64,
    pub gal_b_deg: f64,
    pub ecl_lat_deg: f64,
    pub zone: &'static str,
    pub unique: f64,
    pub integration_s: f64,
    pub depth: f64,
    pub hours: f64,
    pub captured: f64,
    pub opposition_month: u32,
}

#[derive(Serialize)]
pub struct SkyCell {
    pub ra_deg: f64,
    pub dec_deg: f64,
    pub ra_width_deg: f64,
    pub dec_height_deg: f64,
    pub prior: f64,
    pub residual: f64,
    pub unique: f64,
}

#[derive(Serialize)]
pub struct Report {
    pub generated_by: &'static str,
    pub seed: u64,
    pub draws_per_prior: usize,
    pub telescope: SpaceTelescope,
    pub tiers_s: Vec<f64>,
    pub tier_depths_ecliptic: Vec<f64>,
    pub budgets_h: Vec<f64>,
    pub policies: Vec<PolicyRow>,
    pub reference_hours: f64,
    pub reference_captured: f64,
    pub reference_area_deg2: f64,
    pub reference_fields: f64,
    pub reference_share_of_unique: f64,
    pub reference_share_of_residual: f64,
    /// Sky Rubin also covers, taken to the same marginal return.
    pub extension_hours: f64,
    pub extension_captured: f64,
    pub extension_area_deg2: f64,
    pub zones: Vec<ZoneRow>,
    pub hours_by_month: Vec<f64>,
    pub window_half_width_days: f64,
    pub robustness: Vec<RobustnessRow>,
    pub frontier: Vec<(f64, f64)>,
    pub tiles: Vec<TileRow>,
    pub sky: Vec<SkyCell>,
    pub tno_plan_total: f64,
    pub tno_plan_new: f64,
    pub tno_plan_distant_new: f64,
    pub deep_program_hours: f64,
    pub deep_program_distant_new: f64,
    pub deep_program_tnos_new: f64,
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v[v.len() / 2]
}

fn rubin_crowding(t: &Tile) -> f64 {
    completeness(
        RUBIN_SEEING.star_depth,
        RUBIN_SEEING.mask_radius_arcsec,
        t.gal_l_deg,
        t.gal_b_deg,
    )
}

/// Expected TNOs in a tile at `depth`, and those Rubin will not also find.
fn tno_yield(t: &Tile, scope: &SpaceTelescope, depth: f64) -> (f64, f64) {
    let all = tno::expected_tnos(t.area_deg2, t.ecl_lat_deg, depth) * t.crowding(scope, depth);
    let share = tno::rubin_share(
        depth,
        in_rubin_footprint(t.dec_deg, t.ecl_lat_deg),
        rubin_crowding(t),
    );
    (all, all * (1.0 - share))
}

/// Objective of the reference campaign: Planet Nine, in the sky Rubin cannot
/// search. Tiles Rubin also covers are left to the optional extension.
fn core_gain(t: &Tile, scope: &SpaceTelescope, depth: f64) -> f64 {
    if zone_of(t) == Zone::RubinOverlap {
        0.0
    } else {
        planet_nine_gain(t, scope, depth)
    }
}

/// Objective of the deep programme: distant TNOs Rubin will not find.
fn distant_new_gain(t: &Tile, scope: &SpaceTelescope, depth: f64) -> f64 {
    tno::DISTANT_FRACTION * tno_yield(t, scope, depth).1
}

fn policy_row(draws: &[Draw], tiling: &Tiles, scope: &SpaceTelescope, policy: Policy) -> PolicyRow {
    let sky = score_sky(draws, tiling, scope, policy);
    let f = Frontier::build(&sky, scope);
    PolicyRow {
        policy: policy.rubin.label(),
        cassini_phase: policy.cassini_phase,
        found_by_ground: sky.found_by_ground,
        residual: sky.residual_total,
        conceded_to_rubin: sky.conceded_to_rubin,
        unique: sky.unique_total,
        captured: BUDGETS_H.iter().map(|&h| f.captured_at(h)).collect(),
    }
}

fn zone_rows(
    plan: &Plan,
    sky: &Sky,
    draws: &[Draw],
    tiling: &Tiles,
    scope: &SpaceTelescope,
) -> Vec<ZoneRow> {
    // Tier of every planned tile, by sky index, to find the captured draws.
    let depth_of: BTreeMap<usize, (Zone, f64)> = plan
        .allocations
        .iter()
        .map(|a| {
            let t = &sky.tiles[a.tile];
            (t.index, (zone_of(t), a.depth))
        })
        .collect();
    let mut captured_draws: BTreeMap<Zone, Vec<(f64, f64)>> = BTreeMap::new();
    for d in draws {
        if let Some(&(zone, depth)) = depth_of.get(&tiling.index(d.ra_deg, d.dec_deg)) {
            if d.v_mag <= depth && d.p_ground < 0.5 {
                captured_draws
                    .entry(zone)
                    .or_default()
                    .push((d.dist_au, d.v_mag));
            }
        }
    }

    Zone::ALL
        .iter()
        .filter_map(|&zone| {
            let members: Vec<_> = plan
                .allocations
                .iter()
                .filter(|a| zone_of(&sky.tiles[a.tile]) == zone)
                .collect();
            if members.is_empty() {
                return None;
            }
            let tiles: Vec<&Tile> = members.iter().map(|a| &sky.tiles[a.tile]).collect();
            let area: f64 = tiles.iter().map(|t| t.area_deg2).sum();
            // RA range as the smallest arc: unwrap about the mean direction.
            let (s, c) = tiles.iter().fold((0.0, 0.0), |(s, c), t| {
                let r = t.ra_deg.to_radians();
                (s + r.sin(), c + r.cos())
            });
            let mid = s.atan2(c).to_degrees();
            let off: Vec<f64> = tiles
                .iter()
                .map(|t| (t.ra_deg - mid + 540.0).rem_euclid(360.0) - 180.0)
                .collect();
            let lo = off.iter().cloned().fold(f64::MAX, f64::min);
            let hi = off.iter().cloned().fold(f64::MIN, f64::max);
            let mut months = [0.0; 12];
            for (a, t) in members.iter().zip(&tiles) {
                months[opposition_month(t.ecl_lon_deg) as usize - 1] += a.hours;
            }
            let total_h: f64 = months.iter().sum();
            let named: Vec<&str> = (0..12)
                .map(|k| (k + 6) % 12)
                .filter(|&k| months[k] > 0.08 * total_h)
                .map(|k| MONTHS[k])
                .collect();
            let drawn = captured_draws.remove(&zone).unwrap_or_default();
            let dist = median(drawn.iter().map(|d| d.0).collect());
            let (tnos, tnos_new) = members.iter().zip(&tiles).fold((0.0, 0.0), |acc, (a, t)| {
                let y = tno_yield(t, scope, a.depth);
                (acc.0 + y.0, acc.1 + y.1)
            });
            Some(ZoneRow {
                zone: zone.label(),
                why: zone.why(),
                area_deg2: area,
                fields: scope.fields_for(area),
                hours: total_h,
                captured: members.iter().map(|a| a.captured).sum(),
                median_integration_s: median(members.iter().map(|a| a.integration_s).collect()),
                median_depth: median(members.iter().map(|a| a.depth).collect()),
                ra_range_h: (
                    (mid + lo - 2.0).rem_euclid(360.0) / 15.0,
                    (mid + hi + 2.0).rem_euclid(360.0) / 15.0,
                ),
                dec_range_deg: (
                    tiles.iter().map(|t| t.dec_deg).fold(f64::MAX, f64::min) - 1.0,
                    tiles.iter().map(|t| t.dec_deg).fold(f64::MIN, f64::max) + 1.0,
                ),
                opposition_months: named.join(" "),
                median_dist_au: dist,
                median_v: median(drawn.iter().map(|d| d.1).collect()),
                min_baseline_hr: scope.min_motion_arcsec / sky_rate_arcsec_per_hr(dist, 0.0),
                tnos,
                tnos_new,
            })
        })
        .collect()
}

/// Run everything at the published sample size.
pub fn build() -> Report {
    build_with(N_DRAWS)
}

/// Run everything with `n_draws` draws per prior.
pub fn build_with(n_draws: usize) -> Report {
    let scope = SpaceTelescope::default();
    let tiling = Tiles::new();
    let draws = sample(PriorKind::Bb21, n_draws, SEED);

    let mut policies = Vec::new();
    for rubin in [
        RubinPolicy::Ignore,
        RubinPolicy::Wfd,
        RubinPolicy::WfdPlusNes,
    ] {
        policies.push(policy_row(
            &draws,
            &tiling,
            &scope,
            Policy {
                rubin,
                cassini_phase: false,
            },
        ));
    }
    policies.push(policy_row(
        &draws,
        &tiling,
        &scope,
        Policy {
            rubin: RubinPolicy::WfdPlusNes,
            cassini_phase: true,
        },
    ));

    let sky = score_sky(&draws, &tiling, &scope, Policy::BASELINE);
    let frontier = Frontier::build_for(&sky, &scope, core_gain);
    let reference_hours = (frontier.hours_at_rate(STOP_RATE_PER_HOUR) / 50.0).round() * 50.0;
    let plan = frontier.plan_for(&sky, &scope, reference_hours, core_gain);

    // The optional extension: sky Rubin also covers, taken to the same
    // marginal return.
    let open = Frontier::build(&sky, &scope);
    let open_plan = open.plan(&sky, &scope, open.hours_at_rate(STOP_RATE_PER_HOUR));
    let (extension_hours, extension_captured, extension_area) = open_plan
        .allocations
        .iter()
        .filter(|a| zone_of(&sky.tiles[a.tile]) == Zone::RubinOverlap)
        .fold((0.0, 0.0, 0.0), |acc, a| {
            (
                acc.0 + a.hours,
                acc.1 + a.captured,
                acc.2 + sky.tiles[a.tile].area_deg2,
            )
        });

    let robustness = PriorKind::ALL
        .iter()
        .map(|&kind| {
            let other = if kind == PriorKind::Bb21 {
                sky.clone()
            } else {
                score_sky(
                    &sample(kind, n_draws, SEED + 1),
                    &tiling,
                    &scope,
                    Policy::BASELINE,
                )
            };
            let own = Frontier::build_for(&other, &scope, core_gain).captured_at(plan.hours);
            RobustnessRow {
                prior: kind.label(),
                unique: other.unique_total,
                reference_plan: evaluate(&plan, &sky, &other, &scope),
                own_plan: own,
            }
        })
        .collect();

    let zones = zone_rows(&plan, &sky, &draws, &tiling, &scope);
    let mut hours_by_month = vec![0.0; 12];
    let tiles: Vec<TileRow> = plan
        .allocations
        .iter()
        .map(|a| {
            let t = &sky.tiles[a.tile];
            let month = opposition_month(t.ecl_lon_deg);
            hours_by_month[month as usize - 1] += a.hours;
            TileRow {
                ra_deg: t.ra_deg,
                dec_deg: t.dec_deg,
                ra_width_deg: t.ra_width_deg,
                dec_height_deg: t.dec_height_deg,
                gal_l_deg: t.gal_l_deg,
                gal_b_deg: t.gal_b_deg,
                ecl_lat_deg: t.ecl_lat_deg,
                zone: zone_of(t).label(),
                unique: t.unique,
                integration_s: a.integration_s,
                depth: a.depth,
                hours: a.hours,
                captured: a.captured,
                opposition_month: month,
            }
        })
        .collect();

    let mut tiles = tiles;
    tiles.sort_by(|a: &TileRow, b: &TileRow| {
        (b.captured / b.hours)
            .partial_cmp(&(a.captured / a.hours))
            .unwrap()
    });

    let (tno_total, tno_new) = plan.allocations.iter().fold((0.0, 0.0), |acc, a| {
        let y = tno_yield(&sky.tiles[a.tile], &scope, a.depth);
        (acc.0 + y.0, acc.1 + y.1)
    });

    // The deep programme: same telescope, same hours, spent on distant TNOs
    // Rubin will not find.
    let deep_frontier = Frontier::build_for(&sky, &scope, distant_new_gain);
    let deep_hours = (0.5 * reference_hours / 50.0).round() * 50.0;
    let deep = deep_frontier.plan_for(&sky, &scope, deep_hours, distant_new_gain);

    let step = (frontier.curve.len() / 400).max(1);
    Report {
        generated_by: "p9-space-strategy",
        seed: SEED,
        draws_per_prior: n_draws,
        tiers_s: TIERS_S.to_vec(),
        tier_depths_ecliptic: TIERS_S.iter().map(|&t| scope.depth(t, 0.0)).collect(),
        budgets_h: BUDGETS_H.to_vec(),
        policies,
        reference_hours: plan.hours,
        reference_captured: plan.captured,
        reference_area_deg2: plan.area_deg2,
        reference_fields: plan.n_fields,
        reference_share_of_unique: plan.captured / sky.unique_total,
        reference_share_of_residual: plan.captured / sky.residual_total,
        extension_hours,
        extension_captured,
        extension_area_deg2: extension_area,
        zones,
        hours_by_month,
        window_half_width_days: window_half_width_days(0.7),
        robustness,
        frontier: frontier
            .curve
            .iter()
            .step_by(step)
            .cloned()
            .filter(|p| p.0 <= 12_000.0)
            .collect(),
        tiles,
        sky: sky
            .tiles
            .iter()
            .filter(|t| t.prior > 0.0)
            .map(|t| SkyCell {
                ra_deg: t.ra_deg,
                dec_deg: t.dec_deg,
                ra_width_deg: t.ra_width_deg,
                dec_height_deg: t.dec_height_deg,
                prior: t.prior,
                residual: t.residual,
                unique: t.unique,
            })
            .collect(),
        tno_plan_total: tno_total,
        tno_plan_new: tno_new,
        tno_plan_distant_new: tno::DISTANT_FRACTION * tno_new,
        deep_program_hours: deep.hours,
        deep_program_distant_new: deep.captured,
        deep_program_tnos_new: deep.captured / tno::DISTANT_FRACTION,
        telescope: scope,
    }
}

// ---- tables -----------------------------------------------------------------

fn pct(x: f64) -> String {
    format!("{:.1}%", 100.0 * x)
}

pub fn tables(r: &Report) -> String {
    let mut s = String::new();
    writeln!(
        s,
        "<!-- Generated by `cargo run --release -p p9-space-strategy`. Do not edit. -->\n"
    )
    .unwrap();
    writeln!(s, "## Where the probability is\n").unwrap();
    writeln!(
        s,
        "Brown & Batygin (2021) prior, {} draws.\n",
        r.draws_per_prior
    )
    .unwrap();
    writeln!(s, "| Planning stance | Found from the ground | Left | Conceded to Rubin | Left for this telescope |").unwrap();
    writeln!(s, "|---|---|---|---|---|").unwrap();
    for p in &r.policies {
        let name = if p.cassini_phase {
            format!("{} + Cassini phase", p.policy)
        } else {
            p.policy.to_string()
        };
        writeln!(
            s,
            "| {} | {} | {} | {} | {} |",
            name,
            pct(p.found_by_ground),
            pct(p.residual),
            pct(p.conceded_to_rubin),
            pct(p.unique)
        )
        .unwrap();
    }

    writeln!(s, "\n## What a budget buys\n").unwrap();
    writeln!(
        s,
        "Probability of finding Planet Nine, as a share of the whole prior.\n"
    )
    .unwrap();
    let head: Vec<String> = r.budgets_h.iter().map(|h| format!("{h:.0} h")).collect();
    writeln!(s, "| Planning stance | {} |", head.join(" | ")).unwrap();
    writeln!(s, "|---|{}", "---|".repeat(r.budgets_h.len())).unwrap();
    for p in &r.policies {
        let name = if p.cassini_phase {
            format!("{} + Cassini phase", p.policy)
        } else {
            p.policy.to_string()
        };
        let cells: Vec<String> = p.captured.iter().map(|&c| pct(c)).collect();
        writeln!(s, "| {} | {} |", name, cells.join(" | ")).unwrap();
    }

    writeln!(s, "\n## The reference campaign\n").unwrap();
    writeln!(
        s,
        "{:.0} wall-clock hours, {:.0} deg², {:.0} fields. Captures {} of the prior: {} of what is left for this telescope, {} of everything not yet found.\n",
        r.reference_hours,
        r.reference_area_deg2,
        r.reference_fields,
        pct(r.reference_captured),
        pct(r.reference_share_of_unique),
        pct(r.reference_share_of_residual)
    )
    .unwrap();
    writeln!(
        s,
        "Optional extension into sky Rubin also covers, to the same marginal return: {:.0} h over {:.0} deg² for a further {}.\n",
        r.extension_hours,
        r.extension_area_deg2,
        pct(r.extension_captured)
    )
    .unwrap();
    writeln!(
        s,
        "| Zone | RA | Dec | Area (deg²) | Hours | Visit (s) | Depth (V) | Captured | Opposition |"
    )
    .unwrap();
    writeln!(s, "|---|---|---|---|---|---|---|---|---|").unwrap();
    for z in &r.zones {
        writeln!(
            s,
            "| {} | {:.1}h–{:.1}h | {:+.0}° to {:+.0}° | {:.0} | {:.0} | {:.0} | {:.1} | {} | {} |",
            z.zone,
            z.ra_range_h.0,
            z.ra_range_h.1,
            z.dec_range_deg.0,
            z.dec_range_deg.1,
            z.area_deg2,
            z.hours,
            z.median_integration_s,
            z.median_depth,
            pct(z.captured),
            z.opposition_months
        )
        .unwrap();
    }
    writeln!(s, "\n| Zone | Planet Nine there (median) | Gap between linked visits | TNOs | TNOs Rubin will not have |").unwrap();
    writeln!(s, "|---|---|---|---|---|").unwrap();
    for z in &r.zones {
        writeln!(
            s,
            "| {} | {:.0} AU, V = {:.1} | ≥ {:.1} h | {:.0} | {:.0} |",
            z.zone, z.median_dist_au, z.median_v, z.min_baseline_hr, z.tnos, z.tnos_new
        )
        .unwrap();
    }

    writeln!(s, "\n## Calendar\n").unwrap();
    writeln!(s, "Hours by month of opposition. A field stays within 70% of its peak motion for ±{:.0} days.\n", r.window_half_width_days).unwrap();
    writeln!(s, "| {} |", MONTHS.join(" | ")).unwrap();
    writeln!(s, "|{}", "---|".repeat(12)).unwrap();
    let cells: Vec<String> = r.hours_by_month.iter().map(|h| format!("{h:.0}")).collect();
    writeln!(s, "| {} |", cells.join(" | ")).unwrap();

    writeln!(s, "\n## If the orbit is someone else's\n").unwrap();
    writeln!(s, "The reference plan, scored as if each other solution were the truth, against the best plan for that solution at the same hours.\n").unwrap();
    writeln!(
        s,
        "| Truth | Left for this telescope | Reference plan captures | Best plan captures | Kept |"
    )
    .unwrap();
    writeln!(s, "|---|---|---|---|---|").unwrap();
    for x in &r.robustness {
        writeln!(
            s,
            "| {} | {} | {} | {} | {:.0}% |",
            x.prior,
            pct(x.unique),
            pct(x.reference_plan),
            pct(x.own_plan),
            100.0 * x.reference_plan / x.own_plan
        )
        .unwrap();
    }

    writeln!(s, "\n## Distant objects\n").unwrap();
    writeln!(
        s,
        "The reference campaign finds {:.0} TNOs, {:.0} of which Rubin will not have, {:.1} of those beyond 60 AU.\n",
        r.tno_plan_total, r.tno_plan_new, r.tno_plan_distant_new
    )
    .unwrap();
    writeln!(
        s,
        "The same telescope spending {:.0} h on nothing but distant objects Rubin cannot reach would find {:.0} new TNOs, {:.1} of them beyond 60 AU.",
        r.deep_program_hours, r.deep_program_tnos_new, r.deep_program_distant_new
    )
    .unwrap();
    s
}

// ---- figures ----------------------------------------------------------------

fn reference_curve(f: impl Fn(f64) -> (f64, f64)) -> Vec<(f64, f64)> {
    (0..=720).map(|k| f(k as f64 * 0.5)).collect()
}

fn galactic_curve(b_deg: f64) -> Vec<(f64, f64)> {
    use p9_core::coords::sky::{ecliptic_vec_to_equatorial_deg, galactic_to_ecliptic_matrix};
    reference_curve(|l| {
        let (l, b) = (l.to_radians(), b_deg.to_radians());
        let g = nalgebra::Vector3::new(b.cos() * l.cos(), b.cos() * l.sin(), b.sin());
        ecliptic_vec_to_equatorial_deg(&(galactic_to_ecliptic_matrix() * g))
    })
}

fn ecliptic_curve(beta_deg: f64) -> Vec<(f64, f64)> {
    use p9_core::coords::sky::ecliptic_to_equatorial_deg;
    reference_curve(|lon| ecliptic_to_equatorial_deg(lon, beta_deg))
}

/// The sky map: what is left for this telescope, and the reference plan.
pub fn sky_figure(r: &Report) -> String {
    let (w, h) = (1400.0, 760.0);
    let mut g = Svg::new(w, h);
    let (dec_lo, dec_hi) = (-60.0, 60.0);
    // East to the left: RA 24h at the left edge, 0h at the right.
    let x = Axis {
        d0: 360.0,
        d1: 0.0,
        p0: 90.0,
        p1: 1330.0,
    };
    let y = Axis {
        d0: dec_lo,
        d1: dec_hi,
        p0: 640.0,
        p1: 110.0,
    };
    g.bold(
        90.0,
        42.0,
        "Where a small space telescope should image for Planet Nine",
        24.0,
        svg::FG,
        "start",
    );
    g.text(
        90.0,
        68.0,
        &format!(
            "Shading: probability nobody else will collect. Outlined tiles: the {:.0} h reference campaign ({:.0} deg², captures {}).",
            r.reference_hours,
            r.reference_area_deg2,
            pct(r.reference_captured)
        ),
        14.0,
        svg::MUTED,
        "start",
    );
    g.rect(x.p0, y.p1, x.p1 - x.p0, y.p0 - y.p1, svg::PANEL, 1.0);

    let peak = r
        .sky
        .iter()
        .map(|c| c.unique / (c.ra_width_deg * c.dec_height_deg))
        .fold(0.0, f64::max);
    let cell_rect = |ra: f64, dec: f64, rw: f64, dh: f64| {
        let (x0, x1) = (x.at(ra + rw / 2.0), x.at(ra - rw / 2.0));
        let (y0, y1) = (y.at(dec + dh / 2.0), y.at(dec - dh / 2.0));
        (x0, y0, x1 - x0, y1 - y0)
    };
    for c in &r.sky {
        if c.dec_deg < dec_lo || c.dec_deg > dec_hi {
            continue;
        }
        let f = (c.unique / (c.ra_width_deg * c.dec_height_deg) / peak).sqrt();
        if f < 0.04 {
            continue;
        }
        let (rx, ry, rw, rh) = cell_rect(c.ra_deg, c.dec_deg, c.ra_width_deg, c.dec_height_deg);
        g.rect(
            rx,
            ry,
            rw + 0.4,
            rh + 0.4,
            svg::ORANGE,
            (0.15 + 0.85 * f).min(1.0),
        );
    }
    for t in &r.tiles {
        if t.dec_deg < dec_lo || t.dec_deg > dec_hi {
            continue;
        }
        let tier = r
            .tiers_s
            .iter()
            .position(|&s| s == t.integration_s)
            .unwrap_or(0);
        let (rx, ry, rw, rh) = cell_rect(t.ra_deg, t.dec_deg, t.ra_width_deg, t.dec_height_deg);
        g.outline(
            rx + 0.6,
            ry + 0.6,
            rw - 1.2,
            rh - 1.2,
            svg::TIER_COLOURS[tier],
            1.3,
        );
    }

    let project = |pts: Vec<(f64, f64)>| -> Vec<(f64, f64)> {
        pts.into_iter()
            .filter(|p| p.1 >= dec_lo && p.1 <= dec_hi)
            .map(|p| (x.at(p.0), y.at(p.1)))
            .collect()
    };
    g.path(&project(galactic_curve(0.0)), svg::PURPLE, 1.8, "", 60.0);
    g.path(
        &project(galactic_curve(15.0)),
        svg::PURPLE,
        0.9,
        "4 4",
        60.0,
    );
    g.path(
        &project(galactic_curve(-15.0)),
        svg::PURPLE,
        0.9,
        "4 4",
        60.0,
    );
    g.path(&project(ecliptic_curve(0.0)), svg::GREEN, 1.4, "7 5", 60.0);
    // Rubin: wide-fast-deep limit and the ecliptic spur.
    g.line(x.p0, y.at(12.0), x.p1, y.at(12.0), svg::BLUE, 1.6, "");
    let spur: Vec<(f64, f64)> = ecliptic_curve(10.0)
        .into_iter()
        .filter(|p| p.1 > 12.0)
        .map(|p| (p.0, p.1.min(30.0)))
        .collect();
    g.path(&project(spur), svg::BLUE, 1.6, "", 60.0);
    g.text(
        x.at(225.0),
        y.at(12.0) - 8.0,
        "Rubin main survey: south of +12°",
        12.0,
        svg::BLUE,
        "middle",
    );
    g.text(
        x.at(20.0),
        y.at(30.0) - 8.0,
        "Rubin ecliptic spur",
        12.0,
        svg::BLUE,
        "middle",
    );

    for ra_h in (0..=24).step_by(2) {
        let px = x.at(15.0 * ra_h as f64);
        g.line(px, y.p0, px, y.p0 + 5.0, svg::MUTED, 1.0, "");
        g.text(
            px,
            y.p0 + 20.0,
            &format!("{}h", ra_h % 24),
            12.0,
            svg::FG,
            "middle",
        );
    }
    for dec in (-60..=60).step_by(20) {
        let py = y.at(dec as f64);
        g.line(x.p0 - 5.0, py, x.p0, py, svg::MUTED, 1.0, "");
        g.line(x.p0, py, x.p1, py, svg::MUTED, 0.4, "2 6");
        g.text(
            x.p0 - 10.0,
            py + 4.0,
            &format!("{dec:+}°"),
            12.0,
            svg::FG,
            "end",
        );
    }
    g.text(
        (x.p0 + x.p1) / 2.0,
        y.p0 + 42.0,
        "Right ascension (east to the left)",
        13.0,
        svg::FG,
        "middle",
    );
    // Month of opposition along the top, from the ecliptic longitude at each RA.
    g.text(
        x.p0,
        y.p1 - 26.0,
        "At opposition in:",
        12.0,
        svg::MUTED,
        "start",
    );
    for ra_h in (1..24).step_by(2) {
        let ra = 15.0 * ra_h as f64;
        let (lon, _) = p9_core::coords::sky::equatorial_to_ecliptic_deg(ra, 0.0);
        let m = opposition_month(lon.rem_euclid(360.0));
        g.text(
            x.at(ra),
            y.p1 - 8.0,
            MONTHS[m as usize - 1],
            12.0,
            svg::MUTED,
            "middle",
        );
    }

    // Legend.
    let ly = 712.0;
    let mut lx = 90.0;
    g.text(lx, ly, "Visit length:", 12.0, svg::FG, "start");
    lx += 84.0;
    for (k, s) in r.tiers_s.iter().enumerate() {
        if !r.tiles.iter().any(|t| t.integration_s == *s) {
            continue;
        }
        g.outline(lx, ly - 11.0, 14.0, 12.0, svg::TIER_COLOURS[k], 1.6);
        g.text(
            lx + 20.0,
            ly,
            &format!("{s:.0} s (V {:.1})", r.tier_depths_ecliptic[k]),
            12.0,
            svg::FG,
            "start",
        );
        lx += 128.0;
    }
    lx += 20.0;
    for (label, colour, dash) in [
        ("Galactic plane, |b| = 15°", svg::PURPLE, ""),
        ("ecliptic", svg::GREEN, "7 5"),
        ("Rubin's northern limit", svg::BLUE, ""),
    ] {
        g.line(lx, ly - 4.0, lx + 26.0, ly - 4.0, colour, 1.8, dash);
        g.text(lx + 32.0, ly, label, 12.0, svg::FG, "start");
        lx += 46.0 + 7.2 * label.len() as f64;
    }
    g.text(90.0, 742.0, "Prior: Brown & Batygin (2021). Ground searches: ZTF, DES, Pan-STARRS1 as reproduced in this workspace, with crowding.", 12.0, svg::MUTED, "start");
    g.finish()
}

/// Probability captured against hours, for each planning stance and truth.
pub fn frontier_figure(r: &Report) -> String {
    let (w, h) = (1180.0, 640.0);
    let mut g = Svg::new(w, h);
    let top = r
        .policies
        .iter()
        .flat_map(|p| p.captured.iter().cloned())
        .fold(0.0, f64::max);
    let ymax = (top * 100.0 / 2.0).ceil() * 2.0;
    let x = Axis {
        d0: 0.0,
        d1: 8000.0,
        p0: 90.0,
        p1: 800.0,
    };
    let y = Axis {
        d0: 0.0,
        d1: ymax,
        p0: 560.0,
        p1: 100.0,
    };
    g.bold(
        90.0,
        42.0,
        "What telescope time buys",
        24.0,
        svg::FG,
        "start",
    );
    g.text(
        90.0,
        68.0,
        "Chance of finding Planet Nine (percent of the whole prior) against campaign length.",
        14.0,
        svg::MUTED,
        "start",
    );
    g.rect(x.p0, y.p1, x.p1 - x.p0, y.p0 - y.p1, svg::PANEL, 1.0);
    for k in 0..=8 {
        let hrs = 1000.0 * k as f64;
        g.line(x.at(hrs), y.p0, x.at(hrs), y.p1, svg::MUTED, 0.4, "2 6");
        g.text(
            x.at(hrs),
            y.p0 + 20.0,
            &format!("{hrs:.0}"),
            12.0,
            svg::FG,
            "middle",
        );
    }
    let mut v = 0.0;
    while v <= ymax + 1e-9 {
        g.line(x.p0, y.at(v), x.p1, y.at(v), svg::MUTED, 0.4, "2 6");
        g.text(
            x.p0 - 10.0,
            y.at(v) + 4.0,
            &format!("{v:.0}%"),
            12.0,
            svg::FG,
            "end",
        );
        v += 2.0;
    }
    g.text(
        (x.p0 + x.p1) / 2.0,
        y.p0 + 44.0,
        "Wall-clock hours",
        13.0,
        svg::FG,
        "middle",
    );

    let colours = [svg::RED, svg::ORANGE, svg::TEAL, svg::PURPLE];
    for (p, colour) in r.policies.iter().zip(colours) {
        let mut pts = vec![(x.at(0.0), y.at(0.0))];
        pts.extend(
            r.budgets_h
                .iter()
                .zip(&p.captured)
                .map(|(&hrs, &c)| (x.at(hrs), y.at(100.0 * c))),
        );
        g.path(
            &pts,
            colour,
            2.4,
            if p.cassini_phase { "6 4" } else { "" },
            1e9,
        );
        for pt in pts.iter().skip(1) {
            g.circle(pt.0, pt.1, 3.2, colour);
        }
    }
    // The reference campaign.
    let rx = x.at(r.reference_hours);
    g.line(
        rx,
        y.p0,
        rx,
        y.at(100.0 * r.reference_captured),
        svg::FG,
        1.2,
        "4 4",
    );
    g.circle(rx, y.at(100.0 * r.reference_captured), 5.0, svg::FG);
    g.text(
        rx + 8.0,
        y.at(100.0 * r.reference_captured) + 18.0,
        &format!("reference: {:.0} h", r.reference_hours),
        12.0,
        svg::FG,
        "start",
    );

    let mut ly = 120.0;
    g.bold(830.0, ly, "Planning stance", 13.0, svg::FG, "start");
    for (p, colour) in r.policies.iter().zip(colours) {
        ly += 24.0;
        g.line(
            830.0,
            ly - 4.0,
            858.0,
            ly - 4.0,
            colour,
            2.4,
            if p.cassini_phase { "6 4" } else { "" },
        );
        let name = if p.cassini_phase {
            "... and trust the Cassini phase"
        } else {
            p.policy
        };
        g.text(866.0, ly, name, 12.0, svg::FG, "start");
        ly += 16.0;
        g.text(
            866.0,
            ly,
            &format!("{} left to find", pct(p.unique)),
            11.0,
            svg::MUTED,
            "start",
        );
    }
    ly += 40.0;
    g.bold(830.0, ly, "If the true orbit is...", 13.0, svg::FG, "start");
    for row in &r.robustness {
        ly += 22.0;
        g.text(830.0, ly, row.prior, 12.0, svg::FG, "start");
        ly += 15.0;
        g.text(
            830.0,
            ly,
            &format!(
                "reference plan {} (best {})",
                pct(row.reference_plan),
                pct(row.own_plan)
            ),
            11.0,
            svg::MUTED,
            "start",
        );
    }
    g.finish()
}

/// Hours by month of opposition, stacked by zone.
pub fn calendar_figure(r: &Report) -> String {
    let (w, h) = (1100.0, 520.0);
    let mut g = Svg::new(w, h);
    g.bold(90.0, 42.0, "When to observe", 24.0, svg::FG, "start");
    g.text(
        90.0,
        68.0,
        &format!(
            "Reference campaign hours by month of opposition. Each field is usable for ±{:.0} days around it.",
            r.window_half_width_days
        ),
        14.0,
        svg::MUTED,
        "start",
    );
    let mut by: Vec<[f64; 12]> = vec![[0.0; 12]; Zone::ALL.len()];
    for t in &r.tiles {
        let k = Zone::ALL.iter().position(|z| z.label() == t.zone).unwrap();
        by[k][t.opposition_month as usize - 1] += t.hours;
    }
    let peak = (0..12)
        .map(|m| by.iter().map(|z| z[m]).sum::<f64>())
        .fold(0.0, f64::max);
    let ymax = (peak / 100.0).ceil() * 100.0;
    // July to June, so the winter season reads as one block.
    let order: Vec<usize> = (0..12).map(|k| (k + 6) % 12).collect();
    let x = Axis {
        d0: 0.0,
        d1: 12.0,
        p0: 90.0,
        p1: 800.0,
    };
    let y = Axis {
        d0: 0.0,
        d1: ymax,
        p0: 440.0,
        p1: 100.0,
    };
    g.rect(x.p0, y.p1, x.p1 - x.p0, y.p0 - y.p1, svg::PANEL, 1.0);
    let mut v = 0.0;
    while v <= ymax + 1e-9 {
        g.line(x.p0, y.at(v), x.p1, y.at(v), svg::MUTED, 0.4, "2 6");
        g.text(
            x.p0 - 10.0,
            y.at(v) + 4.0,
            &format!("{v:.0} h"),
            12.0,
            svg::FG,
            "end",
        );
        v += ymax / 4.0;
    }
    let colours = [svg::ORANGE, svg::RED, svg::TEAL, svg::BLUE];
    for (slot, &m) in order.iter().enumerate() {
        let mut base = 0.0;
        for (k, zone) in by.iter().enumerate() {
            if zone[m] > 0.0 {
                let (x0, x1) = (x.at(slot as f64 + 0.12), x.at(slot as f64 + 0.88));
                g.rect(
                    x0,
                    y.at(base + zone[m]),
                    x1 - x0,
                    y.at(base) - y.at(base + zone[m]),
                    colours[k],
                    0.85,
                );
                base += zone[m];
            }
        }
        g.text(
            x.at(slot as f64 + 0.5),
            y.p0 + 20.0,
            MONTHS[m],
            12.0,
            svg::FG,
            "middle",
        );
    }
    // 730 h in a month; the telescope's duty cycle is already in the hours.
    if 730.0 <= ymax {
        g.line(x.p0, y.at(730.0), x.p1, y.at(730.0), svg::MUTED, 1.2, "6 4");
        g.text(
            x.p1 - 6.0,
            y.at(730.0) - 6.0,
            "one month of telescope time",
            11.0,
            svg::MUTED,
            "end",
        );
    }
    let mut ly = 120.0;
    for (zone, colour) in Zone::ALL.iter().zip(colours) {
        if let Some(row) = r.zones.iter().find(|z| z.zone == zone.label()) {
            g.rect(830.0, ly - 11.0, 14.0, 12.0, colour, 0.85);
            g.text(852.0, ly, zone.label(), 12.0, svg::FG, "start");
            ly += 16.0;
            g.text(
                852.0,
                ly,
                &format!(
                    "{:.0} h, {:.0} deg², {}",
                    row.hours,
                    row.area_deg2,
                    pct(row.captured)
                ),
                11.0,
                svg::MUTED,
                "start",
            );
            ly += 26.0;
        }
    }
    g.finish()
}

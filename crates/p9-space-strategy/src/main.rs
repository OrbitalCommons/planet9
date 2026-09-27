//! Writes the strategy dataset, tables and figures.
//!
//!   cargo run --release -p p9-space-strategy
//!
//! Outputs (relative to the repository root):
//!   figures/space_strategy.json
//!   figures/p9_space_strategy_sky.svg
//!   figures/p9_space_strategy_frontier.svg
//!   figures/p9_space_strategy_calendar.svg
//!   reports/space-telescope/TABLES.md

use std::fs;

use p9_space_strategy::report::{build, calendar_figure, frontier_figure, sky_figure, tables};

fn main() {
    let report = build();
    fs::create_dir_all("figures").expect("figures/");
    fs::create_dir_all("reports/space-telescope").expect("reports/space-telescope/");
    fs::write(
        "figures/space_strategy.json",
        serde_json::to_string(&report).expect("serialise"),
    )
    .expect("write dataset");
    fs::write("figures/p9_space_strategy_sky.svg", sky_figure(&report)).expect("write sky");
    fs::write(
        "figures/p9_space_strategy_frontier.svg",
        frontier_figure(&report),
    )
    .expect("write frontier");
    fs::write(
        "figures/p9_space_strategy_calendar.svg",
        calendar_figure(&report),
    )
    .expect("write calendar");
    let t = tables(&report);
    fs::write("reports/space-telescope/TABLES.md", &t).expect("write tables");
    print!("{t}");
}

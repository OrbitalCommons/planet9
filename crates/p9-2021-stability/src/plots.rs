//! SVG plot generation for stability boundary figures.

use std::fmt::Write;

use crate::chirikov::critical_perihelion;
use crate::stability::{classify, StabilityClass};

/// Generate a stability map plot on the (a, q) plane.
pub fn stability_map_plot(
    a_min: f64,
    a_max: f64,
    q_min: f64,
    q_max: f64,
    n_a: usize,
    n_q: usize,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(8000);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Scattered Disk Stability Map</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    let da = (a_max - a_min) / n_a as f64;
    let dq = (q_max - q_min) / n_q as f64;
    let cell_w = plot_w / n_a as f64;
    let cell_h = plot_h / n_q as f64;

    for ia in 0..n_a {
        let a = a_min + (ia as f64 + 0.5) * da;
        for iq in 0..n_q {
            let q = q_min + (iq as f64 + 0.5) * dq;
            let class = classify(a, q);
            let color = match class {
                StabilityClass::Chaotic => "crimson",
                StabilityClass::Marginal => "gold",
                StabilityClass::Stable => "steelblue",
            };

            let x = margin + ia as f64 * cell_w;
            let y = margin + plot_h - (iq + 1) as f64 * cell_h;
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y:.1}" width="{cell_w:.1}" height="{cell_h:.1}" fill="{color}" opacity="0.7"/>"#,
            )
            .unwrap();
        }
    }

    // Critical perihelion curve
    let mut points = Vec::new();
    for ia in 0..=n_a * 2 {
        let a = a_min + ia as f64 * (a_max - a_min) / (n_a * 2) as f64;
        let q_crit = critical_perihelion(a);
        if q_crit >= q_min && q_crit <= q_max {
            let x = margin + ((a - a_min) / (a_max - a_min)) * plot_w;
            let y = margin + plot_h - ((q_crit - q_min) / (q_max - q_min)) * plot_h;
            points.push(format!("{x:.1},{y:.1}"));
        }
    }
    if points.len() > 1 {
        writeln!(
            svg,
            r#"<polyline points="{}" fill="none" stroke="black" stroke-width="2" stroke-dasharray="6,3"/>"#,
            points.join(" "),
        )
        .unwrap();
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Semi-major axis a (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90, 15, {})">Perihelion q (AU)</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Legend
    let lx = margin + plot_w - 110.0;
    let ly = margin + 15.0;
    for (i, (label, color)) in [
        ("Chaotic", "crimson"),
        ("Marginal", "gold"),
        ("Stable", "steelblue"),
    ]
    .iter()
    .enumerate()
    {
        let y_off = ly + i as f64 * 16.0;
        writeln!(
            svg,
            r#"<rect x="{lx}" y="{}" width="10" height="10" fill="{color}" opacity="0.7"/>"#,
            y_off - 8.0,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{y_off}" font-size="10" font-family="sans-serif">{label}</text>"#,
            lx + 14.0,
        )
        .unwrap();
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stability_map_valid_svg() {
        let svg = stability_map_plot(100.0, 500.0, 28.0, 50.0, 20, 20, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Scattered Disk Stability Map"));
    }

    #[test]
    fn test_stability_map_contains_all_colors() {
        let svg = stability_map_plot(100.0, 500.0, 28.0, 60.0, 30, 30, 600, 400);
        assert!(svg.contains("crimson"), "should have chaotic regions");
        assert!(svg.contains("steelblue"), "should have stable regions");
    }
}

//! SVG plot generation for constraint figures.
//!
//! - Figure 2: accepted (a₉, e₉) parameter space per mass
//! - Figure 5: accepted (i₉, ω₉) combinations
//! - Brightness curve along orbit

use std::fmt::Write;

use p9_core::constants::*;

use crate::parameter_grid::GridResult;

/// Generate a heat map of accepted/rejected parameter space in (a₉, e₉).
///
/// Green = accepted, red = rejected, darker = more extreme mass.
pub fn acceptance_heatmap(
    results: &[GridResult],
    mass_filter: f64,
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let filtered: Vec<&GridResult> = results
        .iter()
        .filter(|r| (r.point.mass_earth - mass_filter).abs() < 0.01)
        .collect();

    let mut svg = String::with_capacity(filtered.len() * 100 + 1000);
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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="18" text-anchor="middle" font-size="13" font-family="sans-serif">m = {} M⊕</text>"#,
        width as f64 / 2.0,
        mass_filter,
    )
    .unwrap();

    let a_min = 200.0_f64;
    let a_max = 2000.0_f64;
    let e_min = 0.1_f64;
    let e_max = 0.9_f64;

    let cell_w = plot_w / 19.0;
    let cell_h = plot_h / 9.0;

    for result in &filtered {
        let px = margin + (result.point.a - a_min) / (a_max - a_min) * plot_w;
        let py = margin + (1.0 - (result.point.e - e_min) / (e_max - e_min)) * plot_h;

        let color = if result.accepted { "green" } else { "#ffcccc" };

        let opacity = if result.accepted {
            0.3 + 0.7 * result.clustering_fraction
        } else {
            0.3
        };

        writeln!(
            svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{color}" opacity="{opacity:.2}"/>"#,
            px - cell_w / 2.0,
            py - cell_h / 2.0,
            cell_w,
            cell_h,
        )
        .unwrap();
    }

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">a₉ (AU)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">e₉</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a brightness curve plot (V magnitude vs true anomaly).
pub fn brightness_curve_plot(
    curve: &[(f64, f64, f64)],
    survey_depths: &[(&str, f64)],
    width: u32,
    height: u32,
) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let v_min = curve
        .iter()
        .map(|&(_, _, v)| v)
        .fold(f64::INFINITY, f64::min);
    let v_max = curve
        .iter()
        .map(|&(_, _, v)| v)
        .fold(f64::NEG_INFINITY, f64::max);
    let v_range = v_max - v_min;

    let mut svg = String::with_capacity(curve.len() * 60 + 2000);
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

    // Plot brightness curve
    let mut path = format!("M");
    for (i, &(nu, _, v)) in curve.iter().enumerate() {
        let px = margin + (nu + std::f64::consts::PI) / TWO_PI * plot_w;
        let py = margin + (v - v_min) / v_range * plot_h; // Brighter (lower V) at top
        if i == 0 {
            write!(path, "{px:.1},{py:.1}").unwrap();
        } else {
            write!(path, " L{px:.1},{py:.1}").unwrap();
        }
    }
    writeln!(
        svg,
        r#"<path d="{path}" fill="none" stroke="steelblue" stroke-width="2"/>"#,
    )
    .unwrap();

    // Survey depth lines
    let colors = ["red", "orange", "green", "blue", "purple"];
    for (i, &(name, depth)) in survey_depths.iter().enumerate() {
        if depth >= v_min && depth <= v_max {
            let py = margin + (depth - v_min) / v_range * plot_h;
            let color = colors[i % colors.len()];
            writeln!(
                svg,
                r#"<line x1="{margin}" y1="{py:.1}" x2="{x2:.1}" y2="{py:.1}" stroke="{color}" stroke-dasharray="4,4"/>"#,
                x2 = margin + plot_w,
            )
            .unwrap();
            writeln!(
                svg,
                r#"<text x="{x:.1}" y="{y:.1}" font-size="9" font-family="sans-serif" fill="{color}">{name}</text>"#,
                x = margin + plot_w + 3.0,
                y = py + 3.0,
            )
            .unwrap();
        }
    }

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">True anomaly (rad)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 8.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="12" y="{yc:.1}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90,12,{yc:.1})">V mag</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_grid::GridPoint;

    #[test]
    fn test_acceptance_heatmap() {
        let results = vec![
            GridResult {
                point: GridPoint {
                    mass_earth: 10.0,
                    a: 700.0,
                    e: 0.6,
                    perihelion: 280.0,
                },
                clustering_fraction: 0.8,
                high_perihelion_fraction: 0.2,
                n_survivors: 15,
                n_total: 400,
                accepted: true,
            },
            GridResult {
                point: GridPoint {
                    mass_earth: 10.0,
                    a: 300.0,
                    e: 0.2,
                    perihelion: 240.0,
                },
                clustering_fraction: 0.1,
                high_perihelion_fraction: 0.0,
                n_survivors: 3,
                n_total: 400,
                accepted: false,
            },
        ];

        let svg = acceptance_heatmap(&results, 10.0, 500, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_brightness_curve_plot() {
        let curve = vec![
            (-3.0, 280.0, 22.0),
            (-1.0, 400.0, 24.0),
            (0.0, 280.0, 22.0),
            (1.0, 400.0, 24.0),
            (3.0, 1120.0, 28.0),
        ];
        let depths = crate::detection_limits::survey_depths();
        let svg = brightness_curve_plot(&curve, &depths, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("path"));
    }
}

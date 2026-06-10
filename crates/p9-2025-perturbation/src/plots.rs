//! SVG plot generation for perturbation theory figures.

use std::fmt::Write;

use crate::stability_boundary::BoundaryCurve;

/// Plot the three stability boundary curves on the (a, q) plane.
pub fn stability_boundary_plot(curve: &BoundaryCurve, width: u32, height: u32) -> String {
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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Scattered Disk Stability Boundaries</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    if curve.a_values.is_empty() {
        writeln!(svg, "</svg>").unwrap();
        return svg;
    }

    let a_min = curve.a_values[0];
    let a_max = *curve.a_values.last().unwrap();
    let a_range = (a_max - a_min).max(1.0);

    let q_min = 25.0_f64;
    let q_max = 55.0_f64;
    let q_range = q_max - q_min;

    let to_x = |a: f64| margin + ((a - a_min) / a_range) * plot_w;
    let to_y = |q: f64| margin + plot_h - ((q - q_min) / q_range) * plot_h;

    // Draw the three boundary curves
    let curves: [(&[f64], &str, &str); 3] = [
        (&curve.q_chaotic, "crimson", "Chaotic layer"),
        (&curve.q_diffusion, "darkorange", "1:j diffusion barrier"),
        (&curve.q_comb, "steelblue", "Fine comb"),
    ];

    for (q_values, color, _label) in &curves {
        let points: Vec<String> = curve
            .a_values
            .iter()
            .zip(q_values.iter())
            .filter(|(_, &q)| q >= q_min && q <= q_max)
            .map(|(&a, &q)| format!("{:.1},{:.1}", to_x(a), to_y(q)))
            .collect();

        if points.len() > 1 {
            writeln!(
                svg,
                r#"<polyline points="{}" fill="none" stroke="{color}" stroke-width="2"/>"#,
                points.join(" "),
            )
            .unwrap();
        }
    }

    // Effective boundary (thick black dashed)
    let eff_points: Vec<String> = curve
        .a_values
        .iter()
        .zip(curve.q_effective.iter())
        .filter(|(_, &q)| q >= q_min && q <= q_max)
        .map(|(&a, &q)| format!("{:.1},{:.1}", to_x(a), to_y(q)))
        .collect();

    if eff_points.len() > 1 {
        writeln!(
            svg,
            r#"<polyline points="{}" fill="none" stroke="black" stroke-width="2.5" stroke-dasharray="8,4"/>"#,
            eff_points.join(" "),
        )
        .unwrap();
    }

    // Legend
    let lx = margin + 10.0;
    let ly = margin + 15.0;
    for (i, (_, color, label)) in curves.iter().enumerate() {
        let y_off = ly + i as f64 * 18.0;
        writeln!(
            svg,
            r#"<line x1="{lx}" y1="{y_off}" x2="{}" y2="{y_off}" stroke="{color}" stroke-width="2"/>"#,
            lx + 20.0,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif">{label}</text>"#,
            lx + 25.0,
            y_off + 4.0,
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

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boundary_plot_valid_svg() {
        let curve = BoundaryCurve::compute(100.0, 500.0, 50);
        let svg = stability_boundary_plot(&curve, 700, 500);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("crimson"));
        assert!(svg.contains("darkorange"));
        assert!(svg.contains("steelblue"));
    }

    #[test]
    fn test_boundary_plot_empty_curve() {
        let curve = BoundaryCurve {
            a_values: vec![],
            q_chaotic: vec![],
            q_diffusion: vec![],
            q_comb: vec![],
            q_effective: vec![],
            regimes: vec![],
        };
        let svg = stability_boundary_plot(&curve, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }
}

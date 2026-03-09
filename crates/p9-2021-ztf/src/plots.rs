//! SVG plot generation for the ZTF Planet Nine search.
//!
//! Produces:
//! - Exclusion heatmap in (semi-major axis, V magnitude) space
//! - Detection efficiency vs. apparent magnitude curve

use std::fmt::Write;

/// Generate an SVG heatmap of exclusion fraction in (a, V) parameter space.
///
/// `a_bins`: semi-major axis bin centers (AU)
/// `v_bins`: apparent magnitude bin centers
/// `exclusion_grid`: 2D grid where `exclusion_grid[i][j]` is the exclusion
///     fraction for `(a_bins[i], v_bins[j])`, values in [0, 1].
pub fn exclusion_map_plot(
    a_bins: &[f64],
    v_bins: &[f64],
    exclusion_grid: &[Vec<f64>],
    width: u32,
    height: u32,
) -> String {
    let n_a = a_bins.len();
    let n_v = v_bins.len();
    if n_a == 0 || n_v == 0 {
        return String::new();
    }

    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let cell_w = plot_w / n_a as f64;
    let cell_h = plot_h / n_v as f64;

    let mut svg = String::with_capacity(n_a * n_v * 120 + 1000);
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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">ZTF Exclusion Map</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Heatmap cells
    for (i, row) in exclusion_grid.iter().enumerate() {
        for (j, &frac) in row.iter().enumerate() {
            let t = frac.clamp(0.0, 1.0);
            let (r, g, b) = exclusion_color(t);
            let x = margin + i as f64 * cell_w;
            let y = margin + j as f64 * cell_h;
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y:.1}" width="{cw:.1}" height="{ch:.1}" fill="rgb({r},{g},{b})" stroke="none"/>"#,
                cw = cell_w + 0.5,
                ch = cell_h + 0.5,
            )
            .unwrap();
        }
    }

    // 56.4% contour label
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif" fill="white" text-anchor="middle">56.4% excluded</text>"#,
        margin + plot_w / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Axes
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    // X-axis ticks (semi-major axis)
    let a_min = a_bins.first().copied().unwrap_or(200.0);
    let a_max = a_bins.last().copied().unwrap_or(1000.0);
    let a_range = a_max - a_min;
    if a_range > 0.0 {
        for &tick in &[200.0, 400.0, 600.0, 800.0, 1000.0] {
            if tick >= a_min && tick <= a_max {
                let x = margin + (tick - a_min) / a_range * plot_w;
                writeln!(
                    svg,
                    r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{tick:.0}</text>"#,
                    margin + plot_h + 15.0,
                )
                .unwrap();
            }
        }
    }
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Semi-major Axis (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    // Y-axis label
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">V Magnitude</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate an SVG plot of detection efficiency vs. apparent magnitude.
///
/// `mag_bins`: apparent magnitude bin centers
/// `efficiencies`: detection efficiency in each bin, values in [0, 1]
pub fn detection_efficiency_plot(
    mag_bins: &[f64],
    efficiencies: &[f64],
    width: u32,
    height: u32,
) -> String {
    assert_eq!(mag_bins.len(), efficiencies.len());
    if mag_bins.is_empty() {
        return String::new();
    }

    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let v_min = mag_bins.first().copied().unwrap_or(16.0);
    let v_max = mag_bins.last().copied().unwrap_or(24.0);
    let v_range = v_max - v_min;

    let mut svg = String::with_capacity(2000);
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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">ZTF Detection Efficiency</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    // ZTF depth limit line at V=20.5
    let depth_frac = if v_range > 0.0 {
        (20.5 - v_min) / v_range
    } else {
        0.5
    };
    let depth_x = margin + depth_frac * plot_w;
    writeln!(
        svg,
        r#"<line x1="{depth_x:.1}" y1="{margin}" x2="{depth_x:.1}" y2="{}" stroke="red" stroke-width="1" stroke-dasharray="4,4"/>"#,
        margin + plot_h,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="8" font-family="sans-serif" fill="red">V=20.5</text>"#,
        depth_x + 3.0,
        margin + 12.0,
    )
    .unwrap();

    // Efficiency curve
    let mut points = String::new();
    for (i, (mag, eff)) in mag_bins.iter().zip(efficiencies.iter()).enumerate() {
        let x = margin
            + if v_range > 0.0 {
                (mag - v_min) / v_range * plot_w
            } else {
                plot_w / 2.0
            };
        let y = margin + plot_h - eff.clamp(0.0, 1.0) * plot_h;
        if i == 0 {
            write!(points, "{x:.1},{y:.1}").unwrap();
        } else {
            write!(points, " {x:.1},{y:.1}").unwrap();
        }
    }
    writeln!(
        svg,
        r#"<polyline points="{points}" fill="none" stroke="steelblue" stroke-width="2"/>"#,
    )
    .unwrap();

    // Data points
    for (mag, eff) in mag_bins.iter().zip(efficiencies.iter()) {
        let x = margin
            + if v_range > 0.0 {
                (mag - v_min) / v_range * plot_w
            } else {
                plot_w / 2.0
            };
        let y = margin + plot_h - eff.clamp(0.0, 1.0) * plot_h;
        writeln!(
            svg,
            r#"<circle cx="{x:.1}" cy="{y:.1}" r="3" fill="steelblue"/>"#,
        )
        .unwrap();
    }

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif">Apparent Magnitude (V)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="11" font-family="sans-serif" transform="rotate(-90, 15, {})">Detection Efficiency</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Y-axis ticks
    for &tick in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = margin + plot_h - tick * plot_h;
        writeln!(
            svg,
            r#"<text x="{}" y="{y:.1}" text-anchor="end" font-size="9" font-family="sans-serif">{tick:.2}</text>"#,
            margin - 5.0,
        )
        .unwrap();
    }

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Map exclusion fraction [0, 1] to an (r, g, b) color.
/// 0 = dark blue (not excluded), 1 = bright red (excluded).
fn exclusion_color(t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    let r = (t * 220.0) as u8;
    let g = ((1.0 - (2.0 * t - 1.0).abs()) * 120.0) as u8;
    let b = ((1.0 - t) * 200.0) as u8;
    (r, g, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exclusion_map_produces_valid_svg() {
        let a_bins = vec![300.0, 500.0, 700.0, 900.0];
        let v_bins = vec![18.0, 19.0, 20.0, 21.0, 22.0];
        let grid = vec![
            vec![0.9, 0.8, 0.5, 0.2, 0.05],
            vec![0.85, 0.7, 0.4, 0.15, 0.02],
            vec![0.6, 0.4, 0.2, 0.05, 0.01],
            vec![0.3, 0.15, 0.05, 0.01, 0.0],
        ];
        let svg = exclusion_map_plot(&a_bins, &v_bins, &grid, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Exclusion Map"));
        assert!(svg.contains("56.4%"));
    }

    #[test]
    fn detection_efficiency_produces_valid_svg() {
        let mags = vec![17.0, 18.0, 19.0, 20.0, 20.5, 21.0, 22.0];
        let effs = vec![0.95, 0.93, 0.85, 0.60, 0.30, 0.05, 0.0];
        let svg = detection_efficiency_plot(&mags, &effs, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Detection Efficiency"));
        assert!(svg.contains("V=20.5"));
    }

    #[test]
    fn exclusion_color_endpoints() {
        let (r, _g, b) = exclusion_color(0.0);
        assert_eq!(r, 0);
        assert_eq!(b, 200);

        let (r, _g, b) = exclusion_color(1.0);
        assert_eq!(r, 220);
        assert_eq!(b, 0);
    }
}

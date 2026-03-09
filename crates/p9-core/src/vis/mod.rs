//! SVG visualization utilities for orbital dynamics.
//!
//! Generates simple SVG plots for phase-space portraits and orbital element distributions.

use std::fmt::Write;

/// Generate an SVG phase-space portrait from a Hamiltonian grid.
///
/// Produces a filled contour-style heatmap of H(e, Δϖ) values.
/// `e_vals`: eccentricity grid points
/// `dvarpi_vals`: Δϖ grid points (radians)
/// `portrait`: 2D grid of Hamiltonian values, portrait[i][j] = H(e_i, dvarpi_j)
pub fn phase_portrait_svg(
    e_vals: &[f64],
    dvarpi_vals: &[f64],
    portrait: &[Vec<f64>],
    width: u32,
    height: u32,
) -> String {
    let n_e = e_vals.len();
    let n_dv = dvarpi_vals.len();
    if n_e == 0 || n_dv == 0 {
        return String::new();
    }

    // Find min/max Hamiltonian values for color scaling
    let mut h_min = f64::INFINITY;
    let mut h_max = f64::NEG_INFINITY;
    for row in portrait {
        for &v in row {
            if v < h_min {
                h_min = v;
            }
            if v > h_max {
                h_max = v;
            }
        }
    }
    let h_range = h_max - h_min;

    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let cell_w = plot_w / n_dv as f64;
    let cell_h = plot_h / n_e as f64;

    let mut svg = String::with_capacity(n_e * n_dv * 100);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#,
    )
    .unwrap();

    // Heatmap cells
    for (i, row) in portrait.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            let t = if h_range > 0.0 {
                (val - h_min) / h_range
            } else {
                0.5
            };
            let (r, g, b) = viridis_color(t);
            let x = margin + j as f64 * cell_w;
            // Flip y so eccentricity increases upward
            let y = margin + (n_e - 1 - i) as f64 * cell_h;
            writeln!(
                svg,
                r#"<rect x="{x:.1}" y="{y:.1}" width="{cw:.1}" height="{ch:.1}" fill="rgb({r},{g},{b})" stroke="none"/>"#,
                cw = cell_w + 0.5,
                ch = cell_h + 0.5,
            )
            .unwrap();
        }
    }

    // Axes labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="14" font-family="sans-serif">Δϖ (rad)</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 10.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="14" font-family="sans-serif" transform="rotate(-90, 15, {})">e</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Tick marks on x-axis
    let x_ticks = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    let dv_min = dvarpi_vals
        .first()
        .copied()
        .unwrap_or(-std::f64::consts::PI);
    let dv_max = dvarpi_vals.last().copied().unwrap_or(std::f64::consts::PI);
    let dv_range = dv_max - dv_min;
    for &tick in &x_ticks {
        if tick >= dv_min && tick <= dv_max {
            let frac = (tick - dv_min) / dv_range;
            let x = margin + frac * plot_w;
            let y1 = margin + plot_h;
            let y2 = y1 + 5.0;
            writeln!(
                svg,
                r#"<line x1="{x:.1}" y1="{y1:.1}" x2="{x:.1}" y2="{y2:.1}" stroke="black"/>"#,
            )
            .unwrap();
            writeln!(
                svg,
                r#"<text x="{x:.1}" y="{y3:.1}" text-anchor="middle" font-size="10" font-family="sans-serif">{tick:.0}</text>"#,
                y3 = y2 + 12.0,
            )
            .unwrap();
        }
    }

    // Tick marks on y-axis
    let e_min = e_vals.first().copied().unwrap_or(0.0);
    let e_max = e_vals.last().copied().unwrap_or(1.0);
    let e_range = e_max - e_min;
    for tick_i in 0..=4 {
        let tick = e_min + (tick_i as f64 / 4.0) * e_range;
        let frac = (tick - e_min) / e_range;
        let y = margin + (1.0 - frac) * plot_h;
        let x1 = margin - 5.0;
        writeln!(
            svg,
            r#"<line x1="{x1:.1}" y1="{y:.1}" x2="{margin:.1}" y2="{y:.1}" stroke="black"/>"#,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{x2:.1}" y="{y2:.1}" text-anchor="end" font-size="10" font-family="sans-serif">{tick:.2}</text>"#,
            x2 = margin - 8.0,
            y2 = y + 3.0,
        )
        .unwrap();
    }

    // Border
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate an SVG scatter plot of orbital elements (a vs e, or similar).
pub fn scatter_svg(
    x_data: &[f64],
    y_data: &[f64],
    x_label: &str,
    y_label: &str,
    width: u32,
    height: u32,
) -> String {
    assert_eq!(x_data.len(), y_data.len());
    if x_data.is_empty() {
        return String::new();
    }

    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let x_min = x_data.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = x_data.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_min = y_data.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = y_data.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    let x_range = if (x_max - x_min).abs() < 1e-30 {
        1.0
    } else {
        x_max - x_min
    };
    let y_range = if (y_max - y_min).abs() < 1e-30 {
        1.0
    } else {
        y_max - y_min
    };

    let mut svg = String::with_capacity(x_data.len() * 60 + 500);
    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"#,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<rect width="{width}" height="{height}" fill="white"/>"#,
    )
    .unwrap();

    // Points
    for (x, y) in x_data.iter().zip(y_data.iter()) {
        let px = margin + (x - x_min) / x_range * plot_w;
        let py = margin + (1.0 - (y - y_min) / y_range) * plot_h;
        writeln!(
            svg,
            r#"<circle cx="{px:.1}" cy="{py:.1}" r="2" fill="steelblue" opacity="0.6"/>"#,
        )
        .unwrap();
    }

    // Axes
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="14" font-family="sans-serif">{x_label}</text>"#,
        margin + plot_w / 2.0,
        height as f64 - 10.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{yc}" text-anchor="middle" font-size="14" font-family="sans-serif" transform="rotate(-90, 15, {yc})">{y_label}</text>"#,
        yc = margin + plot_h / 2.0,
    )
    .unwrap();

    // Border
    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Approximate viridis colormap: maps t in [0, 1] to (r, g, b) in [0, 255].
fn viridis_color(t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    // Simplified piecewise-linear approximation of viridis
    let (r, g, b) = if t < 0.25 {
        let s = t / 0.25;
        (
            68.0 + s * (49.0 - 68.0),
            1.0 + s * (54.0 - 1.0),
            84.0 + s * (149.0 - 84.0),
        )
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (
            49.0 + s * (33.0 - 49.0),
            54.0 + s * (145.0 - 54.0),
            149.0 + s * (140.0 - 149.0),
        )
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (
            33.0 + s * (144.0 - 33.0),
            145.0 + s * (201.0 - 145.0),
            140.0 + s * (93.0 - 140.0),
        )
    } else {
        let s = (t - 0.75) / 0.25;
        (
            144.0 + s * (253.0 - 144.0),
            201.0 + s * (231.0 - 201.0),
            93.0 + s * (37.0 - 93.0),
        )
    };
    (r as u8, g as u8, b as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_portrait_svg_produces_valid_svg() {
        let e_vals = vec![0.1, 0.3, 0.5];
        let dvarpi_vals = vec![-3.0, -1.0, 1.0, 3.0];
        let portrait = vec![
            vec![1.0, 2.0, 3.0, 2.0],
            vec![2.0, 4.0, 4.0, 2.0],
            vec![1.0, 3.0, 3.0, 1.0],
        ];
        let svg = phase_portrait_svg(&e_vals, &dvarpi_vals, &portrait, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("rect"));
    }

    #[test]
    fn scatter_svg_produces_valid_svg() {
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = vec![0.1, 0.5, 0.3, 0.8];
        let svg = scatter_svg(&x, &y, "a (AU)", "e", 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("circle"));
    }

    #[test]
    fn viridis_endpoints() {
        let (r0, g0, b0) = viridis_color(0.0);
        assert_eq!(r0, 68);
        assert_eq!(g0, 1);
        assert_eq!(b0, 84);

        let (r1, g1, b1) = viridis_color(1.0);
        assert_eq!(r1, 253);
        assert_eq!(g1, 231);
        assert_eq!(b1, 37);
    }
}

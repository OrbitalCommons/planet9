//! SVG plot generation for IOC self-gravity figures.

use std::fmt::Write;

use crate::vzlk::VzlkPoint;

/// Generate a vZLK phase portrait plot on the (omega, q) plane.
pub fn vzlk_portrait_plot(points: &[VzlkPoint], width: u32, height: u32) -> String {
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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">vZLK Phase Portrait (omega, q)</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    if !points.is_empty() {
        let h_min = points
            .iter()
            .map(|p| p.h_value)
            .fold(f64::INFINITY, f64::min);
        let h_max = points
            .iter()
            .map(|p| p.h_value)
            .fold(f64::NEG_INFINITY, f64::max);
        let h_range = (h_max - h_min).max(1e-30);
        let q_min = points.iter().map(|p| p.q).fold(f64::INFINITY, f64::min);
        let q_max = points.iter().map(|p| p.q).fold(f64::NEG_INFINITY, f64::max);
        let q_range = (q_max - q_min).max(1e-10);
        let omega_max = std::f64::consts::TAU;

        for p in points {
            let x = margin + (p.omega / omega_max) * plot_w;
            let y = margin + plot_h - ((p.q - q_min) / q_range) * plot_h;
            let t = (p.h_value - h_min) / h_range;

            let r = (68.0 + t * 187.0) as u8;
            let g = (1.0 + t * 219.0) as u8;
            let b = (84.0 + t * (167.0 - 84.0 * t)) as u8;

            writeln!(
                svg,
                r#"<circle cx="{x:.1}" cy="{y:.1}" r="1.5" fill="rgb({r},{g},{b})" opacity="0.6"/>"#,
            )
            .unwrap();
        }
    }

    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">omega (rad)</text>"#,
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
    fn test_vzlk_plot_valid_svg() {
        let points = vec![
            VzlkPoint {
                omega: 0.5,
                q: 100.0,
                h_value: -1.0,
            },
            VzlkPoint {
                omega: 1.5,
                q: 200.0,
                h_value: -0.5,
            },
            VzlkPoint {
                omega: 3.0,
                q: 150.0,
                h_value: -0.8,
            },
        ];
        let svg = vzlk_portrait_plot(&points, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_vzlk_plot_empty() {
        let svg = vzlk_portrait_plot(&[], 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }
}

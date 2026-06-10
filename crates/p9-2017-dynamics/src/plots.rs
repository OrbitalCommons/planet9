//! SVG plot generation for dynamical evolution figures.

use std::fmt::Write;

use crate::phase_space::PhasePoint;
use crate::resonance::Resonance;

/// Generate an e vs delta_varpi phase portrait plot.
pub fn phase_portrait_plot(points: &[PhasePoint], width: u32, height: u32) -> String {
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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Phase Portrait: e vs Delta varpi</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
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

        let dv_max = std::f64::consts::TAU;
        let e_min = points.iter().map(|p| p.e).fold(f64::INFINITY, f64::min);
        let e_max = points.iter().map(|p| p.e).fold(f64::NEG_INFINITY, f64::max);
        let e_range = (e_max - e_min).max(1e-10);

        for p in points {
            let x = margin + (p.delta_varpi / dv_max) * plot_w;
            let y = margin + plot_h - ((p.e - e_min) / e_range) * plot_h;
            let t = (p.h_value - h_min) / h_range;

            // Viridis-inspired color mapping
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

    // Axis labels
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Delta varpi (rad)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="15" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90, 15, {})">Eccentricity</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a resonance structure diagram showing key MMR locations.
pub fn resonance_diagram(resonances: &[Resonance], a9: f64, width: u32, height: u32) -> String {
    let margin = 60.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let mut svg = String::with_capacity(4000);
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
        r#"<text x="{}" y="25" text-anchor="middle" font-size="14" font-family="sans-serif" font-weight="bold">Mean-Motion Resonances with Planet Nine</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    writeln!(
        svg,
        r#"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black"/>"#,
    )
    .unwrap();

    let a_min = 100.0_f64;
    let a_max = a9 + 100.0;
    let a_range = a_max - a_min;

    let colors = [
        "steelblue",
        "crimson",
        "forestgreen",
        "darkorange",
        "purple",
        "teal",
        "brown",
    ];

    for (idx, res) in resonances.iter().enumerate() {
        let x_center = margin + ((res.a_nominal - a_min) / a_range) * plot_w;
        let w = (res.delta_a / a_range) * plot_w * 2.0;
        let color = colors[idx % colors.len()];

        writeln!(
            svg,
            r#"<rect x="{:.1}" y="{margin}" width="{:.1}" height="{plot_h}" fill="{color}" opacity="0.2"/>"#,
            x_center - w / 2.0,
            w,
        )
        .unwrap();

        writeln!(
            svg,
            r#"<line x1="{x_center:.1}" y1="{margin}" x2="{x_center:.1}" y2="{}" stroke="{color}" stroke-width="1.5"/>"#,
            margin + plot_h,
        )
        .unwrap();

        writeln!(
            svg,
            r#"<text x="{x_center:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{}:{}</text>"#,
            margin - 5.0,
            res.p,
            res.q,
        )
        .unwrap();
    }

    // Axis label
    writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif">Semi-major axis (AU)</text>"#,
        width as f64 / 2.0,
        margin + plot_h + 35.0,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamiltonian::SecularHamiltonianParams;
    use crate::phase_space::{compute_phase_portrait, PhasePortraitConfig};
    use crate::resonance::key_resonances_700;

    #[test]
    fn test_phase_portrait_plot_valid_svg() {
        let config = PhasePortraitConfig::low_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_phase_portrait(&config, &params);
        let svg = phase_portrait_plot(&points, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_resonance_diagram_valid_svg() {
        let res = key_resonances_700();
        let svg = resonance_diagram(&res, 700.0, 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("1:1"));
    }

    #[test]
    fn test_phase_portrait_empty() {
        let svg = phase_portrait_plot(&[], 600, 400);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
    }
}

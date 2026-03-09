//! SVG plot generation for Brown & Batygin (2021) figures.
//!
//! - Corner-plot style showing semi-major axis vs eccentricity posterior
//! - Reference population scatter of semi-major axis vs V magnitude

use std::fmt::Write;

use crate::posterior::P9Posterior;
use crate::reference_population::ReferenceP9;

/// Generate a corner-plot style SVG showing a vs e posterior contours.
///
/// Draws the 1-sigma and 2-sigma ellipses of the joint (a, e) posterior,
/// with the median marked.
pub fn posterior_plot(posterior: &P9Posterior, width: u32, height: u32) -> String {
    let margin = 70.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let a_min = 200.0_f64;
    let a_max = 700.0_f64;
    let e_min = 0.0_f64;
    let e_max = 0.6_f64;

    let a_to_x = |a: f64| margin + (a - a_min) / (a_max - a_min) * plot_w;
    let e_to_y = |e: f64| margin + plot_h - (e - e_min) / (e_max - e_min) * plot_h;

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

    // Title
    writeln!(
        svg,
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Planet Nine Posterior: a vs e (Brown &amp; Batygin 2021)</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r##"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="#f8f8f8" stroke="black"/>"##,
    )
    .unwrap();

    // Grid lines
    for a_tick in [200, 300, 400, 500, 600, 700] {
        let x = a_to_x(a_tick as f64);
        writeln!(
            svg,
            r##"<line x1="{x:.1}" y1="{margin}" x2="{x:.1}" y2="{}" stroke="#ddd" stroke-width="0.5"/>"##,
            margin + plot_h,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{a_tick}</text>"#,
            margin + plot_h + 15.0,
        )
        .unwrap();
    }
    for e_tick_10 in [0, 1, 2, 3, 4, 5, 6] {
        let e_val = e_tick_10 as f64 / 10.0;
        let y = e_to_y(e_val);
        writeln!(
            svg,
            r##"<line x1="{margin}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="#ddd" stroke-width="0.5"/>"##,
            margin + plot_w,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{y:.1}" text-anchor="end" font-size="9" font-family="sans-serif" dy="3">{:.1}</text>"#,
            margin - 5.0,
            e_val,
        )
        .unwrap();
    }

    // 2-sigma contour ellipse
    let cx = a_to_x(posterior.a.median);
    let cy = e_to_y(posterior.e.median);
    let rx_2sig =
        2.0 * (posterior.a.sigma_upper + posterior.a.sigma_lower) / 2.0 / (a_max - a_min) * plot_w;
    let ry_2sig =
        2.0 * (posterior.e.sigma_upper + posterior.e.sigma_lower) / 2.0 / (e_max - e_min) * plot_h;
    writeln!(
        svg,
        r#"<ellipse cx="{cx:.1}" cy="{cy:.1}" rx="{rx_2sig:.1}" ry="{ry_2sig:.1}" fill="steelblue" fill-opacity="0.15" stroke="steelblue" stroke-width="1" stroke-dasharray="4,4"/>"#,
    )
    .unwrap();

    // 1-sigma contour ellipse
    let rx_1sig =
        (posterior.a.sigma_upper + posterior.a.sigma_lower) / 2.0 / (a_max - a_min) * plot_w;
    let ry_1sig =
        (posterior.e.sigma_upper + posterior.e.sigma_lower) / 2.0 / (e_max - e_min) * plot_h;
    writeln!(
        svg,
        r#"<ellipse cx="{cx:.1}" cy="{cy:.1}" rx="{rx_1sig:.1}" ry="{ry_1sig:.1}" fill="steelblue" fill-opacity="0.3" stroke="steelblue" stroke-width="1.5"/>"#,
    )
    .unwrap();

    // Median point
    writeln!(
        svg,
        r#"<circle cx="{cx:.1}" cy="{cy:.1}" r="4" fill="red" stroke="black" stroke-width="1"/>"#,
    )
    .unwrap();

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
        r#"<text x="15" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90, 15, {})">Eccentricity e</text>"#,
        margin + plot_h / 2.0,
        margin + plot_h / 2.0,
    )
    .unwrap();

    // Legend
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="10" font-family="sans-serif" fill="steelblue">Median: a = {:.0} AU, e = {:.2}</text>"#,
        margin + 10.0,
        margin + 15.0,
        posterior.a.median,
        posterior.e.median,
    )
    .unwrap();

    writeln!(svg, "</svg>").unwrap();
    svg
}

/// Generate a scatter plot of semi-major axis vs V magnitude for the reference population.
///
/// Shows how brightness varies across the posterior, helping identify
/// the most detectable regions of parameter space.
pub fn reference_population_plot(population: &[ReferenceP9], width: u32, height: u32) -> String {
    let margin = 70.0_f64;
    let plot_w = width as f64 - 2.0 * margin;
    let plot_h = height as f64 - 2.0 * margin;

    let a_min = 200.0_f64;
    let a_max = 800.0_f64;
    let v_min = 18.0_f64;
    let v_max = 28.0_f64;

    let a_to_x = |a: f64| margin + (a - a_min) / (a_max - a_min) * plot_w;
    let v_to_y = |v: f64| margin + (v - v_min) / (v_max - v_min) * plot_h;

    let mut svg = String::with_capacity(4000 + population.len() * 120);
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
        r#"<text x="{}" y="20" text-anchor="middle" font-size="13" font-family="sans-serif">Reference Population: a vs V magnitude</text>"#,
        width as f64 / 2.0,
    )
    .unwrap();

    // Plot frame
    writeln!(
        svg,
        r##"<rect x="{margin}" y="{margin}" width="{plot_w}" height="{plot_h}" fill="#f8f8f8" stroke="black"/>"##,
    )
    .unwrap();

    // Grid lines
    for a_tick in [200, 300, 400, 500, 600, 700, 800] {
        let x = a_to_x(a_tick as f64);
        writeln!(
            svg,
            r##"<line x1="{x:.1}" y1="{margin}" x2="{x:.1}" y2="{}" stroke="#ddd" stroke-width="0.5"/>"##,
            margin + plot_h,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{x:.1}" y="{}" text-anchor="middle" font-size="9" font-family="sans-serif">{a_tick}</text>"#,
            margin + plot_h + 15.0,
        )
        .unwrap();
    }
    for v_tick in [18, 20, 22, 24, 26, 28] {
        let y = v_to_y(v_tick as f64);
        writeln!(
            svg,
            r##"<line x1="{margin}" y1="{y:.1}" x2="{}" y2="{y:.1}" stroke="#ddd" stroke-width="0.5"/>"##,
            margin + plot_w,
        )
        .unwrap();
        writeln!(
            svg,
            r#"<text x="{}" y="{y:.1}" text-anchor="end" font-size="9" font-family="sans-serif" dy="3">{v_tick}</text>"#,
            margin - 5.0,
        )
        .unwrap();
    }

    // Detection limit line (LSST ~24.5)
    let det_y = v_to_y(24.5);
    writeln!(
        svg,
        r#"<line x1="{margin}" y1="{det_y:.1}" x2="{}" y2="{det_y:.1}" stroke="red" stroke-width="1" stroke-dasharray="6,3"/>"#,
        margin + plot_w,
    )
    .unwrap();
    writeln!(
        svg,
        r#"<text x="{}" y="{}" font-size="9" font-family="sans-serif" fill="red">LSST limit ~24.5</text>"#,
        margin + plot_w - 90.0,
        det_y - 5.0,
    )
    .unwrap();

    // Data points
    let max_points = population.len().min(2000);
    let step = if population.len() > max_points {
        population.len() / max_points
    } else {
        1
    };

    for obj in population.iter().step_by(step) {
        let x = a_to_x(obj.a);
        let y = v_to_y(obj.v_magnitude);

        if x >= margin && x <= margin + plot_w && y >= margin && y <= margin + plot_h {
            writeln!(
                svg,
                r#"<circle cx="{x:.1}" cy="{y:.1}" r="2" fill="steelblue" fill-opacity="0.4"/>"#,
            )
            .unwrap();
        }
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
        r#"<text x="15" y="{}" text-anchor="middle" font-size="12" font-family="sans-serif" transform="rotate(-90, 15, {})">V magnitude (fainter down)</text>"#,
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
    use crate::posterior::mcmc_2021_posterior;
    use crate::reference_population::generate_reference_population;
    use rand::SeedableRng;

    #[test]
    fn test_posterior_plot_valid_svg() {
        let posterior = mcmc_2021_posterior();
        let svg = posterior_plot(&posterior, 600, 500);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Planet Nine Posterior"));
        assert!(svg.contains("Semi-major axis"));
        assert!(svg.contains("Eccentricity"));
    }

    #[test]
    fn test_reference_population_plot_valid_svg() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let pop = generate_reference_population(200, &mut rng);
        let svg = reference_population_plot(&pop, 600, 500);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("</svg>"));
        assert!(svg.contains("Reference Population"));
        assert!(svg.contains("LSST"));
        assert!(svg.contains("circle"));
    }
}

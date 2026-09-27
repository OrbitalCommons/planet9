//! A small SVG writer for the strategy figures (dark theme, opaque backdrop).

use std::fmt::Write;

pub const BG: &str = "#1a1b26";
pub const PANEL: &str = "#16171f";
pub const FG: &str = "#c0caf5";
pub const MUTED: &str = "#565f89";
pub const BLUE: &str = "#7aa2f7";
pub const GREEN: &str = "#9ece6a";
pub const RED: &str = "#f7768e";
pub const ORANGE: &str = "#e0af68";
pub const PURPLE: &str = "#bb9af7";
pub const TEAL: &str = "#7dcfff";

/// Colours of the exposure tiers, shallow to deep.
pub const TIER_COLOURS: [&str; 7] = [
    "#7dcfff", "#7aa2f7", "#9ece6a", "#e0af68", "#ff9e64", "#f7768e", "#bb9af7",
];

pub struct Svg {
    body: String,
    pub width: f64,
    pub height: f64,
}

impl Svg {
    pub fn new(width: f64, height: f64) -> Self {
        let mut body = String::new();
        write!(
            body,
            r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" font-family="Helvetica, Arial, sans-serif"><rect width="{width}" height="{height}" fill="{BG}"/>"#
        )
        .unwrap();
        Self {
            body,
            width,
            height,
        }
    }

    pub fn rect(&mut self, x: f64, y: f64, w: f64, h: f64, fill: &str, opacity: f64) {
        write!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="{fill}" fill-opacity="{opacity:.3}"/>"#
        )
        .unwrap();
    }

    pub fn outline(&mut self, x: f64, y: f64, w: f64, h: f64, stroke: &str, width: f64) {
        write!(
            self.body,
            r#"<rect x="{x:.2}" y="{y:.2}" width="{w:.2}" height="{h:.2}" fill="none" stroke="{stroke}" stroke-width="{width}"/>"#
        )
        .unwrap();
    }

    #[allow(clippy::too_many_arguments)]
    pub fn line(
        &mut self,
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        stroke: &str,
        width: f64,
        dash: &str,
    ) {
        write!(
            self.body,
            r#"<line x1="{x1:.2}" y1="{y1:.2}" x2="{x2:.2}" y2="{y2:.2}" stroke="{stroke}" stroke-width="{width}" stroke-dasharray="{dash}"/>"#
        )
        .unwrap();
    }

    /// A polyline through `pts`; consecutive points further apart than
    /// `break_px` in x start a new segment (RA wrap).
    pub fn path(
        &mut self,
        pts: &[(f64, f64)],
        stroke: &str,
        width: f64,
        dash: &str,
        break_px: f64,
    ) {
        let mut d = String::new();
        let mut last: Option<(f64, f64)> = None;
        for &(x, y) in pts {
            let pen = match last {
                Some((lx, _)) if (x - lx).abs() <= break_px => 'L',
                _ => 'M',
            };
            write!(d, "{pen}{x:.2},{y:.2} ").unwrap();
            last = Some((x, y));
        }
        write!(
            self.body,
            r#"<path d="{d}" fill="none" stroke="{stroke}" stroke-width="{width}" stroke-dasharray="{dash}" stroke-linejoin="round"/>"#
        )
        .unwrap();
    }

    pub fn circle(&mut self, x: f64, y: f64, r: f64, fill: &str) {
        write!(
            self.body,
            r#"<circle cx="{x:.2}" cy="{y:.2}" r="{r}" fill="{fill}"/>"#
        )
        .unwrap();
    }

    /// Text anchored `start`, `middle` or `end`.
    pub fn text(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str, anchor: &str) {
        let s = s
            .replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;");
        write!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-size="{size}" fill="{fill}" text-anchor="{anchor}">{s}</text>"#
        )
        .unwrap();
    }

    pub fn bold(&mut self, x: f64, y: f64, s: &str, size: f64, fill: &str, anchor: &str) {
        let s = s
            .replace('&', "&amp;")
            .replace('<', "&lt;")
            .replace('>', "&gt;");
        write!(
            self.body,
            r#"<text x="{x:.2}" y="{y:.2}" font-size="{size}" font-weight="bold" fill="{fill}" text-anchor="{anchor}">{s}</text>"#
        )
        .unwrap();
    }

    pub fn finish(mut self) -> String {
        self.body.push_str("</svg>\n");
        self.body
    }
}

/// A linear data → pixel axis.
#[derive(Debug, Clone, Copy)]
pub struct Axis {
    pub d0: f64,
    pub d1: f64,
    pub p0: f64,
    pub p1: f64,
}

impl Axis {
    pub fn at(&self, d: f64) -> f64 {
        self.p0 + (d - self.d0) / (self.d1 - self.d0) * (self.p1 - self.p0)
    }
}

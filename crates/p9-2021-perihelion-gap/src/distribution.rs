//! Perihelion histogramming and the gap (dip) statistic.
//!
//! Oldroyd & Trujillo (2021) quantify the perihelion gap as a relative
//! *deficit* of objects at intermediate q (≈ 50–65 AU) compared with the two
//! flanking populations (the low-q scattering ETNOs and the high-q detached
//! IOCs). Here we bin a perihelion sample and compute a dip statistic: the
//! ratio of the count *density* (objects per AU) inside the gap window to the
//! mean density of the two flanks. A value < 1 is a deficit; the smaller, the
//! deeper the gap.

/// A perihelion histogram on a uniform grid.
#[derive(Debug, Clone)]
pub struct Histogram {
    pub edges: Vec<f64>,
    pub counts: Vec<usize>,
}

impl Histogram {
    /// Bin `values` into `nbins` uniform bins spanning [`lo`, `hi`].
    /// Values outside the range are dropped (the sample is pre-restricted to
    /// the distant range by the caller).
    pub fn new(values: &[f64], lo: f64, hi: f64, nbins: usize) -> Self {
        assert!(nbins >= 1 && hi > lo);
        let width = (hi - lo) / nbins as f64;
        let mut counts = vec![0usize; nbins];
        for &v in values {
            if v < lo || v >= hi {
                continue;
            }
            let mut idx = ((v - lo) / width) as usize;
            if idx >= nbins {
                idx = nbins - 1;
            }
            counts[idx] += 1;
        }
        let edges = (0..=nbins).map(|i| lo + i as f64 * width).collect();
        Self { edges, counts }
    }

    pub fn bin_width(&self) -> f64 {
        self.edges[1] - self.edges[0]
    }

    /// Total objects binned.
    pub fn total(&self) -> usize {
        self.counts.iter().sum()
    }
}

/// Count of values with q in the half-open window [`lo`, `hi`).
pub fn count_in_window(values: &[f64], lo: f64, hi: f64) -> usize {
    values.iter().filter(|&&q| q >= lo && q < hi).count()
}

/// Count density (objects per AU) in window [`lo`, `hi`).
pub fn density_in_window(values: &[f64], lo: f64, hi: f64) -> f64 {
    assert!(hi > lo);
    count_in_window(values, lo, hi) as f64 / (hi - lo)
}

/// The dip statistic for a perihelion sample.
#[derive(Debug, Clone, Copy)]
pub struct DipStatistic {
    /// Density (per AU) inside the gap window.
    pub gap_density: f64,
    /// Density (per AU) in the low-q (scattering) flank.
    pub low_flank_density: f64,
    /// Density (per AU) in the high-q (detached) flank.
    pub high_flank_density: f64,
    /// Mean of the two flank densities.
    pub flank_density: f64,
    /// gap_density / flank_density (< 1 ⇒ deficit).
    pub dip_ratio: f64,
}

/// Compute the dip statistic given a gap window [`gap_lo`, `gap_hi`) and flank
/// windows abutting it on each side, each `flank_width` AU wide.
///
/// The low flank is [gap_lo − flank_width, gap_lo); the high flank is
/// [gap_hi, gap_hi + flank_width). Comparing *densities* (per AU) makes the
/// statistic insensitive to the differing widths of the gap and flanks.
pub fn dip_statistic(values: &[f64], gap_lo: f64, gap_hi: f64, flank_width: f64) -> DipStatistic {
    let gap_density = density_in_window(values, gap_lo, gap_hi);
    let low_flank_density = density_in_window(values, gap_lo - flank_width, gap_lo);
    let high_flank_density = density_in_window(values, gap_hi, gap_hi + flank_width);
    let flank_density = 0.5 * (low_flank_density + high_flank_density);
    let dip_ratio = if flank_density > 0.0 {
        gap_density / flank_density
    } else {
        f64::INFINITY
    };
    DipStatistic {
        gap_density,
        low_flank_density,
        high_flank_density,
        flank_density,
        dip_ratio,
    }
}

/// Locate the deepest deficit: the interior bin with the lowest count that is
/// a genuine *valley* between two populations — there is a strictly
/// higher-count bin somewhere to its left AND somewhere to its right. Returns
/// the q at the centre of that minimum bin (ties broken toward the left).
///
/// Using "higher somewhere on each side" rather than "both immediate
/// neighbours higher" makes the locator robust to a noisy, multi-bin-wide
/// deficit (a real bimodal q distribution rarely has a single empty bin with
/// strictly taller immediate neighbours on both sides). Returns `None` if no
/// such valley exists (e.g. a flat or monotone distribution).
pub fn locate_gap_center(hist: &Histogram) -> Option<f64> {
    let n = hist.counts.len();
    if n < 3 {
        return None;
    }
    let mut best: Option<(usize, usize)> = None; // (count, idx)
    for i in 1..n - 1 {
        let c = hist.counts[i];
        let higher_left = hist.counts[..i].iter().any(|&x| x > c);
        let higher_right = hist.counts[i + 1..].iter().any(|&x| x > c);
        if higher_left && higher_right {
            match best {
                Some((bc, _)) if bc <= c => {}
                _ => best = Some((c, i)),
            }
        }
    }
    best.map(|(_, i)| 0.5 * (hist.edges[i] + hist.edges[i + 1]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn histogram_bins_and_totals() {
        let v = vec![41.0, 45.0, 48.0, 70.0, 76.0, 81.0];
        let h = Histogram::new(&v, 40.0, 90.0, 5); // 10-AU bins
        assert_eq!(h.total(), 6);
        assert_relative_eq!(h.bin_width(), 10.0);
        // bin0 [40,50): 41,45,48 -> 3; bin1 [50,60): 0; bin3 [70,80): 70,76 -> 2
        assert_eq!(h.counts[0], 3);
        assert_eq!(h.counts[1], 0);
        assert_eq!(h.counts[3], 2);
    }

    #[test]
    fn density_is_per_au() {
        let v = vec![52.0, 55.0, 60.0];
        // 3 objects over a 15-AU window -> 0.2 per AU.
        assert_relative_eq!(density_in_window(&v, 50.0, 65.0), 3.0 / 15.0);
    }

    #[test]
    fn dip_ratio_below_one_for_a_deficit() {
        // Dense flanks, empty gap.
        let mut v = vec![];
        v.extend([40.0, 42.0, 44.0, 46.0, 48.0]); // low flank
        v.extend([66.0, 68.0, 70.0, 72.0, 74.0]); // high flank
        v.push(57.0); // one lonely gap object
        let d = dip_statistic(&v, 50.0, 65.0, 10.0);
        assert!(d.dip_ratio < 1.0, "dip_ratio = {}", d.dip_ratio);
        assert!(d.gap_density < d.flank_density);
    }

    #[test]
    fn locate_gap_finds_the_empty_interior_bin() {
        let v = vec![41.0, 45.0, 48.0, 67.0, 70.0, 76.0];
        let h = Histogram::new(&v, 40.0, 90.0, 5); // bins centred 45,55,65,75,85
        let c = locate_gap_center(&h).expect("a gap exists");
        // The [50,60) bin is empty and flanked -> centre 55.
        assert_relative_eq!(c, 55.0);
    }

    // ---- Headline reproduction: the gap in the OBSERVED distant-TNO sample ----

    #[test]
    fn observed_sample_has_a_perihelion_deficit_at_50_to_65() {
        use crate::published::{GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
        use crate::sample::observed_perihelia;

        let qs = observed_perihelia();
        let d = dip_statistic(&qs, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, 12.0);
        // The real, vetted sample shows a strong deficit at q = 50-65 AU:
        // the count density there is far below both flanks.
        assert!(
            d.dip_ratio < 0.3,
            "observed sample dip_ratio = {} (expected a deep deficit)",
            d.dip_ratio
        );
        assert!(d.gap_density < d.low_flank_density);
        assert!(d.gap_density < d.high_flank_density);
    }

    #[test]
    fn observed_gap_center_matches_published_window() {
        use crate::published::{GAP_Q_CENTER_AU, GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
        use crate::sample::observed_perihelia;

        let qs = observed_perihelia();
        // 13 bins of 5 AU over [35, 100). Bin centres: 37.5, 42.5, ... 97.5.
        let h = Histogram::new(&qs, 35.0, 100.0, 13);
        let center = locate_gap_center(&h).expect("observed gap present");
        // The located minimum-density interior bin sits inside the published
        // 50-65 AU window, within a stated tolerance of its centre (57.5 AU).
        assert!(
            center >= GAP_Q_LOW_AU && center <= GAP_Q_HIGH_AU,
            "observed gap center {center} outside [{GAP_Q_LOW_AU}, {GAP_Q_HIGH_AU}]"
        );
        assert!(
            (center - GAP_Q_CENTER_AU).abs() <= 5.0,
            "observed gap center {center} vs published {GAP_Q_CENTER_AU}"
        );
    }

    #[test]
    fn high_e_subset_also_shows_the_gap() {
        use crate::published::{GAP_ECCENTRICITY_FLOOR, GAP_Q_HIGH_AU, GAP_Q_LOW_AU};
        use crate::sample::high_e_perihelia;

        // The paper notes the gap is clearest at e >= 0.65.
        let qs = high_e_perihelia(GAP_ECCENTRICITY_FLOOR);
        let d = dip_statistic(&qs, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, 12.0);
        assert!(
            d.dip_ratio < 0.5,
            "high-e subset dip_ratio = {}",
            d.dip_ratio
        );
    }
}

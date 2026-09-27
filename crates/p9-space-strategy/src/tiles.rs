//! An equal-area tiling of the whole sky: 2° declination bands, each cut
//! into right-ascension cells whose count scales with cos(Dec), so every tile
//! is ~4 deg² from the equator to the poles. Tiles are the planning unit:
//! probability is accumulated per tile, exposure is allocated per tile.

use p9_core::constants::{DEG2RAD, RAD2DEG};

/// Declination height of every band (deg).
pub const BAND_DEG: f64 = 2.0;

#[derive(Debug, Clone)]
struct Band {
    dec_lo: f64,
    n_ra: usize,
    first: usize,
}

/// The tiling.
#[derive(Debug, Clone)]
pub struct Tiles {
    bands: Vec<Band>,
    n: usize,
}

impl Default for Tiles {
    fn default() -> Self {
        Self::new()
    }
}

impl Tiles {
    pub fn new() -> Self {
        let n_bands = (180.0 / BAND_DEG).round() as usize;
        let mut bands = Vec::with_capacity(n_bands);
        let mut first = 0;
        for k in 0..n_bands {
            let dec_lo = -90.0 + k as f64 * BAND_DEG;
            let dec_c = dec_lo + 0.5 * BAND_DEG;
            let n_ra = ((360.0 * (dec_c * DEG2RAD).cos() / BAND_DEG).round() as usize).max(1);
            bands.push(Band {
                dec_lo,
                n_ra,
                first,
            });
            first += n_ra;
        }
        Self { bands, n: first }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    fn band_of(&self, dec_deg: f64) -> &Band {
        let k = (((dec_deg + 90.0) / BAND_DEG).floor() as isize)
            .clamp(0, self.bands.len() as isize - 1) as usize;
        &self.bands[k]
    }

    fn band_containing(&self, idx: usize) -> &Band {
        let k = self.bands.partition_point(|b| b.first <= idx) - 1;
        &self.bands[k]
    }

    /// Tile index of a direction (degrees).
    pub fn index(&self, ra_deg: f64, dec_deg: f64) -> usize {
        let b = self.band_of(dec_deg);
        let w = 360.0 / b.n_ra as f64;
        let k = ((ra_deg.rem_euclid(360.0) / w).floor() as usize).min(b.n_ra - 1);
        b.first + k
    }

    /// Centre (RA, Dec) of a tile in degrees.
    pub fn centre(&self, idx: usize) -> (f64, f64) {
        let b = self.band_containing(idx);
        let w = 360.0 / b.n_ra as f64;
        (
            (idx - b.first) as f64 * w + 0.5 * w,
            b.dec_lo + 0.5 * BAND_DEG,
        )
    }

    /// (RA width, Dec height) of a tile in degrees.
    pub fn extent(&self, idx: usize) -> (f64, f64) {
        let b = self.band_containing(idx);
        (360.0 / b.n_ra as f64, BAND_DEG)
    }

    /// Solid angle of a tile (deg²).
    pub fn area_deg2(&self, idx: usize) -> f64 {
        let b = self.band_containing(idx);
        let (lo, hi) = (b.dec_lo * DEG2RAD, (b.dec_lo + BAND_DEG) * DEG2RAD);
        360.0 * RAD2DEG * (hi.sin() - lo.sin()) / b.n_ra as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tiles_cover_the_sphere_with_near_equal_area() {
        let t = Tiles::new();
        let total: f64 = (0..t.len()).map(|k| t.area_deg2(k)).sum();
        assert!((total - 41_252.96).abs() < 1.0, "total = {total}");
        for k in 0..t.len() {
            let (_, dec) = t.centre(k);
            if dec.abs() < 80.0 {
                let a = t.area_deg2(k);
                assert!((3.6..4.4).contains(&a), "tile {k} at dec {dec}: {a} deg2");
            }
        }
    }

    #[test]
    fn index_and_centre_round_trip() {
        let t = Tiles::new();
        for k in (0..t.len()).step_by(37) {
            let (ra, dec) = t.centre(k);
            assert_eq!(t.index(ra, dec), k);
        }
        assert_eq!(t.index(359.999, 0.5), t.index(-0.001, 0.5));
    }
}

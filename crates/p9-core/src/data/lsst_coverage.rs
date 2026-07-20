//! Loader for the rubin-watch LSST coverage map
//! (`rubin_watch/lsst_coverage.json`, schema v1 sparse encoding per
//! `rubin_watch/design/07-schemas.md`).
//!
//! Pure parsing/validation — consumers (p9-search-hull, p9-viability,
//! p9-rubin-watch) decide what depth model to apply. Pixel indices are
//! HEALPix RING at the file's nside, produced by healpy in the builder and
//! matched by [`crate::coords::healpix::ang2pix_ring`].

use std::collections::BTreeMap;

pub const SCHEMA_VERSION: u64 = 1;

#[derive(Debug, Clone)]
pub struct CoverageMap {
    pub nside: u32,
    pub generated_utc: String,
    pub first_night: String,
    pub last_night: String,
    pub source: String,
    pub n_reconstructed_visits: u64,
    /// Band → single-visit 5σ fiducial depth (mag).
    pub fiducial_depth: BTreeMap<String, f64>,
    /// Band → sparse per-pixel coverage (aligned, pixel-sorted).
    pub bands: BTreeMap<String, BandCoverage>,
    /// Pixels flagged |b| < 10° (crowding-limited DIA).
    pub crowding_pixels: Vec<usize>,
    /// Pixels whose templates include survey-year data (year-2 hazard).
    pub template_epoch_pixels: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct BandCoverage {
    pub pixels: Vec<usize>,
    pub n_visits: Vec<u32>,
    pub last_visit_mjd: Vec<f64>,
}

impl BandCoverage {
    /// Visits at a pixel (0 = never observed). Pixels are sorted, so this is
    /// a binary search.
    pub fn n_visits_at(&self, pix: usize) -> u32 {
        match self.pixels.binary_search(&pix) {
            Ok(i) => self.n_visits[i],
            Err(_) => 0,
        }
    }
}

impl CoverageMap {
    pub fn parse(json: &str) -> Result<Self, String> {
        let v: serde_json::Value =
            serde_json::from_str(json).map_err(|e| format!("coverage json: {e}"))?;
        let sv = v["schema_version"].as_u64().ok_or("no schema_version")?;
        if sv != SCHEMA_VERSION {
            return Err(format!("coverage schema_version {sv} != {SCHEMA_VERSION}"));
        }
        let ordering = v["healpix"]["ordering"].as_str().ok_or("no ordering")?;
        if ordering != "ring" {
            return Err(format!("unsupported healpix ordering {ordering:?}"));
        }
        let nside = v["healpix"]["nside"].as_u64().ok_or("no nside")? as u32;
        if !nside.is_power_of_two() {
            return Err(format!("nside {nside} not a power of two"));
        }
        let npix = crate::coords::healpix::npix(nside);

        let fiducial_depth: BTreeMap<String, f64> = v["band_fiducial_depth"]
            .as_object()
            .ok_or("no band_fiducial_depth")?
            .iter()
            .filter_map(|(k, x)| x.as_f64().map(|d| (k.clone(), d)))
            .collect();

        let mut bands = BTreeMap::new();
        for (band, b) in v["bands"].as_object().ok_or("no bands")? {
            let usize_list = |key: &str| -> Result<Vec<usize>, String> {
                b[key]
                    .as_array()
                    .ok_or_else(|| format!("band {band}: no {key}"))?
                    .iter()
                    .map(|x| {
                        x.as_u64()
                            .map(|u| u as usize)
                            .ok_or_else(|| format!("{band}.{key}"))
                    })
                    .collect()
            };
            let pixels = usize_list("pixels")?;
            let n_visits: Vec<u32> = usize_list("n_visits")?
                .into_iter()
                .map(|u| u as u32)
                .collect();
            let last_visit_mjd: Vec<f64> = b["last_visit_mjd"]
                .as_array()
                .ok_or_else(|| format!("band {band}: no last_visit_mjd"))?
                .iter()
                .map(|x| x.as_f64().ok_or_else(|| format!("{band}.mjd")))
                .collect::<Result<_, _>>()?;
            if pixels.len() != n_visits.len() || pixels.len() != last_visit_mjd.len() {
                return Err(format!("band {band}: unaligned sparse arrays"));
            }
            if !pixels.windows(2).all(|w| w[0] < w[1]) {
                return Err(format!("band {band}: pixels not strictly sorted"));
            }
            if pixels.last().is_some_and(|&p| p >= npix) {
                return Err(format!("band {band}: pixel index out of range"));
            }
            if n_visits.iter().any(|&n| n == 0) {
                return Err(format!("band {band}: zero-visit pixel in sparse list"));
            }
            if !fiducial_depth.contains_key(band) {
                return Err(format!("band {band}: no fiducial depth"));
            }
            bands.insert(
                band.clone(),
                BandCoverage {
                    pixels,
                    n_visits,
                    last_visit_mjd,
                },
            );
        }

        let flag_list = |key: &str| -> Vec<usize> {
            v["flags"][key]
                .as_array()
                .map(|a| {
                    a.iter()
                        .filter_map(|x| x.as_u64().map(|u| u as usize))
                        .collect()
                })
                .unwrap_or_default()
        };

        Ok(CoverageMap {
            nside,
            generated_utc: v["generated_utc"].as_str().unwrap_or("").to_string(),
            first_night: v["window"]["first_night"]
                .as_str()
                .unwrap_or("")
                .to_string(),
            last_night: v["window"]["last_night"].as_str().unwrap_or("").to_string(),
            source: v["source"].as_str().unwrap_or("").to_string(),
            n_reconstructed_visits: v["n_reconstructed_visits"].as_u64().unwrap_or(0),
            fiducial_depth,
            bands,
            crowding_pixels: flag_list("crowding_pixels"),
            template_epoch_pixels: flag_list("template_epoch_pixels"),
        })
    }

    pub fn load(path: &std::path::Path) -> Result<Self, String> {
        let raw = std::fs::read_to_string(path).map_err(|e| format!("read {path:?}: {e}"))?;
        Self::parse(&raw)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny() -> String {
        r#"{
          "schema_version": 1,
          "generated_utc": "t",
          "healpix": {"nside": 64, "ordering": "ring"},
          "window": {"first_night": "2026-07-13", "last_night": "2026-07-14"},
          "source": "test",
          "n_reconstructed_visits": 3,
          "band_fiducial_depth": {"r": 24.3},
          "bands": {"r": {"pixels": [10, 20, 30], "n_visits": [1, 3, 5],
                          "last_visit_mjd": [61234.1, 61234.2, 61234.3]}},
          "flags": {"template_epoch_pixels": [], "crowding_pixels": [20]}
        }"#
        .to_string()
    }

    #[test]
    fn parses_and_queries() {
        let m = CoverageMap::parse(&tiny()).unwrap();
        assert_eq!(m.nside, 64);
        let r = &m.bands["r"];
        assert_eq!(r.n_visits_at(20), 3);
        assert_eq!(r.n_visits_at(21), 0);
        assert_eq!(m.crowding_pixels, vec![20]);
        assert_eq!(m.fiducial_depth["r"], 24.3);
    }

    #[test]
    fn rejects_malformed() {
        for (bad, why) in [
            (
                tiny().replace("\"schema_version\": 1", "\"schema_version\": 2"),
                "version",
            ),
            (tiny().replace("\"ring\"", "\"nested\""), "ordering"),
            (tiny().replace("[10, 20, 30]", "[30, 20, 10]"), "sorted"),
            (tiny().replace("[1, 3, 5]", "[1, 3]"), "aligned"),
            (tiny().replace("[1, 3, 5]", "[1, 0, 5]"), "zero-visit"),
            (tiny().replace("[10, 20, 30]", "[10, 20, 49152]"), "range"),
        ] {
            assert!(CoverageMap::parse(&bad).is_err(), "{why}");
        }
    }
}

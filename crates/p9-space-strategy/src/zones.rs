//! The four parts of the sky the plan divides into, named by why a space
//! telescope is needed there.

use serde::Serialize;

use crate::field::Tile;
use crate::prior::{RUBIN_NES_BETA_MAX_DEG, RUBIN_NES_DEC_MAX_DEG, RUBIN_WFD_DEC_MAX_DEG};

/// |b| inside which crowding defeats the ground surveys (deg).
pub const PLANE_BAND_DEG: f64 = 15.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
pub enum Zone {
    /// Where the predicted orbit crosses the Milky Way opposite the Galactic
    /// centre (Taurus–Gemini–Orion, RA 4–8 h), on the aphelion side.
    AnticentreCrossing,
    /// Where it crosses toward the Galactic centre (Ophiuchus–Sagittarius,
    /// RA 16–20 h), on the perihelion side.
    CentreCrossing,
    /// Clean sky north of Rubin's reach.
    NorthOfRubin,
    /// Clean sky Rubin also covers; only what is too faint for Rubin is left.
    RubinOverlap,
}

impl Zone {
    pub const ALL: [Zone; 4] = [
        Zone::AnticentreCrossing,
        Zone::CentreCrossing,
        Zone::NorthOfRubin,
        Zone::RubinOverlap,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Zone::AnticentreCrossing => "Anticentre crossing",
            Zone::CentreCrossing => "Galactic-centre crossing",
            Zone::NorthOfRubin => "North of Rubin",
            Zone::RubinOverlap => "Rubin overlap",
        }
    }

    pub fn why(self) -> &'static str {
        match self {
            Zone::AnticentreCrossing => {
                "ground surveys are crowded out; the planet dwells here, far and faint"
            }
            Zone::CentreCrossing => {
                "ground surveys and Rubin are crowded out; the planet would be near and bright"
            }
            Zone::NorthOfRubin => "beyond Rubin's declination limit and off its ecliptic spur",
            Zone::RubinOverlap => "only draws too faint for Rubin remain",
        }
    }
}

/// Is a direction inside Rubin's footprint (wide-fast-deep or ecliptic spur)?
pub fn in_rubin_footprint(dec_deg: f64, ecl_lat_deg: f64) -> bool {
    dec_deg <= RUBIN_WFD_DEC_MAX_DEG
        || (dec_deg <= RUBIN_NES_DEC_MAX_DEG && ecl_lat_deg.abs() <= RUBIN_NES_BETA_MAX_DEG)
}

pub fn zone_of(tile: &Tile) -> Zone {
    if tile.gal_b_deg.abs() < PLANE_BAND_DEG {
        // Galactic longitude within 90° of the centre, or of the anticentre.
        if tile.gal_l_deg.to_radians().cos() > 0.0 {
            Zone::CentreCrossing
        } else {
            Zone::AnticentreCrossing
        }
    } else if in_rubin_footprint(tile.dec_deg, tile.ecl_lat_deg) {
        Zone::RubinOverlap
    } else {
        Zone::NorthOfRubin
    }
}

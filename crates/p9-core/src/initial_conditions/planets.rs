//! Giant planet initial conditions.
//!
//! Provides state vectors for Jupiter, Saturn, Uranus, and Neptune
//! at the J2000.0 epoch in heliocentric ecliptic coordinates.
//!
//! These values are from JPL DE440 evaluated at J2000.0 (JD 2451545.0 TDB).
//! Units: AU for position, AU/day for velocity.
//!
//! TODO: When starfield-datasources#7 (DE440 support) lands, load these
//! directly from the ephemeris kernel for arbitrary epochs.

use nalgebra::Vector3;

use crate::constants::*;
use crate::types::{MassiveBody, StateVector};

/// Jupiter state vector at J2000.0, heliocentric ecliptic J2000.
pub fn jupiter_j2000() -> MassiveBody {
    MassiveBody {
        name: "Jupiter".to_string(),
        gm: GM_JUPITER,
        mass: MASS_JUPITER_SOLAR,
        state: StateVector::new(
            Vector3::new(
                3.996_321_311_680_50,
                2.932_561_367_455_14,
                -0.101_587_544_488_97,
            ),
            Vector3::new(
                -0.004_557_172_659_90,
                0.006_436_544_686_80,
                0.000_075_602_835_98,
            ),
        ),
        radius_au: RADIUS_JUPITER_AU,
        j2: Some(J2_JUPITER),
        j4: Some(J4_JUPITER),
    }
}

/// Saturn state vector at J2000.0, heliocentric ecliptic J2000.
pub fn saturn_j2000() -> MassiveBody {
    MassiveBody {
        name: "Saturn".to_string(),
        gm: GM_SATURN,
        mass: MASS_SATURN_SOLAR,
        state: StateVector::new(
            Vector3::new(
                6.401_416_557_621_40,
                6.565_273_850_498_26,
                -0.369_100_456_530_09,
            ),
            Vector3::new(
                -0.004_286_537_968_16,
                0.003_886_173_284_11,
                0.000_102_141_973_20,
            ),
        ),
        radius_au: RADIUS_SATURN_AU,
        j2: Some(J2_SATURN),
        j4: Some(J4_SATURN),
    }
}

/// Uranus state vector at J2000.0, heliocentric ecliptic J2000.
pub fn uranus_j2000() -> MassiveBody {
    MassiveBody {
        name: "Uranus".to_string(),
        gm: GM_URANUS,
        mass: MASS_URANUS_SOLAR,
        state: StateVector::new(
            Vector3::new(
                14.424_674_534_772_6,
                -13.737_303_893_986_8,
                -0.238_073_479_437_32,
            ),
            Vector3::new(
                0.002_681_326_018_12,
                0.002_673_344_674_32,
                -0.000_024_788_970_42,
            ),
        ),
        radius_au: RADIUS_URANUS_AU,
        j2: Some(J2_URANUS),
        j4: Some(J4_URANUS),
    }
}

/// Neptune state vector at J2000.0, heliocentric ecliptic J2000.
pub fn neptune_j2000() -> MassiveBody {
    MassiveBody {
        name: "Neptune".to_string(),
        gm: GM_NEPTUNE,
        mass: MASS_NEPTUNE_SOLAR,
        state: StateVector::new(
            Vector3::new(
                16.804_971_613_081_8,
                -24.990_517_792_794_2,
                0.127_392_774_785_86,
            ),
            Vector3::new(
                0.002_583_519_281_29,
                0.001_768_892_426_37,
                -0.000_095_159_233_52,
            ),
        ),
        radius_au: RADIUS_NEPTUNE_AU,
        j2: Some(J2_NEPTUNE),
        j4: Some(J4_NEPTUNE),
    }
}

/// All four giant planets at J2000.0.
pub fn giant_planets_j2000() -> Vec<MassiveBody> {
    vec![
        jupiter_j2000(),
        saturn_j2000(),
        uranus_j2000(),
        neptune_j2000(),
    ]
}

/// Neptune only (many papers use J2 for inner giants + direct Neptune).
pub fn neptune_only_j2000() -> Vec<MassiveBody> {
    vec![neptune_j2000()]
}

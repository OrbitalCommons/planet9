//! Dynamical groups of the paper's Table 1, defined on perihelion distance
//! `p` (AU), eccentricity and inclination, and the group census.
//!
//! Two conventions the Table leaves implicit are fixed here:
//!
//! * **Weights.** The flyby is run with a constant tracer surface density
//!   for resolution and post-processed "by assigning different masses to
//!   the particles to model the actual mass density distribution"
//!   (Section 2); Pfalzner et al. (2018) state the assumed profile,
//!   Σ ∝ 1/r. Every tracer therefore carries the weight
//!   [`surface_density_weight`] of its pre-flyby radius.
//! * **Denominator.** The Table's rows cover 0.1 < p < 100 AU, so `N` is
//!   the weight of bound tracers with p ≤ 100 AU; the unperturbed outer disc
//!   beyond that is not counted (the paper compares only the region
//!   "sufficiently covered by observations").
//!
//! The hot belt is read as the 30.1 < p < 48 AU, i < 35° orbits that are
//! not cold (e ≥ 0.15 or i ≥ 10°): the Table lists it as e > 0.15 with
//! 10° ≤ i ≤ 35°, which would leave eccentric low-inclination belt orbits —
//! a large fraction of the post-flyby belt — in no family at all.

use p9_core::constants::RAD2DEG;
use p9_core::types::OrbitalElements;

/// Post-processing weight of a tracer that formed at `r0_au` in a
/// constant-surface-density run, for a Σ ∝ 1/r disc (Pfalzner et al. 2018).
pub fn surface_density_weight(r0_au: f64) -> f64 {
    1.0 / r0_au
}

/// Outer perihelion edge of the region the Table covers (AU).
pub const CENSUS_P_MAX: f64 = 100.0;

/// Table 1 dynamical families.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum DynamicalGroup {
    /// p ∈ [30.1, 48], e < 0.15, i < 10°
    ColdKb,
    /// p ∈ [30.1, 48], i < 35°, not cold (e ≥ 0.15 or i ≥ 10°)
    HotKb,
    /// p ∈ [48, 60], e > 0.24, i < 35°
    Detached,
    /// p ∈ [60, 100], e > 0.6, i < 35°
    SednaLike,
    /// p ∈ [30.1, 60], 35° ≤ i ≤ 90°
    Inclined,
    /// p ∈ [30.1, 60], i > 90°
    Retrograde,
    /// p ∈ [0.1, 29.5]: injected into the planetary region
    Inner,
}

/// Table order.
pub const GROUPS: [DynamicalGroup; 7] = [
    DynamicalGroup::ColdKb,
    DynamicalGroup::HotKb,
    DynamicalGroup::Detached,
    DynamicalGroup::SednaLike,
    DynamicalGroup::Inclined,
    DynamicalGroup::Retrograde,
    DynamicalGroup::Inner,
];

impl DynamicalGroup {
    pub fn index(self) -> usize {
        GROUPS.iter().position(|&g| g == self).unwrap()
    }

    pub fn label(self) -> &'static str {
        match self {
            DynamicalGroup::ColdKb => "cold KB",
            DynamicalGroup::HotKb => "hot KB",
            DynamicalGroup::Detached => "detached",
            DynamicalGroup::SednaLike => "Sedna-like",
            DynamicalGroup::Inclined => "inclined",
            DynamicalGroup::Retrograde => "retrograde",
            DynamicalGroup::Inner => "inner",
        }
    }
}

/// Classify a bound orbit; `None` when it falls in none of the Table 1
/// boxes (e.g. p ∈ (29.5, 30.1), p > 100 AU, or a low-e orbit with
/// 10° ≤ i ≤ 35°).
pub fn classify(el: &OrbitalElements) -> Option<DynamicalGroup> {
    let p = el.a * (1.0 - el.e);
    let e = el.e;
    let i = el.i * RAD2DEG;
    if (0.1..=29.5).contains(&p) {
        return Some(DynamicalGroup::Inner);
    }
    if p < 30.1 {
        return None;
    }
    if p <= 60.0 && i > 90.0 {
        return Some(DynamicalGroup::Retrograde);
    }
    if p <= 60.0 && i >= 35.0 {
        return Some(DynamicalGroup::Inclined);
    }
    if p <= 48.0 {
        return Some(if e < 0.15 && i < 10.0 {
            DynamicalGroup::ColdKb
        } else {
            DynamicalGroup::HotKb
        });
    }
    if p <= 60.0 {
        return (e > 0.24).then_some(DynamicalGroup::Detached);
    }
    if p <= 100.0 {
        return (e > 0.6).then_some(DynamicalGroup::SednaLike);
    }
    None
}

/// Weighted group census over a population.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Census {
    /// Total weight in each group, Table order.
    pub weights: [f64; 7],
    /// Weight of bound orbits with p ≤ [`CENSUS_P_MAX`] (the Table's N).
    pub total: f64,
    /// Unweighted number of orbits behind `total`.
    pub n: usize,
}

impl Census {
    /// Census of `(pre-flyby radius, elements)` pairs with Σ ∝ 1/r weights.
    pub fn of<'a>(bodies: impl IntoIterator<Item = (f64, &'a OrbitalElements)>) -> Self {
        let mut weights = [0.0; 7];
        let mut total = 0.0;
        let mut n = 0;
        for (r0, el) in bodies {
            if el.a * (1.0 - el.e) > CENSUS_P_MAX {
                continue;
            }
            let w = surface_density_weight(r0);
            total += w;
            n += 1;
            if let Some(g) = classify(el) {
                weights[g.index()] += w;
            }
        }
        Self { weights, total, n }
    }

    pub fn weight(&self, g: DynamicalGroup) -> f64 {
        self.weights[g.index()]
    }

    /// Table 1 "N_x / N".
    pub fn fraction(&self, g: DynamicalGroup) -> f64 {
        if self.total <= 0.0 {
            0.0
        } else {
            self.weight(g) / self.total
        }
    }

    pub fn fractions(&self) -> [f64; 7] {
        let mut out = [0.0; 7];
        for (k, g) in GROUPS.iter().enumerate() {
            out[k] = self.fraction(*g);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    fn el(p: f64, e: f64, i_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: p / (1.0 - e),
            e,
            i: i_deg * DEG2RAD,
            omega_big: 0.0,
            omega: 0.0,
            mean_anomaly: 0.0,
        }
    }

    #[test]
    fn table_boxes() {
        assert_eq!(classify(&el(40.0, 0.05, 3.0)), Some(DynamicalGroup::ColdKb));
        assert_eq!(classify(&el(40.0, 0.4, 20.0)), Some(DynamicalGroup::HotKb));
        assert_eq!(
            classify(&el(55.0, 0.5, 20.0)),
            Some(DynamicalGroup::Detached)
        );
        assert_eq!(
            classify(&el(76.0, 0.85, 12.0)),
            Some(DynamicalGroup::SednaLike)
        );
        assert_eq!(
            classify(&el(45.0, 0.5, 60.0)),
            Some(DynamicalGroup::Inclined)
        );
        assert_eq!(
            classify(&el(35.0, 0.6, 110.0)),
            Some(DynamicalGroup::Retrograde)
        );
        assert_eq!(classify(&el(20.0, 0.3, 5.0)), Some(DynamicalGroup::Inner));
        assert_eq!(classify(&el(40.0, 0.05, 20.0)), Some(DynamicalGroup::HotKb));
        assert_eq!(classify(&el(40.0, 0.4, 5.0)), Some(DynamicalGroup::HotKb));
        assert_eq!(classify(&el(55.0, 0.1, 5.0)), None);
        assert_eq!(classify(&el(120.0, 0.9, 5.0)), None);
    }
}

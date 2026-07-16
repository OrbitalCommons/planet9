//! The massive, eccentric, apsidally-aligned debris disk as a set of confocal
//! eccentric rings (Sefilian & Touma 2019, arXiv:1804.06859).
//!
//! The disk spans the paper's 40-750 AU region, carries a total mass M_disk of
//! order a few-to-ten Earth masses, and is internally apsidally aligned: every ring
//! shares the disk's longitude of perihelion `varpi_disk` and follows a
//! prescribed eccentricity profile e_d(a). Such a disk is what the paper shows
//! can SHEPHERD test ETNOs into apsidal anti-alignment without a planet.

use p9_core::units::{au, radians, solar_masses, Angle, Length, Mass};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Reference / published values from Sefilian & Touma (2019), kept as labelled
/// constants (never silently baked into the computation).
pub mod reference {
    /// Inner edge of the modeled disk (AU): the paper's a_in = 40 AU
    /// (Sefilian & Touma 2019 Sec. 2; a previous version carried 100 AU
    /// under this "published values" label — false provenance).
    pub const A_INNER_AU: f64 = 40.0;
    /// Outer edge of the modeled disk (AU): the paper's a_out = 750 AU
    /// (previously 1000 AU).
    pub const A_OUTER_AU: f64 = 750.0;
    /// Fiducial disk mass: the paper explores a few to ~10 Earth masses; its
    /// headline shepherding case uses ~10 M_earth (Sec. 3-4).
    pub const M_DISK_EARTH_FIDUCIAL: f64 = 10.0;
    /// Disk longitude of perihelion in the paper's frame (deg). The test
    /// particles confine ANTI-aligned, i.e. near varpi_disk + 180 deg.
    pub const VARPI_DISK_DEG: f64 = 0.0;
    /// Secular precession period induced by a ~10 M_earth disk on a ~250 AU
    /// test particle is of order 1e8-1e9 yr (paper Sec. 3, their Fig. 2-3
    /// libration periods). Documented order of magnitude, not a point value.
    pub const PRECESSION_PERIOD_LO_YR: f64 = 1.0e8;
    pub const PRECESSION_PERIOD_HI_YR: f64 = 1.0e9;
}

/// Power-law indices for the disk surface-density and eccentricity profiles.
/// Sefilian & Touma use Sigma ~ a^{-p} with an apsidally-aligned eccentricity
/// e_d(a). Defaults give a mildly declining surface density and a moderate,
/// roughly flat disk eccentricity, the regime where shepherding is strongest.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DiskProfile {
    /// Inner edge (AU).
    pub a_in: f64,
    /// Outer edge (AU).
    pub a_out: f64,
    /// Total disk mass in solar masses.
    pub mass_solar: f64,
    /// Surface-density power-law index p in Sigma ~ a^{-p}.
    pub sigma_index: f64,
    /// Disk eccentricity at the reference radius.
    pub e_disk: f64,
    /// Disk longitude of perihelion (radians); all rings share it.
    pub varpi_disk: f64,
    /// Softening as a fraction of the local semi-major axis (finite ring
    /// width / disk thickness). Regularizes the alpha -> 1 Laplace divergence.
    pub softening_frac: f64,
    /// Number of rings discretizing the continuous disk.
    pub n_rings: usize,
}

impl DiskProfile {
    /// Fiducial Sefilian & Touma disk: 100-1000 AU, given mass in Earth masses,
    /// Sigma ~ a^{-2}, e_disk = 0.2, apse along varpi_disk.
    pub fn fiducial(mass_earth: f64) -> Self {
        Self {
            a_in: reference::A_INNER_AU,
            a_out: reference::A_OUTER_AU,
            mass_solar: mass_earth * p9_core::constants::EARTH_MASS_SOLAR,
            sigma_index: 2.0,
            e_disk: 0.2,
            varpi_disk: reference::VARPI_DISK_DEG * p9_core::constants::DEG2RAD,
            // Fractional softening (~ disk scale height / radial smearing). A
            // dynamically warm (e ~ 0.2) debris disk has a finite thickness of
            // tens of percent; 0.2 is the value at which a 10 M_earth disk
            // gives the paper's ~1e8-1e9 yr secular timescale at a = 250 AU
            // (see `secular` tests and REPRODUCTION_NOTES residual).
            softening_frac: 0.2,
            n_rings: 200,
        }
    }

    pub fn mass_earth(&self) -> f64 {
        self.mass_solar / p9_core::constants::EARTH_MASS_SOLAR
    }

    /// Inner edge as a typed [`Length`].
    pub fn inner_radius(&self) -> Length {
        au(self.a_in)
    }

    /// Outer edge as a typed [`Length`].
    pub fn outer_radius(&self) -> Length {
        au(self.a_out)
    }

    /// Total disk mass as a typed [`Mass`].
    pub fn mass(&self) -> Mass {
        solar_masses(self.mass_solar)
    }

    /// Disk longitude of perihelion as a typed [`Angle`].
    pub fn varpi_disk_angle(&self) -> Angle {
        radians(self.varpi_disk)
    }

    /// A single confocal eccentric ring of the discretized disk.
    fn ring(&self, k: usize) -> Ring {
        // Log-spaced ring semi-major axes (the disk spans a decade in a).
        let frac = (k as f64 + 0.5) / self.n_rings as f64;
        let ln_a = self.a_in.ln() + frac * (self.a_out.ln() - self.a_in.ln());
        let a = ln_a.exp();
        Ring {
            a,
            e: self.e_disk,
            varpi: self.varpi_disk,
            mass_solar: self.ring_mass(k),
        }
    }

    /// Mass of ring k from the Sigma ~ a^{-p} profile, normalized so the rings
    /// sum to `mass_solar`. With log-spaced a, the annulus mass is
    /// dM = 2 pi Sigma(a) a da, and a da = a^2 d(ln a) with d(ln a) uniform, so
    /// the unnormalized weight is a^{2 - p}.
    fn ring_mass(&self, k: usize) -> f64 {
        let frac = (k as f64 + 0.5) / self.n_rings as f64;
        let ln_a = self.a_in.ln() + frac * (self.a_out.ln() - self.a_in.ln());
        let a = ln_a.exp();
        let w = a.powf(2.0 - self.sigma_index);
        let total: f64 = (0..self.n_rings)
            .map(|m| {
                let f = (m as f64 + 0.5) / self.n_rings as f64;
                let la = self.a_in.ln() + f * (self.a_out.ln() - self.a_in.ln());
                la.exp().powf(2.0 - self.sigma_index)
            })
            .sum();
        self.mass_solar * w / total
    }

    /// All rings.
    pub fn rings(&self) -> Vec<Ring> {
        (0..self.n_rings).map(|k| self.ring(k)).collect()
    }
}

/// A confocal eccentric ring: orbit-averaged mass smeared along a Keplerian
/// ellipse of semi-major axis `a`, eccentricity `e`, apse `varpi`.
#[derive(Debug, Clone, Copy)]
pub struct Ring {
    pub a: f64,
    pub e: f64,
    pub varpi: f64,
    pub mass_solar: f64,
}

impl Ring {
    /// Ring semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }

    /// Ring apse as a typed [`Angle`].
    pub fn varpi_angle(&self) -> Angle {
        radians(self.varpi)
    }

    /// Ring mass as a typed [`Mass`].
    pub fn mass(&self) -> Mass {
        solar_masses(self.mass_solar)
    }
}

/// Disk surface density Sigma(a) (solar masses / AU^2), for documentation/plots.
pub fn surface_density(profile: &DiskProfile, a: f64) -> f64 {
    if a < profile.a_in || a > profile.a_out {
        return 0.0;
    }
    // Normalize integral 2 pi int Sigma a da = M_disk over [a_in, a_out].
    let p = profile.sigma_index;
    let norm = if (p - 2.0).abs() < 1e-9 {
        (profile.a_out / profile.a_in).ln()
    } else {
        (profile.a_out.powf(2.0 - p) - profile.a_in.powf(2.0 - p)) / (2.0 - p)
    };
    let sigma0 = profile.mass_solar / (2.0 * PI * norm);
    sigma0 * a.powf(-p)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_ring_masses_sum_to_total() {
        let d = DiskProfile::fiducial(10.0);
        let total: f64 = d.rings().iter().map(|r| r.mass_solar).sum();
        assert_relative_eq!(total, d.mass_solar, epsilon = 1e-12 * d.mass_solar);
    }

    #[test]
    fn test_typed_accessors_match_f64() {
        use uom::si::angle::radian;
        use uom::si::length::astronomical_unit;
        use uom::si::mass::kilogram;
        let d = DiskProfile::fiducial(10.0);
        assert_relative_eq!(
            d.inner_radius().get::<astronomical_unit>(),
            d.a_in,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            d.outer_radius().get::<astronomical_unit>(),
            d.a_out,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            d.mass().get::<kilogram>(),
            d.mass_solar * p9_core::units::SOLAR_MASS_KG,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            d.varpi_disk_angle().get::<radian>(),
            d.varpi_disk,
            epsilon = 1e-12
        );
        let r = d.rings()[0];
        assert_relative_eq!(
            r.semi_major_axis().get::<astronomical_unit>(),
            r.a,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            r.mass().get::<kilogram>(),
            r.mass_solar * p9_core::units::SOLAR_MASS_KG,
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_mass_earth_roundtrip() {
        let d = DiskProfile::fiducial(7.5);
        assert_relative_eq!(d.mass_earth(), 7.5, epsilon = 1e-12);
    }

    #[test]
    fn test_rings_span_disk_and_share_apse() {
        let d = DiskProfile::fiducial(10.0);
        let rings = d.rings();
        assert!(rings.first().unwrap().a >= d.a_in);
        assert!(rings.last().unwrap().a <= d.a_out);
        for r in &rings {
            assert_eq!(r.varpi, d.varpi_disk);
            assert_eq!(r.e, d.e_disk);
        }
    }

    #[test]
    fn test_surface_density_integrates_to_mass() {
        let d = DiskProfile::fiducial(10.0);
        // 2 pi int_{a_in}^{a_out} Sigma a da == M_disk (numerical check).
        let n = 20000;
        let mut s = 0.0;
        let da = (d.a_out - d.a_in) / n as f64;
        for k in 0..n {
            let a = d.a_in + (k as f64 + 0.5) * da;
            s += 2.0 * PI * surface_density(&d, a) * a * da;
        }
        assert_relative_eq!(s, d.mass_solar, epsilon = 1e-3 * d.mass_solar);
    }
}

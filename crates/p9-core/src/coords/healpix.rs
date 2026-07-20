//! Minimal HEALPix RING-scheme ang2pix (Górski et al. 2005), enough to index
//! the rubin-watch coverage maps without a healpix dependency.
//!
//! Standard equatorial/polar-cap formulas; validated against `healpy`
//! (canonical chealpix) on 20 fixture directions including poles, the
//! equator and the RA seam (tests below).

/// Number of pixels for an nside.
pub fn npix(nside: u32) -> usize {
    12 * (nside as usize) * (nside as usize)
}

/// RING-scheme pixel index for an equatorial direction (degrees).
pub fn ang2pix_ring(nside: u32, ra_deg: f64, dec_deg: f64) -> usize {
    let ns = nside as i64;
    // z via the colatitude, mirroring healpy's theta path (the coverage
    // builder uses healpy; astropy-healpix has a different boundary
    // convention at the exact equator and is deliberately NOT the
    // reference).
    let z = (std::f64::consts::FRAC_PI_2 - dec_deg.to_radians()).cos();
    let za = z.abs();
    let phi = ra_deg.rem_euclid(360.0).to_radians();
    let tt = phi / std::f64::consts::FRAC_PI_2; // in [0, 4)

    if za <= 2.0 / 3.0 {
        // Equatorial region.
        let temp1 = ns as f64 * (0.5 + tt);
        let temp2 = ns as f64 * z * 0.75;
        let jp = (temp1 - temp2) as i64; // ascending-edge line index
        let jm = (temp1 + temp2) as i64; // descending-edge line index
        let ir = ns + 1 + jp - jm; // ring counter, 1-based within region
        let kshift = 1 - (ir & 1);
        let t1 = jp + jm - ns + kshift + 1;
        let mut ip = t1 / 2;
        ip %= 4 * ns;
        let ncap = 2 * ns * (ns - 1);
        (ncap + (ir - 1) * 4 * ns + ip) as usize
    } else {
        // Polar caps.
        let tp = tt - tt.floor();
        let tmp = ns as f64 * (3.0 * (1.0 - za)).sqrt();
        let jp = (tp * tmp) as i64;
        let jm = ((1.0 - tp) * tmp) as i64;
        let ir = jp + jm + 1; // ring number from the pole
        let ip = ((tt * ir as f64) as i64) % (4 * ir);
        if z > 0.0 {
            (2 * ir * (ir - 1) + ip) as usize
        } else {
            (12 * ns * ns - 2 * ir * (ir + 1) + ip) as usize
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Ground truth from healpy (nside 64, ring), including poles, the
    /// equator, and the RA wrap seam.
    const FIXTURES: [(f64, f64, usize); 20] = [
        (313.289713, 75.150981, 829),
        (103.254195, 34.695831, 10697),
        (217.133334, -21.641199, 33818),
        (279.91227, 24.246767, 14663),
        (257.786867, -23.680609, 34359),
        (329.536843, 28.783217, 12906),
        (309.741714, -26.136142, 35420),
        (330.565546, 34.414919, 10603),
        (9.571584, 77.286182, 545),
        (157.40928, 76.526331, 643),
        (174.579958, 49.971546, 5828),
        (23.455511, 55.651411, 4336),
        (2.025378, 24.608384, 14209),
        (299.023717, 6.257313, 22100),
        (0.0, 89.99, 0),
        (359.999, -89.99, 49151),
        (90.0, 0.0, 24512),
        (180.0, 0.001, 24320),
        (45.0, 41.810315, 7843),
        (270.0, -41.810315, 41024),
    ];

    #[test]
    fn matches_astropy_healpix_fixtures() {
        for &(ra, dec, expect) in &FIXTURES {
            let got = ang2pix_ring(64, ra, dec);
            assert_eq!(got, expect, "ang2pix(64, {ra}, {dec})");
        }
    }

    #[test]
    fn every_pixel_index_in_range() {
        let n = npix(64);
        let mut hit = std::collections::HashSet::new();
        for i in 0..73 {
            for j in 0..37 {
                let ra = i as f64 * 4.93;
                let dec = -90.0 + j as f64 * 4.87;
                let p = ang2pix_ring(64, ra, dec.clamp(-89.999, 89.999));
                assert!(p < n);
                hit.insert(p);
            }
        }
        assert!(
            hit.len() > 1000,
            "grid should scatter widely: {}",
            hit.len()
        );
    }
}

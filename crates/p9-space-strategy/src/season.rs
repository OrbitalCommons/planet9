//! When a field can be observed: a distant body is best searched near
//! opposition, where it is brightest, the sky behind it is darkest of
//! sunlight, and the reflex of Earth's motion moves it fastest.

/// Julian date of J2000.0.
const J2000_JD: f64 = 2_451_545.0;
/// Days in a tropical year.
const YEAR: f64 = 365.2422;

/// Mean ecliptic longitude of the Sun (deg) at Julian date `jd`.
pub fn sun_longitude_deg(jd: f64) -> f64 {
    (280.460 + 0.985_647_4 * (jd - J2000_JD)).rem_euclid(360.0)
}

/// Day of year (0–365) on which a field at ecliptic longitude `lon_deg` is at
/// opposition.
pub fn opposition_day_of_year(lon_deg: f64) -> f64 {
    // Sun at lon + 180°. Its longitude is 280.46° on 1 January.
    ((lon_deg + 180.0 - 280.46).rem_euclid(360.0)) / 360.0 * YEAR
}

/// Calendar month (1–12) of opposition.
pub fn opposition_month(lon_deg: f64) -> u32 {
    const CUMULATIVE: [f64; 12] = [
        31.0, 59.25, 90.25, 120.25, 151.25, 181.25, 212.25, 243.25, 273.25, 304.25, 334.25, 366.0,
    ];
    let d = opposition_day_of_year(lon_deg);
    CUMULATIVE.iter().position(|&c| d < c).unwrap_or(11) as u32 + 1
}

pub const MONTHS: [&str; 12] = [
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];

/// Days either side of opposition during which a field's reflex motion stays
/// above `fraction` of its opposition value.
pub fn window_half_width_days(fraction: f64) -> f64 {
    fraction.clamp(0.0, 1.0).acos().to_degrees() / 360.0 * YEAR
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sun_is_at_the_equinox_in_march() {
        // 20 March 2026 ≈ JD 2461119.5. The mean longitude trails the true
        // Sun by the equation of centre (~2° in March).
        let l = sun_longitude_deg(2_461_119.5);
        assert!(!(3.0..=357.0).contains(&l), "sun longitude {l}");
    }

    #[test]
    fn opposition_seasons() {
        // Ecliptic longitude 90° (RA 6h, Gemini/Taurus) is opposite the Sun
        // at the December solstice; longitude 270° in June.
        assert_eq!(opposition_month(90.0), 12);
        assert_eq!(opposition_month(270.0), 6);
        assert_eq!(opposition_month(70.0), 12);
        assert_eq!(opposition_month(45.0), 11);
    }

    #[test]
    fn the_window_is_about_three_months_wide() {
        // Motion stays above 70% of its peak for ±46 days.
        assert!((window_half_width_days(0.7) - 46.2).abs() < 0.5);
    }
}

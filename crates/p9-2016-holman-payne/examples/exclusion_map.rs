//! Dump the Holman & Payne (2016) mass–distance exclusion map and favored sky
//! position.
//!
//! Run with: `cargo run -p p9-2016-holman-payne --example exclusion_map`

use p9_2016_holman_payne::{
    exclusion::ExclusionMap, favored_sky_position, published, range_residual_m,
};
use p9_core::data::ephemeris_constraint::brown_batygin_orbit;

fn main() {
    let p = brown_batygin_orbit();
    println!(
        "P9 (Brown & Batygin): a={} AU, e={}, m={} M⊕; Cassini precision ~{} m",
        p.a,
        p.e,
        p.mass_earth,
        published::CASSINI_RANGE_PRECISION_M
    );

    println!("\nMass–distance exclusion boundary (max allowed mass at each a):");
    println!("{:>8}  {:>16}", "a[AU]", "max allowed[M⊕]");
    let map = ExclusionMap::build(&p, 300.0, 1500.0, 13);
    for (a, m) in map.a_au.iter().zip(map.max_allowed_mass_earth.iter()) {
        println!("{a:>8.0}  {m:>16.3}");
    }

    println!("\nRange residual vs true anomaly (nominal 10 M⊕, a = 700 AU):");
    println!("{:>5}  {:>14}", "ν°", "residual[m]");
    for k in 0..18 {
        let nu = ((k * 20) as f64).to_radians();
        println!(
            "{:>5.0}  {:>14.3e}",
            k as f64 * 20.0,
            range_residual_m(&p, nu)
        );
    }

    let sky = favored_sky_position(&p, 1440);
    println!(
        "\nFavored sky position: RA = {:.1}°, Dec = {:.1}° (ν = {:.1}°, r = {:.0} AU)",
        sky.ra_deg,
        sky.dec_deg,
        sky.true_anomaly.to_degrees(),
        sky.distance_au
    );
    println!(
        "Published preferred region: RA ≈ {:.0}°, Dec ≈ {:.0}° (±{:.0}°); separation {:.0}°",
        published::PREFERRED_RA_DEG,
        published::PREFERRED_DEC_DEG,
        published::PREFERRED_SKY_HALF_EXTENT_DEG,
        sky.separation_deg(published::PREFERRED_RA_DEG, published::PREFERRED_DEC_DEG)
    );
}

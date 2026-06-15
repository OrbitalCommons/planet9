//! Dump the Fienga et al. (2016) range-perturbation curve over P9 true anomaly.
//!
//! Run with: `cargo run -p p9-2016-cassini-ranging --example range_curve`

use p9_2016_cassini_ranging::{
    favored_true_anomaly, perturbation::prefit_amplitude, range_perturbation_amplitude,
};
use p9_core::data::ephemeris_constraint::brown_batygin_orbit;
use p9_core::types::position_at_true_anomaly as p9_position_at_true_anomaly;

fn main() {
    let p = brown_batygin_orbit();
    println!(
        "P9 (B&B / Fienga Table 1): a={} AU, e={}, i={:.0}°, ω={:.0}°, Ω={:.0}°, m={} M⊕",
        p.a,
        p.e,
        p.i.to_degrees(),
        p.omega.to_degrees(),
        p.omega_big.to_degrees(),
        p.mass_earth
    );
    println!(
        "{:>5}  {:>8}  {:>12}  {:>12}",
        "ν°", "r[AU]", "post-fit[km]", "pre-fit[km]"
    );
    for k in 0..36 {
        let nu = ((k * 10) as f64).to_radians();
        let post = range_perturbation_amplitude(&p, nu);
        let pre = prefit_amplitude(&p, nu);
        let r = p9_position_at_true_anomaly(&p, nu).distance;
        println!(
            "{:>5.0}  {:>8.1}  {:>12.4e}  {:>12.4e}",
            k as f64 * 10.0,
            r,
            post,
            pre
        );
    }
    let fav = favored_true_anomaly(&p, 1440).to_degrees();
    let r_fav = p9_position_at_true_anomaly(&p, fav.to_radians()).distance;
    println!("\nfavored ν = {fav:.1}° (r = {r_fav:.0} AU); published 117.8°, ~600 AU");
}

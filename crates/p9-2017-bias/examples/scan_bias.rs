//! Parameter scan for the simplified bias model: reports the Monte Carlo
//! clustering p-values for a grid of bias-model assumptions, used to check
//! how far the defaults sit from the paper's p_combined = 2.5e-4.

use p9_2017_bias::bias_function::BiasParams;
use p9_2017_bias::clustering_test::monte_carlo_clustering_test_with_params;
use p9_2017_bias::kbo_sample::paper_sample_a230;
use p9_core::constants::DEG2RAD;

fn main() {
    let kbos = paper_sample_a230();
    // (sigma_anti, sigma_center, floor_anti, floor_center, south_penalty)
    for (s_a, s_c, f_a, f_c, south) in [
        (15.0, 30.0, 0.35, 0.02, 0.5),
        (15.0, 30.0, 0.35, 0.02, 1.0),
        (15.0, 40.0, 0.5, 0.02, 0.5),
        (15.0, 40.0, 0.5, 0.02, 0.3),
        (20.0, 45.0, 0.5, 0.0, 0.4),
        (15.0, 35.0, 0.4, 0.0, 0.5),
    ] {
        let params = BiasParams {
            galactic_lon_sigmas: [s_a * DEG2RAD, s_c * DEG2RAD],
            galactic_lon_floors: [f_a, f_c],
            south_penalty: south,
            ..BiasParams::default()
        };
        let r = monte_carlo_clustering_test_with_params(&kbos, 200_000, 42, &params);
        println!(
            "s=({s_a},{s_c}) f=({f_a},{f_c}) south={south} -> p_varpi={:.2e} p_omega={:.2e} p_comb={:.2e}",
            r.p_varpi, r.p_omega, r.p_combined
        );
    }
}

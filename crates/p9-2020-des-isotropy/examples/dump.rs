use p9_2020_des_isotropy::analysis::run_battery;
use p9_2020_des_isotropy::selection::NullModel;

fn main() {
    for null in [NullModel::Uniform, NullModel::DesSelection] {
        println!("=== null = {null:?} ===");
        let cells = run_battery(null, 2020, 8000);
        for c in &cells {
            let p = &c.p;
            println!(
                "{:>6} {:>6} n={} R={:.3} rayl={:.3} kuip={:.3} wat={:.3} ber={:.3}",
                c.case, c.angle, p.n, p.r_bar, p.rayleigh_p, p.kuiper_p, p.watson_p, p.beran_p
            );
        }
    }
}

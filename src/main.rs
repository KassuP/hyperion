use hyperion::euler::euler;
use hyperion::calculator::*;
use hyperion::plot::*;
use std::f64::consts::PI;

fn main() {
    let theta0 :f64= 0.0;
    let dttheta = 1.0e-5;
    let dtheta = dttheta + theta0;
    let e :f64= 0.3;
    let dt :f64= 200.0;
    let t_end :f64= 1.0e6;
    let de :f64= 0.002;

    // euler( theta0, e, dt, t_end,"src/data/hyperion.csv");
    // // euler( dtheta, e, dt, t_end,"src/data/dhyperion.csv");
    // plot_hyperion_theta_cir();
    // plot_hyperion_omega();
    // plot_hyperion_delta_theta("src/data/hyperion.csv","src/data/dhyperion.csv");
    // plot_hyperion_double_theta("src/data/hyperion.csv","src/data/dhyperion.csv");
    // plot_hyperion_double_omega("src/data/hyperion.csv","src/data/dhyperion.csv");
    // let lyapunov = plot_hyperion_improved_Lyapunov("src/data/hyperion.csv","src/data/dhyperion.csv").unwrap();
    // println!("最终Lyapunov指数: {:.6e}", lyapunov);
    plot_hyperion_lyapunov_e(&euler_scan(theta0,dtheta,dt,t_end,de));
}

use std::fs::File;
use std::io::Write;

pub fn euler(
    mut theta: f64,  // 初始角位置 θ0
    e: f64,          // 偏心率（不需要mut，因为不会改变）
    dt: f64,         // 时间步长 Δt
    t_end: f64,      // 总时间
    file: &str
) {
    let g: f64 = 6.67430e-11;     // 万有引力常数
    let m_sat: f64 = 5.683e26;   // 例如土星质量
    let a: f64 = 1.5e8;
    
    // 声明为可变变量
    let mut w: f64 = 1e-4;      // 初始角速度（与Python代码一致）
    let mut vx: f64 = 0.0;
    let mut vy: f64 = ((g * m_sat) * (1.0 + e) / (a * (1.0 - e))).sqrt(); // 近地点速度
    let mut y: f64 = 0.0; 
    let mut x: f64 = a * (1.0 - e); // 近地点位置

    // 创建输出文件
    std::fs::create_dir_all("data").unwrap();
    let mut output_file = File::create(file).expect("无法创建输出文件");

    // 写入表头
    writeln!(output_file, "t,x,y,vx,vy,r,theta,omega").unwrap();

    let mut t = 0.0;

    while t <= t_end {
        let r = (x * x + y * y).sqrt();

        // ---- Euler-Cromer 更新 ----
        // 更新速度
        vx -= g * m_sat * x / (r.powi(3)) * dt;
        vy -= g * m_sat * y / (r.powi(3)) * dt;

        // 更新位置
        x += vx * dt;
        y += vy * dt;

        // 更新角速度
        w -= (3.0 * g * m_sat / r.powi(5))
            * (x * theta.sin() - y * theta.cos())
            * (x * theta.cos() + y * theta.sin())
            * dt;

        // 更新角位置
        theta += w * dt;

        // 写入数据
        writeln!(output_file, "{},{},{},{},{},{},{},{}",
            t, x, y, vx, vy, r, theta, w
        ).unwrap();

        t += dt;
    }

    println!("✅ 计算完成，数据已保存到 {}", file);
}
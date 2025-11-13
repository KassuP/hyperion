use std::fs::File;
use std::io::{BufReader, BufRead};
use std::error::Error;
use std::f64::consts::PI;

pub fn plot_hyperion_improved_Lyapunov(file1: &str, file2: &str) -> Result<f64, Box<dyn std::error::Error>> {
    // 读取CSV文件
    let reader1 = BufReader::new(File::open(file1)?);
    let reader2 = BufReader::new(File::open(file2)?);

    // -------------------- 工具函数 --------------------
    /// 将角度限制到 [0, 2π)
    fn wrap_to_2pi(angle: f64) -> f64 {
        let mut wrapped = angle % (2.0 * PI);
        if wrapped < 0.0 {
            wrapped += 2.0 * PI;
        }
        wrapped
    }

    /// 实现类似 np.unwrap 的逻辑，使角度连续
    fn unwrap_angles(angles: &[f64]) -> Vec<f64> {
        if angles.is_empty() {
            return vec![];
        }
        let mut unwrapped = Vec::with_capacity(angles.len());
        unwrapped.push(angles[0]);
        for i in 1..angles.len() {
            let mut delta = angles[i] - angles[i - 1];
            if delta > PI {
                delta -= 2.0 * PI;
            } else if delta < -PI {
                delta += 2.0 * PI;
            }
            unwrapped.push(unwrapped[i - 1] + delta);
        }
        unwrapped
    }

    /// 线性回归计算斜率
    fn linear_regression(x: &[f64], y: &[f64]) -> Result<(f64, f64), Box<dyn Error>> {
        if x.len() != y.len() || x.len() < 2 {
            return Err("Insufficient data for linear regression".into());
        }
        
        let n = x.len() as f64;
        let sum_x: f64 = x.iter().sum();
        let sum_y: f64 = y.iter().sum();
        let sum_xy: f64 = x.iter().zip(y.iter()).map(|(&a, &b)| a * b).sum();
        let sum_x2: f64 = x.iter().map(|&a| a * a).sum();
        
        let denominator = n * sum_x2 - sum_x * sum_x;
        if denominator.abs() < 1e-15 {
            return Err("Cannot compute linear regression: denominator too small".into());
        }
        
        let slope = (n * sum_xy - sum_x * sum_y) / denominator;
        let intercept = (sum_y - slope * sum_x) / n;
        
        Ok((slope, intercept))
    }

    // -------------------- 数据读取 --------------------
    let mut theta1_positions = Vec::new();
    let mut theta2_positions = Vec::new();
    let mut t_positions = Vec::new();

    // 读取文件1
    let mut lines = reader1.lines();
    lines.next(); // 跳过表头
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 7 {
            let theta = parts[6].parse::<f64>()?;
            let t = parts[0].parse::<f64>()?;
            theta1_positions.push(wrap_to_2pi(theta));
            t_positions.push(t);
        }
    }

    // 读取文件2
    let mut lines = reader2.lines();
    lines.next(); // 跳过表头
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 7 {
            let theta = parts[6].parse::<f64>()?;
            theta2_positions.push(wrap_to_2pi(theta));
        }
    }

    // -------------------- unwrap + 差值 --------------------
    let theta1_u = unwrap_angles(&theta1_positions);
    let theta2_u = unwrap_angles(&theta2_positions);
    let delta_theta_positions: Vec<f64> = theta1_u
        .iter()
        .zip(theta2_u.iter())
        .map(|(a, b)| (a - b).abs())
        .collect();

    // -------------------- 线性回归方法计算Lyapunov指数 --------------------
    
    // 检查数据长度
    if t_positions.len() != delta_theta_positions.len() || t_positions.len() < 2 {
        return Err("Insufficient data for Lyapunov calculation".into());
    }

    let n = t_positions.len();
    
    // 第一步：设置初始阈值并创建mask
    let mut low_th = 1e-10;
    let mut high_th = 20.0;
    let mut mask: Vec<bool> = delta_theta_positions.iter()
        .map(|&x| x > low_th && x < high_th)
        .collect();
    
    // 检查有效点数
    let mut valid_count = mask.iter().filter(|&&x| x).count();
    // println!("初始阈值筛选有效点数: {}", valid_count);
    
    // 第二步：如果点数太少，放宽阈值
    if valid_count < 10 {
        println!("放宽阈值...");
        low_th = 1e-12;
        high_th = 1e-1;
        mask = delta_theta_positions.iter()
            .map(|&x| x > low_th && x < high_th)
            .collect();
        valid_count = mask.iter().filter(|&&x| x).count();
        // println!("放宽后阈值筛选有效点数: {}", valid_count);
    }
    
    // 第三步：如果仍然太少，使用中间20%
    let final_mask = if valid_count < 10 {
        println!("使用中间20%数据...");
        let start = (0.4 * n as f64) as usize;
        let end = (0.6 * n as f64) as usize;
        (0..n).map(|i| i >= start && i < end).collect()
    } else {
        mask
    };
    
    // 提取筛选后的数据
    let t_filtered: Vec<f64> = t_positions.iter().zip(&final_mask)
        .filter_map(|(&t_val, &mask_val)| if mask_val { Some(t_val) } else { None })
        .collect();
    
    let ln_delta_filtered: Vec<f64> = delta_theta_positions.iter().zip(&final_mask)
        .filter_map(|(&dt_val, &mask_val)| {
            if mask_val && dt_val > 0.0 { 
                Some(dt_val.ln()) 
            } else { 
                None 
            }
        })
        .collect();
    
    // println!("最终用于线性回归的数据点数: {}", t_filtered.len());
    
    // 线性回归计算Lyapunov指数
    let lyapunov_exponent = if t_filtered.len() >= 2 {
        match linear_regression(&t_filtered, &ln_delta_filtered) {
            Ok((slope, _intercept)) => {
                // println!("线性回归斜率 (λ): {:.6e}", slope);
                slope
            }
            Err(e) => {
                return Err(format!("线性回归失败: {}", e).into());
            }
        }
    } else {
        return Err("线性回归数据点不足".into());
    };

    // -------------------- 输出结果 --------------------
    // println!("\n=== Lyapunov指数计算结果 (线性回归方法) ===");
    // println!("使用的数据筛选方法: {}", 
    //     if valid_count >= 10 { "阈值筛选" } else { "中间20%" });
    // println!("有效数据点数: {}", t_filtered.len());
    // println!("Lyapunov指数 λ: {:.6e}", lyapunov_exponent);
    
    // 判断系统的性质
    // if lyapunov_exponent > 1e-6 {
    //     println!("系统性质: 混沌 (λ > 0)");
    // } else if lyapunov_exponent < -1e-6 {
    //     println!("系统性质: 稳定 (λ < 0)");
    // } else {
    //     println!("系统性质: 中性 (λ ≈ 0)");
    // }
    
    // 单位转换：如果需要每天的Lyapunov指数
    let lyapunov_per_day = lyapunov_exponent * 86400.0; // 秒/天 = 86400
    // println!("Lyapunov指数 (每天): {:.6e}", lyapunov_per_day);
    // println!("============================");

    Ok(lyapunov_exponent)
}


pub fn lyapunov_from_vectors(
    data1: &Vec<(f64, f64)>,  // 系统1: (t, θ)
    data2: &Vec<(f64, f64)>,  // 系统2: (t, θ)
) -> Result<f64, Box<dyn std::error::Error>> {
    use std::f64::consts::PI;

    // -------------------- 工具函数 --------------------
    fn wrap_to_2pi(angle: f64) -> f64 {
        let mut wrapped = angle % (2.0 * PI);
        if wrapped < 0.0 {
            wrapped += 2.0 * PI;
        }
        wrapped
    }

    fn unwrap_angles(angles: &[f64]) -> Vec<f64> {
        if angles.is_empty() {
            return vec![];
        }
        let mut unwrapped = Vec::with_capacity(angles.len());
        unwrapped.push(angles[0]);
        for i in 1..angles.len() {
            let mut delta = angles[i] - angles[i - 1];
            if delta > PI {
                delta -= 2.0 * PI;
            } else if delta < -PI {
                delta += 2.0 * PI;
            }
            unwrapped.push(unwrapped[i - 1] + delta);
        }
        unwrapped
    }

    fn linear_regression(x: &[f64], y: &[f64]) -> Result<(f64, f64), Box<dyn std::error::Error>> {
        if x.len() != y.len() || x.len() < 2 {
            return Err("Insufficient data for linear regression".into());
        }
        let n = x.len() as f64;
        let sum_x: f64 = x.iter().sum();
        let sum_y: f64 = y.iter().sum();
        let sum_xy: f64 = x.iter().zip(y.iter()).map(|(&a, &b)| a * b).sum();
        let sum_x2: f64 = x.iter().map(|&a| a * a).sum();

        let denominator = n * sum_x2 - sum_x * sum_x;
        if denominator.abs() < 1e-15 {
            return Err("Cannot compute linear regression: denominator too small".into());
        }
        let slope = (n * sum_xy - sum_x * sum_y) / denominator;
        let intercept = (sum_y - slope * sum_x) / n;
        Ok((slope, intercept))
    }

    // -------------------- 数据准备 --------------------
    if data1.len() != data2.len() || data1.len() < 2 {
        return Err("Data length mismatch or insufficient points".into());
    }

    let t_positions: Vec<f64> = data1.iter().map(|(t, _)| *t).collect();
    let theta1_positions: Vec<f64> = data1.iter().map(|(_, theta)| wrap_to_2pi(*theta)).collect();
    let theta2_positions: Vec<f64> = data2.iter().map(|(_, theta)| wrap_to_2pi(*theta)).collect();

    // unwrap 角度
    let theta1_u = unwrap_angles(&theta1_positions);
    let theta2_u = unwrap_angles(&theta2_positions);

    // 计算 Δθ
    let delta_theta_positions: Vec<f64> = theta1_u
        .iter()
        .zip(theta2_u.iter())
        .map(|(a, b)| (a - b).abs())
        .collect();

    // -------------------- 数据筛选与拟合 --------------------
    let n = t_positions.len();
    let low_th = 1e-10;
    let high_th = 20.0;
    let mask: Vec<bool> = delta_theta_positions
        .iter()
        .map(|&x| x > low_th && x < high_th)
        .collect();

    let valid_count = mask.iter().filter(|&&x| x).count();
    let final_mask = if valid_count < 10 {
        let start = (0.4 * n as f64) as usize;
        let end = (0.6 * n as f64) as usize;
        (0..n).map(|i| i >= start && i < end).collect::<Vec<_>>()
    } else {
        mask
    };

    let t_filtered: Vec<f64> = t_positions
        .iter()
        .zip(&final_mask)
        .filter_map(|(&t_val, &mask_val)| if mask_val { Some(t_val) } else { None })
        .collect();

    let ln_delta_filtered: Vec<f64> = delta_theta_positions
        .iter()
        .zip(&final_mask)
        .filter_map(|(&dt_val, &mask_val)| {
            if mask_val && dt_val > 0.0 {
                Some(dt_val.ln())
            } else {
                None
            }
        })
        .collect();

    if t_filtered.len() < 2 {
        return Err("Insufficient data after filtering for regression".into());
    }

    let (slope, _intercept) = linear_regression(&t_filtered, &ln_delta_filtered)?;
    Ok(slope)
}

pub fn euler_base(
    mut theta: f64,  // 初始角位置 θ0
    dt: f64,         // 时间步长 Δt
    e : f64,
    t_end: f64,      // 总时间
)-> Vec<(f64, f64)> {
    let g: f64 = 6.67430e-11;     // 万有引力常数
    let m_sat: f64 = 5.683e26;   // 例如土星质量
    let a: f64 = 1.5e8;
    
    // 声明为可变变量
    let mut w: f64 = 1e-4;      // 初始角速度（与Python代码一致）
    let mut vx: f64 = 0.0;
    let mut vy: f64 = ((g * m_sat) * (1.0 + e) / (a * (1.0 - e))).sqrt(); // 近地点速度
    let mut y: f64 = 0.0; 
    let mut x: f64 = a * (1.0 - e); // 近地点位置


    let mut t = 0.0;
    let mut data = Vec::new();

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

        data.push((t, theta));
        // 更新角位置
        theta += w * dt;
        t += dt;
    }
    data
}

pub fn euler_improved(
    mut theta: f64,  // 初始角位置 θ0
    mut dtheta: f64, // 初始微小角度差 Δθ0
    dt: f64,         // 时间步长 Δt
    t_end: f64,      // 总时间
    e: f64,
) -> Result<(f64, f64), Box<dyn std::error::Error>> {
    let data1 = euler_base(theta, dt, e, t_end);
    let data2 = euler_base(theta + dtheta, dt, e, t_end);
    let lambda = lyapunov_from_vectors(&data1, &data2)?;
    Ok((e, lambda))
}

pub fn euler_scan(    
    mut theta: f64,  // 初始角位置 θ0
    mut dtheta: f64, // 初始微小角度差 Δθ0
    dt: f64,         // 时间步长 Δt
    t_end: f64,      // 总时间
    de :f64,
)-> Vec<(f64, f64)>{
    let mut results = Vec::new();  // 存储 (e, λ) 的结果
    let mut e = 0.0;
    while e < 0.9 {
        let (eccentricity, lyapunov) = euler_improved(theta, dtheta, dt, t_end, e).unwrap();
        println!("偏心率: {:.4}, Lyapunov指数: {:.6e}", eccentricity, lyapunov);
        results.push((eccentricity, lyapunov));
        e += de;
    }
    results
}
use plotly::common::{Line,Mode};
use plotly::{Plot, Scatter};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::f64::consts::PI;
use crate::calculator::*;

pub fn plot_hyperion_position(l: f64) -> Result<(), Box<dyn std::error::Error>> {
    // 读取CSV文件
    let path = "src/data/hyperion.csv";
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    // 存储坐标数据
    let mut x_positions = Vec::new();
    let mut y_positions = Vec::new();
    
    // 跳过表头
    let mut lines = reader.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let x: f64 = parts[1].parse()?;
            let y: f64 = parts[2].parse()?;
            let theta: f64 = parts[6].parse()?;
            
            // 计算新位置: (x + l*cos(theta), y + l*sin(theta))
            let x_new = x + l * theta.cos();
            let y_new = y + l * theta.sin();
            // let x_new = x ;
            // let y_new = y ;
            
            x_positions.push(x_new);
            y_positions.push(y_new);
        }
    }
    
    // 创建图表
    let trace = Scatter::new(x_positions, y_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0)) 
        .name("Hyperion Position");
    
    let mut plot = Plot::new();
    plot.add_trace(trace);
    
    // 设置图表布局
    plot.set_layout(
        plotly::Layout::new()
            .title("Hyperion Position".into())
            .x_axis(plotly::layout::Axis::new().title("X Position".into()))
            .y_axis(plotly::layout::Axis::new().title("Y Position".into()))
    );
    
    // 保存为HTML文件
    let output_path = "src/plot/hyperion_position_plot.html";
    plot.write_html(output_path);
    
    println!("✅ 图表已保存到 {}", output_path);
    Ok(())
}


pub fn plot_hyperion_theta_cir() -> Result<(), Box<dyn std::error::Error>> {
    // 读取CSV文件
    let path = "src/data/hyperion.csv";
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    // 存储坐标数据
    let mut theta_plot = Vec::new();
    let mut t_plot = Vec::new();
    
    // 跳过表头
    let mut lines = reader.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let theta_rad: f64 = parts[6].parse()?;
            let t: f64 = parts[0].parse()?;
            
            // 弧度转角度
            let mut theta_deg = theta_rad * 180.0 / PI;
            
            // 循环边界 [-360, 360]
            while theta_deg > 360.0 {
                theta_deg -= 720.0;
            }
            while theta_deg < -360.0 {
                theta_deg += 720.0;
            }
            
            theta_plot.push(theta_deg);
            t_plot.push(t);
        }
    }
    
    // 创建图表
    let trace = Scatter::new(t_plot, theta_plot)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0))
        .name("Hyperion Theta");
    
    let mut plot = Plot::new();
    plot.add_trace(trace);
    
    // 设置图表布局
    plot.set_layout(
        plotly::Layout::new()
            .title("Hyperion Theta".into())
            .x_axis(plotly::layout::Axis::new().title("t".into()))
            .y_axis(plotly::layout::Axis::new().title("Theta (deg)".into()))
    );
    
    // 保存为HTML文件
    let output_path = "src/plot/hyperion_theta_plot.html";
    plot.write_html(output_path);
    
    println!("✅ 图表已保存到 {}", output_path);
    Ok(())
}

pub fn plot_hyperion_omega() -> Result<(), Box<dyn std::error::Error>> {
    // 读取CSV文件
    let path = "src/data/hyperion.csv";
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    
    // 存储坐标数据
    let mut omega_positions = Vec::new();
    let mut t_positions = Vec::new();
    
    // 跳过表头
    let mut lines = reader.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let x: f64 = parts[7].parse()?;
            let y: f64 = parts[0].parse()?;
            
            omega_positions.push(x);
            t_positions.push(y);
        }
    }
    
    // 创建图表
    let trace = Scatter::new(t_positions, omega_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0)) 
        .name("Hyperion Omega");
    
    let mut plot = Plot::new();
    plot.add_trace(trace);
    
    // 设置图表布局
    plot.set_layout(
        plotly::Layout::new()
            .title("Hyperion Position".into())
            .x_axis(plotly::layout::Axis::new().title("t".into()))
            .y_axis(plotly::layout::Axis::new().title("Omaga".into()))
    );
    
    // 保存为HTML文件
    let output_path = "src/plot/hyperion_omega_plot.html";
    plot.write_html(output_path);
    
    println!("✅ 图表已保存到 {}", output_path);
    Ok(())
}

pub fn plot_hyperion_double_theta(file1:&str, file2:&str) -> Result<(), Box<dyn std::error::Error>> {
    // 读取CSV文件
    let file1 = File::open(file1)?;
    let reader1 = BufReader::new(file1);
    
    let file2 = File::open(file2)?;
    let reader2 = BufReader::new(file2);

    // 存储坐标数据
    let mut theta1_positions = Vec::new();
    let mut theta2_positions = Vec::new();
    let mut t_positions = Vec::new();
    
    // 跳过表头------------------------------
    let mut lines = reader1.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let x: f64 = parts[6].parse()?;
            let y: f64 = parts[0].parse()?;
            
            theta1_positions.push(x);
            t_positions.push(y);
        }
    }

    // -----------------------------------------
    let mut lines = reader2.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let x: f64 = parts[6].parse()?;
            let y: f64 = parts[0].parse()?;
            
            theta2_positions.push(x);
        }
    }
    
    // 创建图表
    let trace1 = Scatter::new(t_positions.clone(), theta1_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0).color("#20a094ff")) 
        .name("Hyperion theta");

    let trace2 = Scatter::new(t_positions, theta2_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0).color("#ff9705ff").dash(plotly::common::DashType::Dash)) 
        .name("Hyperion dtheta");

    let mut plot = Plot::new();
    plot.add_trace(trace1);
    plot.add_trace(trace2);
    
    // 设置图表布局
    plot.set_layout(
        plotly::Layout::new()
            .title("Hyperion Position".into())
            .x_axis(plotly::layout::Axis::new().title("t".into()))
            .y_axis(plotly::layout::Axis::new().title("theta".into()))
    );
    
    // 保存为HTML文件
    let output_path = "src/plot/hyperion_double_theta_plot.html";
    plot.write_html(output_path);
    
    println!("✅ 图表已保存到 {}", output_path);
    Ok(())
}

// pub fn plot_hyperion_delta_theta(file1:&str, file2:&str) -> Result<(), Box<dyn std::error::Error>> {
//     // 读取CSV文件
//     let file1 = File::open(file1)?;
//     let reader1 = BufReader::new(file1);
    
//     let file2 = File::open(file2)?;
//     let reader2 = BufReader::new(file2);

//     // 存储坐标数据
//     let mut theta1_positions = Vec::new();
//     let mut theta2_positions = Vec::new();
//     let mut t_positions = Vec::new();
    
//     // 辅助函数：将角度归一化到 [0, 2π]
//     fn wrap_to_2pi(angle: f64) -> f64 {
//         let mut wrapped = angle % (2.0 * PI);
//         if wrapped < 0.0 {
//             wrapped += 2.0 * PI;
//         }
//         wrapped
//     }
//     // 跳过表头------------------------------
//     let mut lines = reader1.lines();
//     lines.next(); // 跳过第一行表头
    
//     // 解析每一行数据
//     for line in lines {
//         let line = line?;
//         let parts: Vec<&str> = line.split(',').collect();
        
//         if parts.len() >= 7 {
//             let x: f64 = parts[6].parse()?;
//             let y: f64 = parts[0].parse()?;
            
//             theta1_positions.push(wrap_to_2pi(x));
//             t_positions.push(y);
//         }
//     }

//     // -----------------------------------------
//     let mut lines = reader2.lines();
//     lines.next(); // 跳过第一行表头
    
//     // 解析每一行数据
//     for line in lines {
//         let line = line?;
//         let parts: Vec<&str> = line.split(',').collect();
        
//         if parts.len() >= 7 {
//             let x: f64 = parts[6].parse()?;
//             let y: f64 = parts[0].parse()?;
            
//             theta2_positions.push(wrap_to_2pi(x));
//         }
//     }

//     let mut delta_theta_positions = Vec::new();
//     for (a, b) in theta1_positions.iter().zip(theta2_positions.iter()) {
//         delta_theta_positions.push((a - b).abs());
//     }

//     // 创建图表
//     let trace = Scatter::new(t_positions, delta_theta_positions)
//         .mode(Mode::Lines)
//         .line(Line::new().width(1.0).color("#35bdafff")) 
//         .name("Hyperion theta");


//     let mut plot = Plot::new();
//     plot.add_trace(trace);
    
//     // 设置图表布局
//     plot.set_layout(
//         plotly::Layout::new()
//             .title("Hyperion Delta theta".into())
//             .x_axis(plotly::layout::Axis::new().title("t".into()))
//             .y_axis(plotly::layout::Axis::new().title("Delta theta".into()))
//     );
    
//     // 保存为HTML文件
//     let output_path = "src/plot/hyperion_delta_theta_plot.html";
//     plot.write_html(output_path);
    
//     println!("✅ 图表已保存到 {}", output_path);
//     Ok(())  
// }
pub fn plot_hyperion_delta_theta(file1: &str, file2: &str) -> Result<(), Box<dyn std::error::Error>> {
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
    let delta_theta_0 = delta_theta_positions[0];
    // -------------------- 绘图 --------------------
    let trace = Scatter::new(t_positions.clone(), delta_theta_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.5).color("#1f77b4"))
        .name("Δθ");

    let lyapunov = plot_hyperion_improved_Lyapunov(file1,file2).unwrap();
    let theoretical_delta_theta: Vec<f64> = t_positions
        .iter()
        .map(|&t| delta_theta_0 * (lyapunov * t).exp())
        .collect();
    let trace_theoretical = Scatter::new(t_positions.clone(), theoretical_delta_theta)
        .mode(Mode::Lines)
        .line(Line::new().width(2.0).color("#ff7f0e").dash(plotly::common::DashType::Dash))
        .name(format!("Lyapunov exponent (λ = {:.4e})", lyapunov));

    let mut plot = Plot::new();
    plot.add_trace(trace);
    plot.add_trace(trace_theoretical);

    // 设置半对数坐标轴（y轴为对数坐标）
    plot.set_layout(
        plotly::Layout::new()
            .title("Hyperion Δθ over time (Semi-log plot)".into())
            .x_axis(plotly::layout::Axis::new().title("t".into()))
            .y_axis(
                plotly::layout::Axis::new()
                    .title("Δθ (rad)".into())
                    .type_(plotly::layout::AxisType::Log)  // 设置y轴为对数坐标
            ),
    );

    let output_path = "src/plot/hyperion_delta_theta_semilog_plot.html";
    plot.write_html(output_path);
    println!("✅ 半对数坐标图表已保存到 {}", output_path);

    Ok(())
}


pub fn plot_hyperion_double_omega(file1:&str, file2:&str) -> Result<(), Box<dyn std::error::Error>> {
    // 读取CSV文件
    let file1 = File::open(file1)?;
    let reader1 = BufReader::new(file1);
    
    let file2 = File::open(file2)?;
    let reader2 = BufReader::new(file2);

    // 存储坐标数据
    let mut theta1_positions = Vec::new();
    let mut theta2_positions = Vec::new();
    let mut t_positions = Vec::new();
    
    // 跳过表头------------------------------
    let mut lines = reader1.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let x: f64 = parts[7].parse()?;
            let y: f64 = parts[0].parse()?;
            
            theta1_positions.push(x);
            t_positions.push(y);
        }
    }

    // -----------------------------------------
    let mut lines = reader2.lines();
    lines.next(); // 跳过第一行表头
    
    // 解析每一行数据
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        
        if parts.len() >= 7 {
            let x: f64 = parts[7].parse()?;
            let y: f64 = parts[0].parse()?;
            
            theta2_positions.push(x);
        }
    }
    
    // 创建图表
    let trace1 = Scatter::new(t_positions.clone(), theta1_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0).color("#10dfcbff")) 
        .name("Hyperion theta");

    let trace2 = Scatter::new(t_positions, theta2_positions)
        .mode(Mode::Lines)
        .line(Line::new().width(1.0).color("#ff8c00ff").dash(plotly::common::DashType::Dash)) 
        .name("Hyperion dtheta");

    let mut plot = Plot::new();
    plot.add_trace(trace1);
    plot.add_trace(trace2);
    
    // 设置图表布局
    plot.set_layout(
        plotly::Layout::new()
            .title("Hyperion Position".into())
            .x_axis(plotly::layout::Axis::new().title("t".into()))
            .y_axis(plotly::layout::Axis::new().title("omega".into()))
    );
    
    // 保存为HTML文件
    let output_path = "src/plot/hyperion_double_omega_plot.html";
    plot.write_html(output_path);
    
    println!("✅ 图表已保存到 {}", output_path);
    Ok(())
}

pub fn plot_hyperion_lyapunov_e(results: &Vec<(f64, f64)>) -> Result<(), Box<dyn std::error::Error>> {
    // 拆分出 e 和 λ 向量，并对 λ 取绝对值
    let e_values: Vec<f64> = results.iter().map(|x| x.0).collect();
    let lambda_values: Vec<f64> = results.iter().map(|x| x.1.abs()).collect(); // ✅ 取绝对值

    // 创建折线图
    let trace = Scatter::new(e_values, lambda_values)
        .name("Lyapunov 指数")
        .mode(plotly::common::Mode::LinesMarkers);

    let mut plot = Plot::new();
    plot.add_trace(trace);

    // 设置标题与轴标签
    plot.set_layout(
        plotly::Layout::new()
            .title("Lyapunov 指数随偏心率变化图".into())
            .x_axis(plotly::layout::Axis::new().title("偏心率 e".into()))
            .y_axis(plotly::layout::Axis::new().title("Lyapunov 指数 |λ|".into())), // ✅ 标注绝对值
    );

    // 保存为HTML文件
    let output_path = "src/plot/hyperion_Lyapunov_e_scan.html";
    plot.write_html(output_path);
    println!("✅ 图表已保存到 {}", output_path);

    Ok(())
}
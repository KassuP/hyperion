import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# === 参数设置 ===
base_dir = os.path.dirname(__file__)
filename = os.path.join(base_dir, "data", "hyperion.csv")

l = 1.0e7     # 偏移距离，可按需要修改
speed_factor = 10000000000.0  # 动画加速倍数，越大越快（推荐 2~10）

# === 1. 读取 CSV 文件 ===
data = pd.read_csv(filename)

x = data["x"].values
y = data["y"].values
theta = data["theta"].values

# 确保角度为弧度制
if np.abs(theta).max() > 2 * np.pi:
    theta = np.radians(theta)

# === 2. 计算新坐标 ===
x_new = x + l * np.cos(theta)
y_new = y + l * np.sin(theta)

# === 3. 设置绘图窗口 ===
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect("auto")  # 不强制比例为1:1
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_title("Particle Motion (adjustable speed)")

# 不自动缩放，显示全轨迹范围
margin = 0.05 * max(x.max() - x.min(), y.max() - y.min())
ax.set_xlim(x.min() - margin, x.max() + margin)
ax.set_ylim(y.min() - margin, y.max() + margin)

# === 4. 初始化绘图对象 ===
line, = ax.plot([], [], 'r-', lw=0.5, label='Trajectory')  # ✅ 线条变细
dot, = ax.plot([], [], 'bo', markersize=4, label='Particle')  # 小点更自然
ax.legend()

# === 5. 更新函数 ===
def update(frame):
    line.set_data(x_new[:frame], y_new[:frame])
    dot.set_data([x_new[frame]], [y_new[frame]])
    return line, dot

# === 6. 创建动画 ===
n_frames = len(x_new)
interval = 4 / speed_factor  # ✅ 加速控制：数值越小动画越快

ani = animation.FuncAnimation(
    fig,
    update,
    frames=n_frames,
    interval=interval,
    blit=True,
    repeat=True
)

plt.show()

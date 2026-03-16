import numpy as np
import matplotlib.pyplot as plt

def plot_closed_shape_from_file(filename, color='blue', label=None, zorder=1):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    # 绘制封闭图形
    plt.plot(x, y, color=color, marker='o', linestyle='-', label=label, zorder=zorder)
    plt.plot([x[-1], x[0]], [y[-1], y[0]], color=color, linestyle='-', zorder=zorder)

def plot_points_from_file(filename, color='red', label=None, size=50, zorder=5):
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    # 绘制散点图，控制点的大小
    plt.scatter(x, y, color=color, marker='o', label=label, s=size, zorder=zorder)

# 创建一个图形并设置尺寸
plt.figure(figsize=(10, 8))

# 绘制封闭形状
plot_closed_shape_from_file('output/FBZ_location.dat', color='blue', label='FBZ')
plot_closed_shape_from_file('output/SBZ_location.dat', color='green', label='SBZ')

# 绘制散点图，指定点的大小
plot_points_from_file('output/DSSF_Mspace_location.dat', color='red', label='disperse', size=20)  # 增加点的大小

# 设置图形标题和坐标标签
plt.title('Combined Plot of Shapes and Points')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
#去掉网格线
plt.grid(False)
plt.axis('equal')  # 保持 x 和 y 轴等比例
plt.legend()  # 显示图例
plt.show()

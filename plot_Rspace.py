import numpy as np
import matplotlib.pyplot as plt

# 从文件加载数据
data1 = np.loadtxt("output/old_Rspace_location.dat")
data2 = np.loadtxt("output/new_Rspace_location.dat")

# 提取ID和坐标
ids1 = data1[:, 0]  # 第一列为ID
x_coords1 = data1[:, 1]  # 第二列为X坐标
y_coords1 = data1[:, 2]  # 第三列为Y坐标

ids2 = data2[:, 0]  # 第一列为ID
x_coords2 = data2[:, 1]  # 第二列为X坐标
y_coords2 = data2[:, 2]  # 第三列为Y坐标



# 创建新图形
plt.figure(figsize=(6.3, 5.3), facecolor='white')

# 绘制第一组坐标
#plt.plot(x_coords1, y_coords1, '.', color='blue', markersize=15)

# 为第一组点标记ID
#for i in range(len(ids1)):
 #   plt.text(x_coords1[i]-0.05, y_coords1[i], str(int(ids1[i])), 
  #           verticalalignment='bottom', horizontalalignment='right', 
   #          color='blue', fontsize=10)
    
#verticalalignment='bottom'：文本的底部与指定的坐标点对齐。
#x_coords1[i]-0.1 使其更加偏左

# 绘制第二组坐标
plt.plot(x_coords2, y_coords2, 'o', markersize=5, color='red', markerfacecolor='red')

# 为第二组点标记ID
for i in range(len(ids2)):
    plt.text(x_coords2[i]+0.05, y_coords2[i], str(int(ids2[i])), 
             verticalalignment='bottom', horizontalalignment='left', 
             color='red', fontsize=10)
    


# 设置坐标轴等比例显示
plt.axis('equal')

# 添加标题和坐标轴标签
plt.xlabel('X')
plt.ylabel('Y')

# 显示图形
plt.show()

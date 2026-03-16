HPhi Input Generator for KJGG' Model (Expert Mode)

依赖
HPhi（3.5.2）：最终执行计算的量子模拟软件包。
C++ 编译器： C++17 标准

简介: 
本代码用于生成量子格点模型模拟软件包 HPhi 专家模式（Expert Mode）所需的全部输入文件。它专门针对沿[111]方向施加外磁场的 Kitaev-Heisenberg-Gamma-Gamma' (KJGG') 模型设计，能够处理任意自旋大小的系统。

主要功能
哈密顿量文件：生成两体相互作用项，涵盖 Kitaev、Heisenberg、Gamma、Gamma' 耦合, 还包括单粒子各向异性。

磁场文件：生成沿 [111] 方向的单体外磁场项。

格林函数文件：自动生成用于计算以下物理量的文件：
单体格林函数
二体格林函数（自旋关联函数）
三体、四体、六体格林函数（例如用于计算手征序参量、标量自旋手性, Wp 等）

我还提供了一些Python脚本:

plot_Rspace.py: 用来绘制格点晶格结构
plot_Mspace.py,plot_Mlines.py: 用来绘制布里渊区动量点
build_DSSF.py, build_DSSF_MPI.py: 计算动态自旋结构因子
Run.py: 批量传参，第一个参数是 '\phi'，第二个参数是 'h'，随后编译，运行。K = cos(\phi Pi), G = -sin(\phi Pi)，其余参数在代码内部设为0。




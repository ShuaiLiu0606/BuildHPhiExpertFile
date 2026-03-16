HPhi Input Generator for KJGG' Model (Expert Mode)

依赖
1. HPhi (3.5.2)：量子格点模型模拟的核心计算引擎
2. C++ 编译器：支持 C++17 标准
3. MKL

项目简介
本项目提供一个专用的输入文件生成器，用于为量子模拟软件包 HPhi 的专家模式（Expert Mode）生成全部必需的输入文件。代码专门针对沿 [111] 晶体学方向施加外磁场的 Kitaev-Heisenberg-Gamma-Gamma' (KJGG') 模型 设计，支持任意自旋大小的系统模拟。

默认配置下，系统自旋设为 S=1/2。如需计算大自旋系统，可在 spinModel 类的构造函数中修改 SpinVal 参数值。

模型参数定义如下：
K=cos(ϕπ), G=−sin(ϕπ), h=h（外磁场强度）, 其余耦合强度均设为 0

主要功能
1. 哈密顿量生成：自动构建两体相互作用项，包括 Kitaev、Heisenberg、Gamma 和 Gamma' 耦合，以及单离子各向异性项。
2. 磁场项生成：生成沿 [111] 方向的单体外磁场项。
3. 单体格林函数
4. 二体格林函数（自旋关联函数）
5. 三体、四体及六体格林函数（用于计算手征序参量、标量自旋手性、Wp 等拓扑量）

附带的 Python 脚本
plot_Rspace.py：绘制实空间晶格结构。

plot_Mspace.py / plot_Mlines.py：绘制布里渊区动量点分布及高对称路径。

build_DSSF.py / build_DSSF_MPI.py：计算动态自旋结构因子（支持串行与 MPI 并行模式）。

Run.py：批量提交任务脚本，接受两个命令行参数：\phi 和 h。脚本将自动编译并运行模拟，其中 
K=cos(ϕπ)，G=−sin(ϕπ)，其余参数在代码内固定为 0。

使用范例 1：手动编译与运行
# 1. 配置构建
cmake -S . -B build

# 2. 编译项目
cmake --build build

# 3. 创建输出目录
mkdir output

# 4. 进入输出目录并运行程序（示例参数：phi=0.0, h=0.0）
cd output/

../build/my_pro 0.0 0.0

使用范例 2：使用 Python 脚本批量提交
# 配置构建,编译项目...
python3 Run.py
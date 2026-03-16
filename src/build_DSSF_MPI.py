import pandas as pd
import numpy as np
import subprocess
import math
import sys
import os

PI = np.pi  # 圆周率

# 确定调用 HPhi 时使用的 MPI 进程数，优先使用总任务数 SLURM_NTASKS，否则使用每节点任务数 SLURM_NTASKS_PER_NODE
if 'SLURM_NTASKS' in os.environ:
    nprocs = int(os.environ['SLURM_NTASKS'])
    print(f"检测到环境变量 SLURM_NTASKS，设置 MPI 进程数为 {nprocs}。", flush=True)
elif 'SLURM_NTASKS_PER_NODE' in os.environ:
    nprocs = int(os.environ['SLURM_NTASKS_PER_NODE'])
    print(f"未检测到 SLURM_NTASKS，使用每节点任务数 SLURM_NTASKS_PER_NODE，设置 MPI 进程数为 {nprocs}。", flush=True)
else:
    sys.stderr.write("错误：未设置 SLURM_NTASKS 或 SLURM_NTASKS_PER_NODE，无法确定 MPI 进程数。\n")
    sys.exit(1)


# 获取参数：L 格点数，和 HPhi 可执行文件路径
if len(sys.argv) == 3:
    L = int(sys.argv[1])
    path_to_HPhi = os.path.abspath(sys.argv[2])
else:
    print("错误\n用法: python build_DSSF.py L path_to_HPhi")
    sys.exit(1)

# 封装 MPI 调用 HPhi 的函数
def run_hphi():
    cmd = ['mpirun', '-np', str(nprocs), path_to_HPhi, '-e', 'namelist.def']
    return subprocess.call(cmd)

# 确保输出目录存在
os.makedirs('output', exist_ok=True)

# 1. 缺少 zvo_energy.dat 时，首次调用 HPhi 计算本征值与本征向量
if not os.path.exists('output/zvo_energy.dat'):
    run_hphi()

# 2. 读取基态能量
energy = 0.0
with open('output/zvo_energy.dat', 'r') as f:
    for line in f:
        parts = line.split()
        if parts and parts[0] == 'Energy':
            energy = float(parts[1])
            break

# 3. 修改 calcmod.def 文件，将 CalcSpec 设为 1
with open('calcmod.def', 'r') as f:
    lines = f.readlines()
for i, line in enumerate(lines):
    if 'CalcSpec' in line:
        parts = line.split()
        parts[-1] = '1'
        lines[i] = ' '.join(parts) + '\n'
with open('calcmod.def', 'w') as f:
    f.writelines(lines)

# 4. 确保 namelist.def 中包含 PairExcitation 行
with open('namelist.def', 'r+') as f:
    content = f.read()
    if 'PairExcitation' not in content:
        f.write('PairExcitation  pair.def\n')

# 5. 更新 modpara.def 参数：LanczosEps、NOmega、OmegaMin、OmegaMax，并添加 OmegaIm
WMin = 0.0 
WMax = 10.01

with open('modpara.def', 'r') as f:
    lines = f.readlines()
for i, line in enumerate(lines):
    parts = line.split()
    if 'LanczosEps' in line:
        parts[-1] = '14'
    elif 'NOmega' in line:
        parts[-1] = '1001'
    elif 'OmegaMax' in line:
        parts[-2] = str(WMax)
    elif 'OmegaMin' in line:
        parts[-2] = str(WMin)
    elif 'OmegaOrg' in line:
        parts[-2] = str(energy)
    lines[i] = '  '.join(parts) + '\n'
with open('modpara.def', 'w') as f:
    f.writelines(lines)
with open('modpara.def', 'a') as f:
    f.write('OmegaIm 0.05\n')

# 6. 读取动量空间 (Mspace) 和实空间 (Rspace) 坐标
Mspace_list = np.loadtxt('DSSF_Mspace_location.dat')

#Mspace_list = np.loadtxt('DSSF_Mline_location.dat')

Rspace_list = np.loadtxt('new_Rspace_location.dat')


# 7. 构建光谱函数，根据轴向 'Sx'、'Sy'、'Sz'
def build_spectrum(axis):
    num = 0
    for qx, qy in Mspace_list:
        # 写入 pair.def 文件
        with open('pair.def', 'w') as f:
            header = f"""===============================
NCisAitCjtAjs      {L*2}
===============================
====== PairExcitation =======
===============================
"""
            f.write(header)
            for j in range(L):
                Rx, Ry = Rspace_list[j, 1], Rspace_list[j, 2]
                phase = qx * Rx + qy * Ry
                wr, wi = math.cos(phase), math.sin(phase)
                if axis == 'Sx':
                    f.write(f"{j} 0 {j} 1 1 {0.5*wr} {0.5*wi}\n")
                    f.write(f"{j} 1 {j} 0 1 {0.5*wr} {0.5*wi}\n")
                elif axis == 'Sy':
                    f.write(f"{j} 0 {j} 1 1 {0.5*wi} {-0.5*wr}\n")
                    f.write(f"{j} 1 {j} 0 1 {-0.5*wi} {0.5*wr}\n")
                elif axis == 'Sz':
                    f.write(f"{j} 0 {j} 0 1 {0.5*wr} {0.5*wi}\n")
                    f.write(f"{j} 1 {j} 1 1 {-0.5*wr} {-0.5*wi}\n")
                    
        print(f"计算 {axis} 方向第 {num} 个 q 点: ({qx:.4f}, {qy:.4f})")
        # 调用 HPhi 进行并行计算
        run_hphi()
        # 复制结果文件
        subprocess.call(['cp', 'output/zvo_DynamicalGreen.dat', f'spectrum_{axis}{num}.dat'])
        num += 1
    # 合并所有 q 点的光谱数据
    with open(f'spectrum_{axis}.dat', 'w') as f_out:
        for i in range(num):
            qx, qy = Mspace_list[i]
            with open(f'spectrum_{axis}{i}.dat') as f_in:
                for line in f_in:
                    parts = line.split()
                    f_out.write(f"{i} {qx:.8f} {qy:.8f} {float(parts[0])} {-float(parts[3]) /PI:.16f}\n")
            f_out.write('\n')
    # 清理临时文件，并保存最终结果
    for i in range(num):
        os.remove(f'spectrum_{axis}{i}.dat')
    subprocess.call(['cp', f'spectrum_{axis}.dat', 'output/spectrum_{axis}.dat'])

# 分别构建 Sx、Sy、Sz 三个方向的光谱
for axis in ['Sx', 'Sy', 'Sz']:
    build_spectrum(axis)

print("所有光谱计算完成。")

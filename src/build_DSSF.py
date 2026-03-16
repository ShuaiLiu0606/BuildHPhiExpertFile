import pandas as pd
import numpy as np
import subprocess
import math
import sys
import os

PI = np.pi

# Set L and path(command) to HPhi
if len(sys.argv) == 3:
    L = int(sys.argv[1])  # 获取第一个命令行参数作为 L 的值
    path_to_HPhi = os.path.abspath(sys.argv[2])  # 第二个参数作为 HPhi 的路径
else:
    print("Error")
    print("Usage: python build_DSSF.py L path_to_HPhi.")
    print("L: integer value representing the lattice size.")
    print("path_to_HPhi: relative or absolute path to HPhi.")
    exit(1)

# first run HPhi for Eigenvalues and Eigenvectors
if not os.path.exists('./output/zvo_energy.dat'):
    subprocess.call([path_to_HPhi, '-e', 'namelist.def'])

# get bground Energy
energy = 0.0
with open('./output/zvo_energy.dat', 'r') as f:
    for line in f:
        line_array = line.split()
        if line_array[0] == 'Energy':
            energy = float(line_array[1])
            break

with open('calcmod.def', 'r') as file:
    lines = file.readlines()

for i, line in enumerate(lines):
    if 'CalcSpec' in line:
        parts = line.split()
        parts[-1] = '     1'
        lines[i] = ' '.join(parts) + '\n'

with open('calcmod.def', 'w') as file:  
    file.writelines(lines)

with open('namelist.def', 'r+') as f:
    lines = f.readlines()
    if not any('PairExcitation' in line for line in lines):
        f.write('PairExcitation  pair.def\n')



WMin = 0.0
WMax = 1.01
frequencyMin = WMin + energy
frequencyMax = WMax + energy

with open('modpara.def', 'r') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if 'Lanczos_max' in line:
        parts = line.split()
        parts[-1] = '     500'
        lines[i] = ' '.join(parts) + '\n'
    elif 'LanczosEps' in line:
        parts = line.split()
        parts[-1] = '     12'
        lines[i] = ' '.join(parts) + '\n'
    elif 'NOmega' in line:
        parts = line.split()
        parts[-1] = '     101'
        lines[i] = ' '.join(parts) + '\n'
    elif 'OmegaMax' in line:
        parts = line.split()
        parts[-2] = str(frequencyMax)
        parts[-1] = '     0'
        lines[i] = ' '.join(parts) + '\n'
    elif 'OmegaMin' in line:
        parts = line.split()
        parts[-2] = str(frequencyMin)
        parts[-1] = '     0'
        lines[i] = ' '.join(parts) + '\n'

with open('modpara.def', 'w') as file:
    file.writelines(lines)
    
with open('modpara.def', 'a') as file:
    file.write('OmegaIm 0.05\n')


#***************Obtain the momentum point and lattice coordinates.
 
Mon_file = 'DSSF_Mline_location.dat'
#Mon_file = 'DSSF_Mspace_location.dat'  
Mspace_list = np.loadtxt(Mon_file)


Rspace_file = 'new_Rspace_location.dat'
Rspace_list = np.loadtxt(Rspace_file)


#********************************DSSF_Sx****************************
num_qx = 0
for row in range(Mspace_list.shape[0]):
    qx = Mspace_list[row, 0]
    qy = Mspace_list[row, 1]

    with open('pair.def', 'w') as f:
        pair_text = """===============================
NCisAitCjtAjs      %i
===============================
====== PairExcitation =======
===============================
""" % (L * 2)
                
        f.write(pair_text)

        for j in range(L):

            Rx = Rspace_list[j, 1]    
            Ry = Rspace_list[j, 2]   
            
            phase = qx * Rx + qy * Ry

            wr = math.cos(phase)
            wi = math.sin(phase)

            f.write(f'{j} 0 {j} 1 1 {0.5*wr} {0.5*wi}\n')
            f.write(f'{j} 1 {j} 0 1 {0.5*wr} {0.5*wi}\n')

    print(f'qoint: {qx,qy}\n')
    subprocess.call([path_to_HPhi, '-e', 'namelist.def'])
    subprocess.call(['cp', 'output/zvo_DynamicalGreen.dat', f'spectrum_Sx{num_qx}.dat'])
    num_qx = num_qx + 1

with open('spectrum_Sx.dat', 'w') as f:
    for i in range(num_qx):
        qx = Mspace_list[i, 0]
        qy = Mspace_list[i, 1]
        with open(f'spectrum_Sx{i}.dat', 'r') as g:
            for line in g:
                line_array = line.split()
                f.write(f'{int(i)} {qx:.8f} {qy:.8f} {float(line_array[0]) - energy:.8f} {-float(line_array[3]) / PI:.16f}\n')
        f.write('\n')

subprocess.call(["mv", "spectrum_Sx.dat", "spectrum_Sx.dat.bak"])
subprocess.call(["rm spectrum_Sx*.dat"], shell=True)
subprocess.call(["mv", "spectrum_Sx.dat.bak", "spectrum_Sx.dat"])
subprocess.call(["cp", "spectrum_Sx.dat", "output/spectrum_Sx.dat"])

#********************************DSSF_Sy****************************
num_qy = 0
for row in range(Mspace_list.shape[0]):
    qx = Mspace_list[row, 0]
    qy = Mspace_list[row, 1]

    with open('pair.def', 'w') as f:
        pair_text = """===============================
NCisAitCjtAjs      %i
===============================
====== PairExcitation =======
===============================
""" % (L * 2)
                
        f.write(pair_text)

        for j in range(L):
            Rx = Rspace_list[j, 1]    
            Ry = Rspace_list[j, 2]   
            
            phase = qx * Rx + qy * Ry
            wr = math.cos(phase)
            wi = math.sin(phase)

            f.write(f'{j} 0 {j} 1 1 {wi*0.5} {-wr*0.5}\n')
            f.write(f'{j} 1 {j} 0 1 {-wi*0.5} {wr*0.5}\n')

    print(f'qoint: {qx,qy}\n')
    subprocess.call([path_to_HPhi, '-e', 'namelist.def'])
    subprocess.call(['cp', 'output/zvo_DynamicalGreen.dat', f'spectrum_Sy{num_qy}.dat'])
    num_qy = num_qy + 1

with open('spectrum_Sy.dat', 'w') as f:
    for i in range(num_qy):
        qx = Mspace_list[i, 0]
        qy = Mspace_list[i, 1]
        with open(f'spectrum_Sy{i}.dat', 'r') as g:
            for line in g:
                line_array = line.split()
                f.write(f'{int(i)} {qx:.8f} {qy:.8f} {float(line_array[0]) - energy:.8f} {-float(line_array[3]) / PI:.16f}\n')
        f.write('\n')

subprocess.call(["mv", "spectrum_Sy.dat", "spectrum_Sy.dat.bak"])
subprocess.call(["rm spectrum_Sy*.dat"], shell=True)
subprocess.call(["mv", "spectrum_Sy.dat.bak", "spectrum_Sy.dat"])
subprocess.call(["cp", "spectrum_Sy.dat", "output/spectrum_Sy.dat"])

#********************************DSSF_Sz****************************
num_qz = 0
for row in range(Mspace_list.shape[0]):
    qx = Mspace_list[row, 0]
    qy = Mspace_list[row, 1]

    with open('pair.def', 'w') as f:
        pair_text = """===============================
NCisAitCjtAjs      %i
===============================
====== PairExcitation =======
===============================
""" % (L * 2)
                
        f.write(pair_text)

        for j in range(L):

            Rx = Rspace_list[j, 1]    
            Ry = Rspace_list[j, 2]   
            
            phase = qx * Rx + qy * Ry

            wr = math.cos(phase)
            wi = math.sin(phase)

            f.write(f'{j} 0 {j} 0 1 {0.5*wr} {0.5*wi}\n')
            f.write(f'{j} 1 {j} 1 1 {-0.5*wr} {-0.5*wi}\n')

    print(f'qoint: {qx,qy}\n')
    subprocess.call([path_to_HPhi, '-e', 'namelist.def'])
    subprocess.call(['cp', 'output/zvo_DynamicalGreen.dat', f'spectrum_Sz{num_qz}.dat'])
    num_qz = num_qz + 1

with open('spectrum_Sz.dat', 'w') as f:
    for i in range(num_qz):
        qx = Mspace_list[i, 0]
        qy = Mspace_list[i, 1]
        with open(f'spectrum_Sz{i}.dat', 'r') as g:
            for line in g:
                line_array = line.split()
                f.write(f'{int(i)} {qx:.8f} {qy:.8f} {float(line_array[0]) - energy:.8f} {-float(line_array[3]) / PI:.16f}\n')
        f.write('\n')

subprocess.call(["mv", "spectrum_Sz.dat", "spectrum_Sz.dat.bak"])
subprocess.call(["rm spectrum_Sz*.dat"], shell=True)
subprocess.call(["mv", "spectrum_Sz.dat.bak", "spectrum_Sz.dat"])
subprocess.call(["cp", "spectrum_Sz.dat", "output/spectrum_Sz.dat"])

#

file_names = ['spectrum_Sx.dat', 'spectrum_Sy.dat', 'spectrum_Sz.dat']
output_file = 'spectrum.dat'

# 初始化用于存储第五列求和的列表
sum_column_5 = []

# 读取每个文件并求和
for file_name in file_names:
    with open(file_name, 'r') as f:
        for i, line in enumerate(f):
            stripped_line = line.strip()
            if not stripped_line:
                continue
            
            values = stripped_line.split()
            
            # 确保每行有至少5列
            if len(values) >= 5:
                if len(sum_column_5) <= i:
                    sum_column_5.append(float(values[4]))  # 第五列在索引为4的位置
                else:
                    sum_column_5[i] += float(values[4])  # 累加第五列的值

# 创建输出文件并写入新数据
with open(file_names[0], 'r') as f_in, open(output_file, 'w') as f_out:
    for i, line in enumerate(f_in):
        stripped_line = line.strip()
        if not stripped_line:
            continue
            
        values = stripped_line.split()
        
        # 确保每行有至少5列
        if len(values) >= 5:
            # 输出格式: 整数值, 前四列值, 第五列求和
            new_line = f"{int(values[0])} {float(values[1]):.8f} {float(values[2]):.8f} {float(values[3]):.8f} {sum_column_5[i]:.16f}\n"
            f_out.write(new_line)

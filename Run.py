import os
import shutil
import subprocess
import numpy as np

# 配置参数
source_dir = os.path.dirname(os.path.abspath(__file__))  # 获取当前脚本所在目录
build_dir = os.path.join(source_dir, "build")            # 构建目录
new_output_dir_pattern = "KG_theta{param1:.3f}PI_h{param2:.3f}"  # 输出文件夹的新命名格式

#param2_values = np.append(np.arange(0.64, 0.80, 0.02), 0.80) 
results_dir = "KG_theta0.000"     # 汇总文件夹
param1_values = [0.00]
param2_values = [0.00]

def build_once():
    """运行 cmake 和 make, 并执行生成的可执行文件"""
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    print("配置 CMake...")
    # 使用正确的 CMake 配置方式
    subprocess.run(["cmake", "-S", ".", "-B", "build"], cwd=source_dir, check=True)

    print("构建项目...")
    subprocess.run(["cmake", "--build", "build"], cwd=source_dir, check=True)
    
    # 检查可执行文件是否生成
    executable_path = os.path.join(build_dir, "my_pro")
    if not os.path.exists(executable_path):
        raise FileNotFoundError(f"编译失败：找不到可执行文件 {executable_path}")
    
    print("构建成功！")

def run_program(param1, param2, output_path):
    """在指定目录下运行 build/my_pro 并保存日志"""
    executable_path = os.path.join(build_dir, "my_pro")
    if not os.path.exists(executable_path):
        raise FileNotFoundError(f"找不到可执行文件 {executable_path}")

    print(f"Running ./my_pro theta={param1:.3f}, h={param2:.3f} ...")
    log_file = os.path.join(output_path, f"log_theta{param1:.3f}_h{param2:.3f}.dat")

    with open(log_file, "w") as f:
        subprocess.run([executable_path, str(param1), str(param2)], cwd=output_path,
            stdout=f, stderr=subprocess.STDOUT, check=True)

# 主程序
success_count = 0
total_count = len(param1_values) * len(param2_values)

# 确保 results 目录存在
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# 编译一次
try:
    build_once()
except Exception as e:
    print(f"编译失败: {e}")
    exit(1)

# 保存参数文件（可选）
para_file = os.path.join(results_dir, "Para.dat")
with open(para_file, "w") as f:
    f.write("# theta     h\n")
    for param1 in param1_values:
        for param2 in param2_values:
            f.write(f"{param1:.6f}  {param2:.6f}\n")

for i, param1 in enumerate(param1_values):
    for j, param2 in enumerate(param2_values):
        print(f"\n=== 处理参数: theta={param1:.3f}, Field={param2:.3f} "
              f"({i*len(param2_values)+j+1}/{total_count}) ===")
        try:
            # 为每个参数组合创建独立的输出目录
            output_name = new_output_dir_pattern.format(param1=param1, param2=param2)
            current_output_dir = os.path.join(results_dir, output_name)
            
            # 清空或创建目录
            if os.path.exists(current_output_dir):
                shutil.rmtree(current_output_dir)
            os.makedirs(current_output_dir)
            
            # 运行程序
            run_program(param1, param2, current_output_dir)
            print(f"✅ 运行完成: theta={param1:.3f}, h={param2:.3f}")
            success_count += 1
            
        except subprocess.CalledProcessError:
            print(f"❌ 运行失败：theta={param1:.3f}, h={param2:.3f}")
        except Exception as e:
            print(f"⚠️ 未知错误: {e}")

print(f"\n=== 全部完成 ===")
print(f"成功处理: {success_count}/{total_count}")
print(f"所有结果已保存在 '{results_dir}' 目录.")

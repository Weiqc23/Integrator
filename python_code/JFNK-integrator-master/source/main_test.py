
import time
import numpy as np
import matplotlib.pyplot as plt
from integrator import jfnk_adaptive
from params import Params
import points
import os
import numpy as np
from scipy.sparse import eye, csr_matrix
from scipy.sparse.linalg import spsolve

from analysis import analyze_corrections, error_analysis, get_error_time_nodes, get_error_time_steps
from plotter import corrections, plot_correction_time_steps, plot_errors
# ODE 的导数函数
# def f_eval(t, y):
#     return -y

start = time.time()

def load_matrices_from_path(A_path, B_path=None):
    """
    从指定路径加载A、B矩阵 - 无pandas版本
    """
    global A_matrix, B_matrix
    
    # 读取A矩阵
    print(f"正在读取A矩阵: {A_path}")
    if os.path.exists(A_path):
        try:
            # 使用numpy读取CSV文件
            A_data = np.loadtxt(A_path, delimiter=',')
            A_matrix = csr_matrix(A_data)
            print(f"A矩阵加载成功: 形状{A_matrix.shape}, 数据类型{A_matrix.dtype}")
        except Exception as e:
            print(f"读取A矩阵时出错: {e}")
            raise
    else:
        raise FileNotFoundError(f"A矩阵文件不存在: {A_path}")
    
    # 读取B矩阵（如果提供）
    if B_path:
        print(f"正在读取B矩阵: {B_path}")
        if os.path.exists(B_path):
            try:
                B_data = np.loadtxt(B_path, delimiter=',')
                B_matrix = csr_matrix(B_data)
                print(f"B矩阵加载成功: 形状{B_matrix.shape}, 数据类型{B_matrix.dtype}")
            except Exception as e:
                print(f"读取B矩阵时出错: {e}")
                # 不抛出异常，允许B矩阵加载失败
        else:
            print(f"警告: B矩阵文件不存在: {B_path}")

def load_csv(filename):
    return np.loadtxt(filename, delimiter=",")


# 节点数
n_nodes = 4

# 加载系统矩阵和初始条件
A_file_path = "C:\\Users\\wei13\\Desktop\\Basic_Code\\Basic_Code\\Ex_AC_System\\problem_dependent\\A.csv"
B_file_path = "C:\\Users\\wei13\\Desktop\\Basic_Code\\Basic_Code\\Ex_AC_System\\problem_dependent\\B.csv"
load_matrices_from_path(A_file_path, B_file_path)
global A_matrix, B_matrix
y0 = load_csv("C:\\Users\\wei13\\Desktop\\Basic_Code\\Basic_Code\\Ex_AC_System\\problem_dependent\\y0.csv")
y0 = y0.reshape(-1)
# print("A矩阵信息： ")
# print(A)
# print("B矩阵信息： ")
# print(B)
# print("y0初始条件信息： ")
# print(y0)

m = A_matrix.shape[0]


def u_func(t):

    peak_amplitude = 735e3 / np.sqrt(3) * np.sqrt(2)
    freq = 60.0
    omega = 2 * np.pi * freq
    phase = 0.0
    u = np.zeros((1, 1))
    u[0, 0] = peak_amplitude * np.sin(omega * t + phase)

    return u

def f_eval(t, y):
    
    y = y.reshape(-1, 1)
    # print("A shape: ", A.shape)
    # print("y shape: ", y.shape)
    # print("B shape: ", B.shape)
    # print("u_func(t) shape: ", u_func(t).shape)
    return (A_matrix @ y + B_matrix.transpose() * u_func(t)).reshape(-1)

node_type = points.GAUSS_RADAU_2A
params = Params(n_nodes=n_nodes, m=m, node_type=node_type)
# print("GAUSS_LEGENDRE")
# print (params.t)
# node_type2 = points.GAUSS_LOBATTO
# params2 = Params(n_nodes=n_nodes, m=m, node_type=node_type2)
# print("GAUSS_LOBATTO")
# print (params2.t)

# node_type3 = points.GAUSS_RADAU
# params3 = Params(n_nodes=n_nodes, m=m, node_type=node_type3)
# print("GAUSS_RADAU")
# print (params3.t)

node_type4 = points.GAUSS_RADAU_2A
params4 = Params(n_nodes=n_nodes, m=m, node_type=node_type4)
np.set_printoptions(precision=15, suppress=True)
print("GAUSS_RADAU_2A")
print("params.t: ", params4.t)
print("params.S: ", params4.S)



# print("params.n_nodes: ", params.n_nodes)
# print("params.node_type: ", params.node_type)

t_init = 0.0
t_final = 0.1

# 初始步长
dt_init = 100e-6

# 容差：控制数值解的误差
tol = 1e-6

# 保持使用从CSV加载的初始条件
# y0 = y0.reshape(-1, 1)  # 确保y0是列向量形式

print("JFNK begin\n")
# JFNK求解ODE
t_all, y_all, Y_all, D_all, h_all = jfnk_adaptive(f_eval, t_init, t_final, dt_init, params, y0,A_matrix, B_matrix, tol,do_print=True)
print("JFNK finish\n")

print("t_all信息： ")
print(t_all)

# print("y_all信息： ")
# print(y_all)
print("y_all.shape: ", y_all.shape)
# print("y_all.row(0)信息： ")
# print(y_all[0, :])

np.set_printoptions(precision=12, suppress=True)

print("y_all.row(row-1)信息： ")
print(y_all[-1, :])

print("y_all.row(row-2)信息： ")
print(y_all[-2, :])


last_row = y_all[-1, :]
second_last_row = y_all[-2, :]

out_array = np.vstack([last_row, second_last_row])

np.savetxt('y_all_output.csv', out_array, fmt='%.12f', delimiter=',')


def plot_results(t_all, y_all):
    plt.rcParams["font.family"] = ["SimHei"]  # 黑体
    plt.rcParams["axes.unicode_minus"] = False  # 负号显示问题
    
    plt.figure(figsize=(10, 6))
    
    plt.plot(t_all, y_all[:, 0], 'b-', linewidth=2, label='第0维解')
    
    plt.title('JFNK_py')
    plt.xlabel('t_all')
    plt.ylabel('y_all(row0)')
    
    plt.legend()
    
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    
    plt.show()

end = time.time()
print("running time: ")
print(end - start, " s")
# 调用绘图函数
plot_results(t_all, y_all)


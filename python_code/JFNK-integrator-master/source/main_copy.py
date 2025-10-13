import numpy as np
import matplotlib.pyplot as plt
from integrator import jfnk_adaptive
from params import Params
import points
from analysis import analyze_corrections, error_analysis, get_error_time_nodes, get_error_time_steps
from plotter import corrections, plot_correction_time_steps, plot_errors
# ODE 的导数函数
def f_eval(t, y):
    return -y

# 节点数
n_nodes = 5

# 求解y的维度
m = 1

node_type = points.GAUSS_LOBATTO

params = Params(n_nodes=n_nodes, m=m, node_type=node_type)

t_init = 0.0
t_final = 1.0

# 初始步长
dt_init = 0.1

# 容差：控制数值解的误差
tol = 1e-6

# 初始条件
y0 = np.array([1.0])

print("JFNK begin\n")
# JFNK求解ODE
t_all, y_all, Y_all, D_all, h_all = jfnk_adaptive(f_eval, t_init, t_final, dt_init, params, y0, tol)
print("JFNK finish\n")

# 准确解
y_exact = np.exp(-t_all)[:, None]

print("\nt:")
print(t_all)

print("Numerical Solution (y_all):")
print(y_all)

print("\nExact Solution (y_exact):")
print(y_exact)

# 分析校正范数
d_norm, log_d_norm, d_norm_rel, log_d_norm_rel = analyze_corrections(D_all, Y_all)
print("Absolute Correction Norms:\n", d_norm)
print("Log10 Absolute Correction Norms:\n", log_d_norm)
print("Relative Correction Norms:\n", d_norm_rel)
print("Log10 Relative Correction Norms:\n", log_d_norm_rel)

# 计算误差
aerr, rerr = error_analysis(y_all, y_exact)
print("Absolute Error:\n", aerr)
print("Relative Error:\n", rerr)

# 计算每个时间节点的误差
y_spline = lambda t: np.exp(-t)
err_nodes = get_error_time_nodes(t_all, y_all, y_spline, do_relative=True)
print("Relative Error at Time Nodes:\n", err_nodes)

# 计算每个时间步的误差
err_steps = get_error_time_steps(t_all, y_all, y_spline, n_nodes=params.n_nodes, do_relative=True)
print("Relative Error at Time Steps:\n", err_steps)

# 绘制数值解和精确解
plt.plot(t_all, y_all, label='Numerical Solution', marker='o')
plt.plot(t_all, y_exact, label='Exact Solution', linestyle='--')
plt.xlabel('Time')
plt.ylabel('Solution')
plt.legend()
plt.title('Numerical Solution vs Exact Solution')
plt.show()

# 使用 plotter.py 中的绘图函数绘制数值解和精确解
plot_errors([y_all, y_exact], ["Numerical Solution", "Exact Solution"], title="Numerical Solution vs Exact Solution", xlabel="Time", ylabel="Solution")

# 使用 plotter.py 中的绘图函数绘制误差分析图
plot_errors([aerr, rerr], ["Absolute Error", "Relative Error"], title="Error Analysis", xlabel="Time", ylabel="Error")

# 使用 plotter.py 中的绘图函数绘制校正分析图
corrections(d_norm, log_d_norm, np.linalg.norm(y_exact, axis=0), do_rerr=False, labels=["Absolute Correction Norm", "Log10 Absolute Correction Norm"], do_legend=True)
# This file was written by Dr. Namdi Brandon
# ORCID: 0000-0001-7050-1538
# June 29, 2018

"""
This module contains code for functions for solving temporal ordinary
differential equations (ODEs)

.. math::
    \\frac{dy(t)}{dt} &= f(t, y(t)) \\\\
    y(0) &= y_0

which has the following solution

.. math::
    y(t) = y_0 +  \int^t_0 \\, f(\\tau, y(\\tau) ) \\, \\mathrm{d}\\tau

The numerical integrators are the following:

* Backward (implicit) Euler method
* Jacobian-Free Newton-Krylov (JFNK) method with uniform time stepping
* Jacobian-Free Newton-Krylov (JFNK) method with adaptive time stepping
* Spectral Deferred Corrections (SDC) (implicit)
* Spectral solution (Gauss collocation formulation) solver

+---------------+-----------------------------------+
| Abbreviations | Meaning                           |
+===============+===================================+
| JFNK          | Jacobian-Free Newton Krylov       |
+---------------+-----------------------------------+
| ODE           | Ordinary differential equation    |
+---------------+-----------------------------------+
| SDC           | Spectral deferred corrections     |
+---------------+-----------------------------------+

.. moduleauthor:: Dr. Namdi Brandon
"""

# ===============================================
# import
# ===============================================
# 导入NumPy库，用于数值计算
import numpy as np
# 从NumPy库中导入线性代数模块
import numpy.linalg as LA

# 从SciPy库的优化模块导入牛顿-克里洛夫求解器
from scipy.optimize import newton_krylov
# 从SciPy库的优化模块导入安德森加速算法
from scipy.optimize import anderson
# 从SciPy库的优化模块导入通用根查找函数
from scipy.optimize import root

# 导入time模块，用于计时
import time

# 导入points模块，该模块来自pyweno代码（作者为Matt Emmett）
# 此模块可能包含节点生成相关的功能
import points

# ===============================================
# constants
# ===============================================

# 牛顿-克里洛夫求解器在残差中的相对容差
NK_FTOL = 1e-1

# 向后欧拉方法的相对容差
BE_TOL = 1e-12

# 谱延迟校正（SDC）方法中校正强度的相对收敛准则
SDC_TOL = 1e-14

# 最大的SDC迭代次数
N_ITER_MAX_SDC = 500

# 最大的JFNK迭代次数
N_ITER_MAX_NEWTON = 50

# 用于调试，自适应步长的最大步数
N_STEPS_MAX_ADAPTIVE = int(1e9)

# 时间步长能增长的最大倍数
ADAPTIVE_SCALER_MAX = 4

# 时间步长能增长的最小倍数
ADAPTIVE_SCALER_MIN = 1.5

# ===============================================
# functions
# ===============================================

def adjust_scaler(x, x_min=ADAPTIVE_SCALER_MIN, x_max=ADAPTIVE_SCALER_MAX):

    """
    For the adaptive time stepping algorithm, this function adjusts the adaptive step size

    .. math::
        x \leftarrow
        \\begin{cases}
            x & \\text{if } x < 1 \\\\
            1 & \\text{if } 1 \\le x \\le x_{min} \\\\
            \\min(x, x_{max}) & \\text{if } x > x_{min}
        \\end{cases}

    :param float x: the ratio in which to grow or shrink the time step
    :param float x_min: the minimum ratio in which to grow a time step
    :param float x_max: the maximum ratio in which to grow a time step

    :return: the ratio in which to grow or shrink the time step
    :rtype: float
    """

    if x > x_min:
        # 上限处理，最大缩放倍数为x_max，避免步长过大
        x = min(x, x_max)

    elif (x >= 1) and (x <= x_min):
        x = 1.0

    return x

#求解dy/dt=f(y,t),y(t0)=y0
def backward_euler(f_eval, t, y0, S, be_tol=BE_TOL, do_print=False):

    """
    This function runs the backward Euler method over all of the nodes :math:`t_i` in an entire time step
    of size :math:`\\Delta{t}`.

     .. math::
        y - \\Delta{t}\\tilde{S}F(y) = rhs

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes over a time step :math:`\\Delta{t}`
    :param numpy.ndarray y0: the initial condition ( length = m)
    :param numpy.ndarray S: the backward Euler integration matrix :math:`\\tilde{S}`
    :return:
    """

    # 开始计时
    start = time.time()

    # 时间节点的数量
    n_nodes = len(t)

    # 系统的规模
    m = len(y0)

    # 时间步内的解，初始化为零矩阵
    y = np.zeros((n_nodes, m))

    #
    # 对整个时间步运行向后欧拉方法
    #
    print
    for i in range(n_nodes):

        h = S[i, i]
        # 如果使用左端点（例如高斯-洛巴托节点）
        if i == 0 and h == 0:
            # print("y",y.shape)
            # print ("y0 ", y0.shape)
            y[0, :] = y0[:]
        else:
            # 如果不使用左端点
            if i == 0 and h != 0:
                rhs = y0[:]
            else:  # 其他所有点
                rhs = y[i - 1, :]

            # 对第i个节点求解向后欧拉方法
            y[i, :] = backward_euler_node(f_eval, t=t[i], h=h, rhs=rhs, x0=rhs, be_tol=be_tol)
        # print("y.row(i): ", i)
        # print(y[i, :])

    # 结束计时
    end = time.time()

    if do_print:
        # print_elapsed_time(start, end)
        print("back_euler runing time: ")
        print(end - start, " s")

    return y

def backward_euler_node(f_eval, t, h, rhs, x0, be_tol=BE_TOL):

    """
    This function solves an ODE system using backward Euler method at a specific time :math:`t`. This code uses
    a general numerical solver in order to do the inversion in the backward Euler method in order
    to find :math:`y_i`.

    .. math::
        y_i - h_i f(t_i, y_{i-1}) = rhs_i

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param float t: the time at a node
    :param float h: the step size
    :param numpy.ndarray rhs: the right hand side of the backward euler system [m x 1]
    :param numpy.ndarray x0: the initial guess for the solution [m x 1]
    :param float be_tol: the tolerance for the backward Euler solver

    :return: the approximate solution [m x 1]
    :rtype: numpy.ndarray
    """

    # 近似解，初始化为零向量
    y       = np.zeros( rhs.shape )

    if h == 0:
        y[:] = rhs[:]
    else:

        # 定义需要求根的系统
        A = lambda x: x - h*f_eval(t, x) - rhs

        # 使用通用根查找函数求解系统
        # print ("y ", y.shape)
        # print ("x ", root(A, x0, tol=be_tol).x.shape)
        # 开始计时
        start = time.time()
        y[:] = root(A, x0, tol=be_tol).x
        
        # print("y: ", y.shape)
        # print(y)
        # print("x0.size(): ", x0.size)
        # print("x0: ", x0)
        # 结束计时
        end = time.time()
        # print("root running: ", end - start, " s")
    
    # print("backward_euler_node")
    # print("t: ", t)
    # print("h: ", h)
    # print("rhs.size(): ", rhs.size)
    # print("rhs: ", rhs)
    # print("x0.size(): ", x0.size)
    # print("x0: ", x0)
    # print("y: ", y)
    
    return y

#判断迭代终止
def convergence_criteria(d, y, tol):

    """
    This function calculates whether or not the correction vector :math:`\\delta` is small enough to
    satisfy the convergence criteria.

    .. math::
        x = \\frac{ \\displaystyle{ \\| \\delta \\| } }{ \\displaystyle{ \\| y \\| } } \\\\
        \\begin{cases}
            x \\le tol & \\text{converged} \\\\
            x > tol & \\text{not converged} \\\\
        \\end{cases}

    :param numpy.ndarray d: the correction vector :math:`\\delta` for a given iteration dimensions (n nodes, size of problem)
    :param numpy.ndarray y: the approximate solution :math:`y` for a given iteration  dimensions (n nodes, size of problem)
    :param float tol: the correction tolerance for the convergence criteria

    :return: a flag indicating whether or not the corrections are small enough to qualify for convegence
    :rtype: bool
    """

    # 计算相对范数
    value = relative_norm(d, y)

    # 如果相对范数小于等于容差，则认为方法收敛
    is_converged = value <= tol

    return is_converged

def jfnk(f_eval, t, y0, y_approx, S, S_p, spectral_radius, n_iter_max_newton=N_ITER_MAX_NEWTON, \
             be_tol=BE_TOL, sdc_tol=SDC_TOL, do_print=False):

    """
    The Jacobian-Free Newton-Krylov (JFNK) method to approximate a solution to the spectral solution

     .. math::
        y - \\Delta{t}SF(y) = y_0

    over one time step of size :math:`\\Delta{t}`. This is done by use a modified version of Newton's method to
    find a calculate a solution

    .. math::
        H(y) = 0

    where :math:`H(y^{[k]}) = \\delta^{[k]}` corresponds to one iteration of the SDC method.

    Given :math:`y^{[0]}`, this method does the following

    1. calculate the initial SDC iterations

        .. math::

            \\begin{cases}
	            \\delta^{[k]} &= H(y^{[k]})  \\text{\indent calculate an SDC correction} \\\\
	            y^{[k+1]} &= y^{[k]} + \\delta^{[k]} \\text{\\indent update the SDC solution}
            \\end{cases}

        until the the solution converges or order convergence has been observed

    2. do the Newton (Jacobian-Free) iterations

        .. math::

            J_{H}(y^{[p]})\\Delta{x} &= -H(y^{[p]})   \\\\
            \\implies J_{H}(y^{[p]})\\Delta{x} &= -\\delta^{[p]}   \\\\

        Set :math:`\\Delta{x} = \\sum^{p-1}_{j=0} c_j \\delta^{[j]}` and solve

        .. math::
            J_{H}(y^{[p]})\\sum^{p-1}_{j=0} c_j \\delta^{[j]} &= -\\delta^{[p]}  \\\\
            \\implies \\sum^{p-1}_{j=0} c_j (\\delta^{[k+1]} - \\delta^{[k]}) &= -\\delta^{[p]} \\\\

        Solve the system for the Jacobian-Free system
            .. math::
                \\begin{cases}
                    Ac &= -\\delta^{[p]} \\\\
                    y &\\leftarrow y^{[p]} + \\sum^{p-1}_{j=0} c_j \\delta^{[j]}
                \\end{cases}

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes over a time step :math:`\\Delta{t}`
    :param numpy.ndarray y0: the initial condition ( length = m)
    :param numpy.ndarray y_approx: the provisional solution (dimensions n time nodes, size of the problem)
    :param numpy.ndarray S: the spectral integration (Gaussian quadrature) matrix, :math:`S`
    :param numpy.ndarray S_p: the backward Euler integration matrix, :math:`\\tilde{S}`
    :param float spectral_radius: the spectral radius of the correction matrix for the extremely stiff case
    :param int n_iter_max_newton: the maximum number of Newton iterations
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param bool do_print: a flag indicating whether or not to print the elapsed time

    :return: the solution, the history of approximations for each iteration, the history of corrections \
    for each iteration

    :rtype: numpy.ndarray (dimensions n time nodes, size of the problem) , list (length number of iterations), \
    list (length of iterations), bool, bool, numpy.ndarray (length of iterations)
    """

    # 开始计时
    start = time.time()

    # 在刚度导致的阶收敛之前进行初始SDC迭代
    y_init, Y_init, D_init, is_converged, is_stiff, ratios \
        = jfnk_initial(f_eval, t, y0, y_approx, S, S_p, spectral_radius, be_tol=be_tol, sdc_tol=sdc_tol)
    # print("jfnk_initial")
    # print("t.size(): ", t.size)
    # print(t)
    # print("y0.size(): ", y0.size)
    # print(y0)
    # print("y_init.size(): ", y_init.size)
    # print(y_init)
    Y, D = list(), list()

    # 存储解
    if is_converged:
        y = y_init

    # 一旦检测到阶收敛，使用JFNK方法
    if (not is_converged) and (is_stiff):
        y, Y, D, is_converged = jfnk_iterations(f_eval, t, y0, y_init, S, S_p, n_iter_max_newton, \
                                     be_tol=be_tol, sdc_tol=sdc_tol)

    # 存储解的历史记录
    Y = [Y_init, Y]
    Y = [subitem for item in Y for subitem in item]

    # 存储校正的历史记录
    D = [D_init, D]
    D = [ subitem for item in D for subitem in item]

    # 停止计时
    end = time.time()

    if do_print:
        print("jfnk run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    # print("jfnk")
    # print("t.size(): ", t.size)
    # print(t)
    # print("y.size(): ", y.size)
    # print(y)
    return y, Y, D, is_converged, is_stiff, ratios

def jfnk_adaptive(f_eval, t_init, t_final, dt_init, p, y0, tol, n_iter_max_newton=N_ITER_MAX_NEWTON, be_tol=BE_TOL, \
                  sdc_tol=SDC_TOL, n_steps_max=N_STEPS_MAX_ADAPTIVE, do_print=False,):
    """
    Run the JFNK with adaptive step sizes from :math:`t \\in [t_{init}, t_{final}]` to calculate an
    approximation to the solution

    .. math::
        y(t_{final}) = y(t_{init}) +  \int^{t_{final}}_{t_{init}} \\, f(\\tau, y(\\tau) ) \\, \\mathrm{d}\\tau

    Such that for each time step the step size :math:`\\Delta{t}` is chosen so that the difference between
    the exact solution :math:`y` and the approximate solution :math:`\\tilde{y}`

    .. math::
        \||y - \\tilde{y} \\|_{\\infty} \\le tol

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param float t_init: the initial time :math:`t_{init}`
    :param float t_final: the final time :math:`t_{final}`
    :param float dt_init: the initial step size :math:`\\Delta{t}_{init}`
    :param params.Params p: the parameters related to the solver
    :param numpy.ndarray y0: the initial condition :math:`y(t_{init})` (length = m)
    :param float tol: the approximated absolute error at each step for the adaptive solution
    :param int n_iter_max_newton: the maximum number of Newton iterations
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param int n_steps_max: the maximum number of steps in the solver
    :param bool do_print: a flag indicating whether or not to print the elapsed time
    
    :param function f_eval: 导数函数 :math:`y' = f(t,y)`
    :param float t_init: 初始时间 :math:`t_{init}`
    :param float t_final: 最终时间 :math:`t_{final}`
    :param float dt_init: 初始步长 :math:`\Delta{t}_{init}`
    
    :param params.Params p: 求解器相关参数对象
    
    :param numpy.ndarray y0: 初始条件 :math:`y(t_{init})` (长度 = m)
    :param float tol: 容差，控制数值解与精确解的最大误差 :math:`\||y - \tilde{y} \|_{\infty} \le tol`
    :param int n_iter_max_newton: 牛顿迭代的最大次数（默认值：N_ITER_MAX_NEWTON）
    :param float be_tol:  backward Euler 求解器的收敛判据（默认值：BE_TOL）
    :param float sdc_tol: SDC 求解器的收敛判据（默认值：SDC_TOL）
    :param int n_steps_max: 求解器的最大步数（默认值：N_STEPS_MAX_ADAPTIVE）
    :param bool do_print: 是否打印计算耗时信息的标志位
    :return: all of the time nodes, the value of the solution at each node, a history of the solution \
    history for each iteration at each time step, history of the deferred correction for each \
    iteration at each time step, the step size for each time step

    :rtype: numpy.ndarray, numpy.ndarray, list, list, numpy.ndarray
    """

    assert t_final > t_init

    assert dt_init <= (t_final - t_init)

    # 开始计时
    start = time.time()

    # 积分步长
    h = dt_init

    # 每个时间步的时间节点数量
    n_nodes = p.n_nodes

    # 节点类型
    node_type = p.node_type

    # 存储所有解、时间节点、解的历史记录和校正的历史记录
    y_all, t_all, Y_all, D_all = list(), list(), list(), list()

    y_all.append(y0)
    t_all.append(np.array(t_init))
    h_all = list()

    # 步数计数器
    j = 0

    # 当前时间步的初始时间
    t0 = t_init

    # 当前时间步的初始条件
    v0 = np.array(y0)

    do_stop = False

    scaler_big_enough = 1.5

    # print("adaptive run")
    while (t0 < t_final) and (not do_stop):

        assert j < n_steps_max, 'Done too many time steps: ' + str(j) + '. Quiting!'

        # print("t0: ",t0)
        # print("h: ", h)
        # print("v0信息： ")
        # print(v0)
        # 使用1步运行JFNK方法
        t1, y1, Y1, D1 = jfnk_uniform(f_eval, t_init=t0, t_final=(t0 + h), n_steps=1, p=p, y0=v0, \
                                      n_iter_max_newton=n_iter_max_newton, be_tol=be_tol, sdc_tol=sdc_tol)

        # print("t1信息： ")
        # print(t1)
        # print("y1信息： ")
        # print(y1)
        
        # ----------------------------------------------------------------------------------------------------
        # 使用2步运行JFNK方法
        # t2, y2, Y2, D2 = jfnk_uniform(f_eval, t_init=t0, t_final=(t0 + h), n_steps=2, p=p, y0=v0, \
        #                               n_iter_max_newton=n_iter_max_newton, be_tol=be_tol, sdc_tol=sdc_tol)
        
        
        # print("t2信息： ")
        # print(t2)
        # print("y2信息： ")
        # print(y2)
        
        # 计算步长缩放因子
        # scaler = step_size_scaler(y1[-1], y2[-1], n_nodes, 2, tol, node_type)
        # print("scaler信息： ")
        # print(scaler)
        # 调整缩放因子
        # scaler = adjust_scaler(scaler, scaler_big_enough)
        # print("scaler信息： ")
        # print(scaler)
        # 仅当步长变小或缩放因子足够大时重新运行积分
        # do_rerun = (scaler < 1) or (scaler > scaler_big_enough)
        # print("do_rerun信息： ")
        # print(do_rerun)
        # 使用新的步长重新运行积分
        # if do_rerun:

            # 使用新的步长
            # h = update_step_size(scaler * h, t0, t_final)
            # print("------------do_rerun--------")
            # print("do_rerun h: ", h)
            # print("do_rerun t0: ", t0)
            # print("do_rerun t0 + h: ", t0 + h)
            # print("do_rerun v0:")
            # print(v0)
            # 使用改进后的时间步重新运行JFNK方法
            # t, y, Y, D = jfnk_uniform(f_eval, t_init=t0, t_final=(t0 + h), n_steps=1, p=p, y0=v0, \
            #                           n_iter_max_newton=n_iter_max_newton, be_tol=be_tol, sdc_tol=sdc_tol)
        # else:
            # 使用原始近似
            # t, y, Y, D = t1, y1, Y1, D1
        # print("t信息： ")
        # print(t)
        # print("y信息： ")
        # print(y)
        
        t, y, Y, D = t1, y1, Y1, D1
        # ----------------------------------------------------------------------------------------------------
        
        # 存储当前时间步的值
        if p.t[-1] == 1:
            y_all.append(y[1:])
            t_all.append(t[1:])
        else:
            y_all.append(y)
            t_all.append(t)

        h_all.append(h)
        Y_all.append(Y[0])
        D_all.append(D[0])

        if j%100==0:
            print("j: ", j)
            print(t0," / ", t_final)
            curTime = time.time()
            print("runing time:")
            print_elapsed_time(start, curTime)
            print("y信息")
            print(y)
            # print("y_all信息")
            # print(y_all)
            print()

        #
        # 更新下一个时间步的初始时间
        #
        t0 = t0 + h
        
        # 更新下一个时间步的初始条件
        v0 = np.array(y[-1])

        # 检查是否到达结束时间
        do_stop = stopping_criteria(t_final, t0)
        # print("do_stop信息： ")
        # print(do_stop)
        # 防止超过最终时间
        # h = update_step_size(h, t0, t_final)
        # print("h信息： ")
        # print(h)
        # 更新计数器
        j = j + 1

    # 存储解和时间
    y_all = np.vstack(y_all)
    t_all = np.hstack(t_all)
    h_all = np.hstack(h_all)

    # 结束计时
    end = time.time()

    if do_print:
        print("jfnk_adaptive run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    return t_all, y_all, Y_all, D_all, h_all

def jfnk_initial(f_eval, t, y0, y_approx, S, S_p, spectral_radius, be_tol=BE_TOL, sdc_tol=SDC_TOL, \
                 n_iter_max=N_ITER_MAX_SDC, do_print=False):

    """
    This function runs initial SDC iterations until order convergence is observed. Given
    :math:`H(y^{[k]}) = \\delta^{[k]}` corresponds to one iteration of the SDC method.

    1. Run 2 initial SDC iterations

    .. math::

		\\begin{cases}
		    \\delta^{[k]} & \\leftarrow H(y^{[k]}) \\\\
		    y^{[k+1]} & \\leftarrow y^{[k]} + \\delta^{[k]},  k=0, 1
		\\end{cases}

    2. Caculate the ratio of the corrections

    .. math::
        r = \\frac{\\| \\delta^{[k-1]} \\|_F}{ \\| \\delta^{[k-2]} \\|_F}

    where :math:`\\| \\cdot \\|_F` is the Frobenius norm.

    If :math:`\\frac{r}{ \\rho(C_s) } > 0.1`, order convergence is observed. Stop the function.

    If :math:`\\frac{r}{ \\rho(C_s) } \\le 0.1`, if not converged, do another SDC iteration and go to step 2.

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes over a time step :math:`\\Delta{t}`
    :param numpy.ndarray y0: the initial condition
    :param numpy.ndarray y_approx: the provisional solution (dimensions n time nodes, size of the problem)
    :param S: the spectral (Gaussian quadrature) integration matrix
    :param S_p: the backward Euler integration matrix
    :param float spectral_radius: the spectral radius of the correction matrix for the extremely stiff case
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param n_iter_max: the maximum number of SDC iterations
    :param bool do_print: a flag indicating whether or not to print the elapsed time

    :return: the approximation, the history of the approximations for each iteration, the history of the \
    deferred corrections for each iteration, a flag indicating whether or not the solution has converged, \
    a flag indicating whether or not the problem is stiff, relative magnitude of consecutive iterations

    :rtype: numpy.ndarray, list, list, bool, bool, numpy.ndarray
    """
    # 开始计时
    start = time.time()

    # 标记是否为刚性系统
    is_stiff = False

    # 计数器
    i = 0

    # 存储连续迭代的校正值比值
    ratios = []

    # print("jfnk_initial i: ", -2)
    # print("y_approx: ")
    # print(y_approx)
    # 运行初始2次迭代
    y_sdc, Y_sdc, D_sdc, is_converged = sdc(f_eval, t, y0, y_approx, S, S_p, n_iter_max_sdc=2, be_tol=be_tol, \
                                            sdc_tol=sdc_tol, do_print=False)
    
    # 打印初始迭代信息
    # print("jfnk_initial i: ", -1)
    # print("y_sdc: ")
    # print(y_sdc)
    
    # 存储之前的近似解
    Y = [x for x in Y_sdc]

    # 存储之前的校正值
    D = [x for x in D_sdc]

    if is_converged:
        ratios = np.array( [] )

    # JFNK初始解的迭代过程
    # print("jfnk_initial")
    while (not is_converged) and (not is_stiff) and (i < n_iter_max):

        # 计算校正值的比值
        ratio = LA.norm(D[-1], ord='fro') / LA.norm(D[-2], ord='fro')

        ratios.append(ratio)

        # 判断是否为刚性系统
        is_stiff = (ratio / spectral_radius) > 0.1
        # 打印当前迭代信息
        # print(i, " ratio: ", ratio)
        # print("spectral radius: ", spectral_radius)
        # print("is_stiff: ", is_stiff)

        if is_stiff:
            if do_print:
                print('Stiff (order convergence detected). Use JFNK.')
        else:
            if do_print:
                print('Non-stiff. Use SDC')

            # 继续进行SDC迭代
            y_sdc, Y_sdc, D_sdc, is_converged = sdc(f_eval, t, y0, y_sdc, S, S_p, n_iter_max_sdc=1, \
                                                    be_tol=be_tol, sdc_tol=sdc_tol, do_print=False)
            # 存储近似解
            Y.append(Y_sdc[0])

            # 存储校正值
            D.append(D_sdc[0])

        # 打印当前迭代信息
        # print("jfnk_initial i: ", i)
        # print("y_sdc: ")
        # print(y_sdc)
        # 更新计数器
        i = i + 1

    if len(ratios) != 0:
        ratios = np.vstack(ratios)

    # 结束计时
    end = time.time()

    if do_print:
        print("jfnk_initial run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")
    
    return y_sdc, Y, D, is_converged, is_stiff, ratios

def jfnk_iterations(f_eval, t, y0, y_init, S, S_p, n_iter_max_newton=N_ITER_MAX_NEWTON, be_tol=BE_TOL, sdc_tol=SDC_TOL, do_print=False):

    """
    This function solves the Newton's method iterations for solving

    .. math::
        H(y) = 0

    where :math:`H(y^{[k]}) = \\delta^{[k]}` corresponds to one iteration of the SDC method.

    1. Run :math:`n+1` SDC iterations.

    .. math::

        \\begin{cases}
		    \\delta^{[k]} & \\leftarrow H(y^{[k]}) \\\\
		    y^{[k+1]} & \\leftarrow y^{[k]} + \\delta^{[k]},  k=0, \\ldots, n
		\\end{cases}

    2. Solve the Newton iteration system without using the Jacobian explicitly

    .. math::

            J_{H}(y^{[n]})\\Delta{x} &= -H(y^{[n]})   \\\\
            \\implies J_{H}(y^{[n]})\\Delta{x} &= -\\delta^{[n]}   \\\\

    Set :math:`\\Delta{x} = \\sum^{n-1}_{j=0} c_j \\delta^{[j]}` and solve

    .. math::
        J_{H}(y^{[n]})\\sum^{n-1}_{j=0} c_j \\delta^{[j]} &= -\\delta^{[n]}  \\\\
        \\implies \\sum^{n-1}_{j=0} c_j (\\delta^{[k+1]} - \\delta^{[k]}) &= -\\delta^{[n]} \\\\

    Solve the system for the Jacobian-Free system

    .. math::
        \\begin{cases}
            Ac &= -\\delta^{[n]} \\\\
            y &\\leftarrow y^{[n]} + \\sum^{n-1}_{j=0} c_j \\delta^{[j]}
        \\end{cases}

    3. If not converged, repeat by going to step 1.

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes in the time step (length n_nodes)
    :param numpy.ndarray y0: the initial solution [size of the problem x 1]
    :param numpy.ndarray y_init: the approximate solution [n_nodes x size of the problem]
    :param numpy.ndarray S: the spectral integration matrix 
    :param numpy.ndarray S_p: the preconditioner integration matrix 
    :param int n_iter_max_newton: the maximum number of Newton iterations
    
    :return: the solution, the history of the solutions for each iteration, the deferred \
    correction for each iteration, a flag indicating whether the procedure converged

    :rtype: numpy.ndarray, list, list, bool
    """
    # 开始计时
    start = time.time()

    # 当前解
    y = np.array( y_init )
    # 打印初始解
    # print("y_init: ")
    # print(y_init)

    # SDC迭代次数
    k_sdc = len(t) + 1

    # 系统的规模
    m = len(y0)

    # 校正值和近似解的历史记录
    D_list, Y_list = list(), list()

    k = 0
    is_converged = False
    # 记录开始时间
    start = time.time()

    # 牛顿迭代过程
    while (not is_converged) and (k < n_iter_max_newton):

        # 进行SDC迭代
        y, Y, D, is_converged = sdc(f_eval, t, y0, y, S, S_p, n_iter_max_sdc=k_sdc,\
                                    be_tol=be_tol, sdc_tol=sdc_tol)

        # 打印当前迭代信息
        
        # print("k: ", k, " / ", n_iter_max_newton)
        # print (k_sdc)
        # 计算并打印运行时间
        # end = time.time()
        # print("jfnk run time: ", end - start, " s")
        # print("is_converged: ", is_converged)
        # print("y: ")
        # print(y)

        if is_converged:
            # 存储SDC扫描的解历史记录
            Y_list.append(Y)
            # 存储校正值的历史记录
            D_list.append(D)
        else:
            # 记录开始时间
            is_converged_start = time.time()
            # print("\n!is_not_converged")
            # 存储SDC扫描的解历史记录（除最后一个）
            Y_list.append(Y[:-1])
            # 存储校正值的历史记录（除最后一个）
            D_list.append(D[:-1])

            # 该解用于克里洛夫迭代
            y = Y[-1]

            #
            # 设置
            #

            # 转置数据
            y = y.T
            D = [x.T for x in D]

            # 构建待求解系统
            A = [ (D[j + 1] - D[j]) for j in range(k_sdc - 1)]

            # 延迟校正向量的“基”
            M = [ D[j] for j in range(k_sdc - 1) ]

            #
            # 牛顿方法
            #

            # 对每个未知数
            for i in range(m):

                # 求解牛顿迭代
                B = [a[i] for a in A]
                B = np.vstack(B).T
                rhs = -D[-1][i]
                rhs = rhs.reshape( (len(rhs), 1) )

                # 求解系统Bc = rhs
                c, res, rank, s = np.linalg.lstsq(B, rhs)

                # 创建更新量dy = Vc
                V = [x[i] for x in M]
                V = np.vstack(V).T
                dy = V.dot(c)

                # 更新
                y[i, :] += dy[:].flatten()

            # 转置回来
            y = y.T
            # 记录结束时间
            # is_converged_end = time.time()
            # 计算并打印运行时间
            # print("is_converged run time: ", is_converged_end - is_converged_start, " s")
            # print("y: ")
            # print(y)
            # print()

        # 更新迭代次数
        k = k + 1

    # 将完整的历史记录存储在一个列表中
    Y_list = [item for subitem in Y_list for item in subitem ]

    D_list = [item for subitem in D_list for item in subitem]
    # 结束计时
    end = time.time()

    if do_print:
        print("jfnk_iterations run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    return y, Y_list, D_list, is_converged

def jfnk_step(f_eval, p, dt, t, y0,  n_iter_max_newton=N_ITER_MAX_NEWTON, be_tol=BE_TOL, \
              sdc_tol=SDC_TOL, do_print=False):
    """
    This function runs everything needed for the JFNK to run for 1 time step (i.e., \
    the backward euler precondition and the JFNK iterations).

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param params.Params p: the parameters related to the solver
    :param float dt: the time step :math:`\\Delta{t}`
    :param numpy.ndarray t: the number of temporal nodes (length number of time nodes)
    :param numpy.ndarray y0: the initial condition (length, size of the problem)
    :param n_iter_max_newton: the maximum number of Newton iteration
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param bool do_print: a flag indicating whether or not to print the elapsed time

    :return: the approximation, the history of the approximations for each iteration, the history of the \
    deferred corrections for each iteration, a flag indicating whether or not the solution has converged, \
    a flag indicating whether or not the problem is stiff, relative magnitude of consecutive iterations

    :rtype: numpy.ndarray, list, list, bool, bool, numpy.ndarray
    """

    # 开始计时
    start = time.time()

    # 谱积分矩阵和预条件积分矩阵
    S, S_p = dt * p.S, dt * p.S_p

    # 谱半径
    spectral_radius = p.spectral_radius

    # 向后欧拉近似解
    y_be = backward_euler(f_eval, t, y0, S_p, be_tol=be_tol)
    # 打印初始迭代信息
    # print("y0: ")
    # print(y0)
    # print("y_be: ")
    # print(y_be)
    # JFNK求解器
    y, Y, D, is_converged, is_stiff, ratios = jfnk(f_eval, t, y0, y_be, S, S_p, spectral_radius, \
                                                   n_iter_max_newton, be_tol=be_tol, sdc_tol=sdc_tol)

    # print("t.size(): ", t.size)
    # print(t)
    # print("y.size(): ", y.size)
    # print(y)
    # 结束计时
    end = time.time()

    if do_print:
        print("jfnk_step run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    return y, Y, D, is_converged, is_stiff, ratios

def jfnk_uniform(f_eval, t_init, t_final, n_steps, p, y0, n_iter_max_newton=N_ITER_MAX_NEWTON, \
                 be_tol=BE_TOL, sdc_tol=SDC_TOL, do_print=False):

    """
    Run the JFNK with uniform step sizes from :math:`t \\in [t_{init}, t_{final}]` to calculate an
    approximation to the solution

    .. math::
        y(t_{final}) = y(t_{init}) +  \int^{t_{final}}_{t_{init}} \\, f(\\tau, y(\\tau) ) \\, \\mathrm{d}\\tau

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param float t_init: the initial time :math:`t_{init}`
    :param float t_final: the final time :math:`t_{final}`
    :param int nsteps: the number of steps
    :param params.Params p: the parameters related to the solver
    :param numpy.ndarray y0: the initial condition :math:`y(t_{init})` (length = m)
    :param int n_iter_max_newton: the maximum number of Newton iterations
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param bool do_print: a flag indicating whether or not to print the elapsed time
    
        :param function f_eval: 导数函数 :math:`y' = f(t,y)`
    :param float t_init: 初始时间 :math:`t_{init}`
    :param float t_final: 最终时间 :math:`t_{final}`
    :param int n_steps: 总步数（将时间区间均匀分成n_steps段）
    :param params.Params p: 求解器相关参数对象
    :param numpy.ndarray y0: 初始条件 :math:`y(t_{init})` (长度 = m)
    :param int n_iter_max_newton: 牛顿迭代的最大次数（默认值：N_ITER_MAX_NEWTON）
    :param float be_tol: backward Euler 求解器的收敛判据（默认值：BE_TOL）
    :param float sdc_tol: SDC 求解器的收敛判据（默认值：SDC_TOL）
    :param bool do_print: 是否打印计算耗时信息的标志位

    :return: all of the time nodes, the value of the solution at each node, a history of the solution \
    history for each iteration at each time step, history of the deferred correction for each \
    iteration at each time step

    :rtype: numpy.ndarray, numpy.ndarray, list, list
    """

    assert t_final > t_init

    # 开始计时
    start = time.time()

    # 积分步长
    dt = (t_final - t_init) / n_steps

    # 初始条件
    v0 = np.array(y0)

    # 归一化的节点，范围在[0, 1]
    if p.t[-1] != 1:
        t_step = np.hstack( [p.t, 1.0] )
    else:
        t_step = np.array(p.t)

    # 将节点缩放至正确的时间步???
    nodes = t_init + dt * t_step
    print("nodes: ", nodes)
    # 存储所有近似解、时间节点、解的历史记录和校正的历史记录
    y_all, t_all, Y_all, D_all = list(), list(), list(), list()

    # 存储初始近似解和时间
    y_all.append(y0)
    t_all.append( np.array(t_init) )

    # 遍历时间步
    for k in range(n_steps):

        # if (n_steps>1):
        #     print("n_steps k: ", k)
        #     print("v0: ")
        #     print(v0)
        #     print("dt: ", dt)
        #     print("nodes: ")
        #     print(nodes)
        # 使用1步运行JFNK求解器
        y, Y, D, is_converged, is_stiff, ratios = jfnk_step(f_eval, p, dt, nodes, v0, n_iter_max_newton, \
                                                            be_tol=be_tol, sdc_tol=sdc_tol)
        # if (n_steps>1):
        #     print("k: ", k)
        #     print("y: ")
        #     print(y)
        
        assert is_converged, 'the JFNK did not converge on time step: ' + str(k) + '. Quitting...'

        # 存储值
        y_all.append(y[1:])
        t_all.append(nodes[1:])

        Y_all.append(Y)
        D_all.append(D)

        # 更新下一个时间步的初始值
        v0 = np.array(y[-1])

        # 更新时间节点
        nodes = nodes + dt

    # 存储值
    y_all = np.vstack(y_all)
    t_all = np.hstack(t_all)

    # 结束计时
    end = time.time()

    # 打印耗时
    if do_print:
        print("jfnk_uniform run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    return t_all, y_all, Y_all, D_all

def print_elapsed_time(start, end):

    """
    This function prints the elapsed time [s] between to time points.

    :param float start: the start time [s]
    :param float end: the end time [s]
    :return:
    """

    print('elapsed time: %.2f[s]' % (end - start))

    return

def relative_norm(top, bot):

    """
    This function calculates the relative ratios between norms of vectors. Given two sets of vectors

    .. math::
        top_{n \\times m}, bot_{n \\times m}

    where :math:`n` is the number of temporal nodes and :math:`m` is the size of the system.

    For each component :math:`j`, calculate the relative norms over the time nodes

    .. math::
        ratio_j = \\frac{ \\| top_j \\|_2 }{\\| bot_j \\|_2}

    Make sure, we avoid division by zero

    .. math::
        x_j =
        \\begin{cases}
            ratio_j & \\text{if } \\|bot_j \\|_2 \\neq 0 \\\\
            \\| top_j \\|_2 & \\text{if } \\| bot_j \\|_2 = 0
        \\end{cases}

    Return the maximum value

    .. math::
        \\| x \\|_{\\infty}

    :param top: the top vector for a given iteration dimensions (n nodes, m size of problem)
    :param bot: the bottom vector for a given iteration  dimensions (n nodes, m size of problem)

    :return: the maximum value of the relative norm between two vectors
    :rtype: float
    """

    # 计算分子向量的范数
    top_norm = LA.norm(top, axis=0)

    # 计算分母向量的范数
    bot_norm = LA.norm(bot, axis=0)

    # 计算相对校正值的大小
    ratio = top_norm / bot_norm

    # 分母不为零的分量的索引
    idx = np.isfinite(ratio)

    if idx.all():
        # 没有除以零的情况
        value = ratio.max()
    else:
        # 至少有一个分量除以零
        # 对于除以零的分量使用绝对校正值
        # 对于不除以零的分量使用相对校正值
        value = max(top_norm[~idx].max(), ratio[idx].max())

    return value

def run_spectral(f_eval, t, y0, S, tol, y_approx, do_print=False, verbose=False):

    start = time.time()

    y_spect = spectral(f_eval, t, y0, S, tol, y_approx=y_approx, verbose=verbose)

    end = time.time()

    if do_print:
        print_elapsed_time(start, end)

    return y_spect

def scaler_lobatto(aerr, n_nodes, k, tol):

    """
    This function calculates :math:`x` the amount the step size should increase \
    or decrease for the proper adaptive time step size for Gauss-Lobatto nodes.

    :param numpy.ndarray aerr: the maximum absolute error between the approximate solution :math:`y_h` \
    and the higher accuracy solution :math:`y_{h/k}`
    :param int n_nodes: the number of temporal nodes within the time step
    :param int k: the amount of mini time steps leading up to :math:`y(\\Delta{t})`
    :param float tol: the approximated absolute error at each step for the adaptive solution

    .. math::
       x =  \\left( tol \\frac{ 1 - (\\frac{1}{k})^{p} }{ \\| y_{h} - y_{h/k}\\|_{\\infty}} \\right)^{1/p}

    where :math:`p = 2 * n - 1` and :math:`n` is the number of Gauss-Lobatto nodes

    :return: the amount the current step size should increase or decrease for the proper adaptive time step
    :rtype: float
    """

    p = 2 * n_nodes - 1
    # print("p: ", p)
    x = tol * (1 - (1 / k)** p) / aerr
    # print("x: ", x)
    scaler = x ** (1 / p)
    # print("scaler_lobatto: ", scaler)
    return scaler

#当前迭代的近似解与理想解之间的偏差
def sdc(f_eval, t, y0, y_old, S, S_p, n_iter_max_sdc=N_ITER_MAX_SDC, be_tol=BE_TOL, 
        sdc_tol=SDC_TOL, do_print=False):

    """
    This function runs the SDC method until convergence.

    .. math::
        \\begin{cases}
		    \\delta^{[k]} & \\leftarrow H(y^{[k]}) \\\\
		    y^{[k+1]} & \\leftarrow y^{[k]} + \\delta^{[k]}
		\\end{cases}

    where :math:`H(y^{[k]}) = \\delta^{[k]}` corresponds to one iteration of the SDC method.

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the temporal nodes (length number of nodes) over the time step of size :math:`\\Delta{t}`
    :param numpy.ndarray y0: the initial condition :math:`y(t_{init})` (length = size of the problem)
    :param numpy.ndarray y_old: the provisional solution (dimensions, size of the problem)
    :param numpy.ndarray S: the spectral integration (Gaussian quadrature) matrix, :math:`S`
    :param numpy.ndarray S_p: the backward Euler integration matrix, :math:`\\tilde{S}`
    :param int n_iter_max_sdc: the maximum number of SDC iterations
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param bool do_print: a flag indicating whether or not to print the elapsed time

    :return: the SDC solution, a history of the approximate solution for each iteration, the history of \
    the deferred correction for each itearation, a flag indicating whether or not the solution has met the \
    convergence criteria

    :rtype: numpy.ndarray, list, list, bool
    """

    # 开始计时
    start = time.time()

    # 每个时间步的导数向量
    F = np.zeros( y_old.shape )

    # 时间节点的数量，系统的规模
    n_nodes, m = F.shape

    # 解的完整历史记录
    Y = list()

    # 校正值的历史记录
    D = list()

    # 初始猜测
    y_sdc = np.array(y_old)

    # 计数器
    k = 0

    # 标记SDC是否收敛
    is_converged = False

    while (not is_converged) and  (k < n_iter_max_sdc):

        # 存储解
        Y.append(y_sdc)

        F[:] = 0
        for i in range(n_nodes):
            F[i] = f_eval(t[i], y_sdc[i])
            # if (n_iter_max_sdc > 2):
            #     print("\nbuild F: i: ", i)
            #     print("t[i]: ", t[i])
            #     print("y_sdc[i]: ", y_sdc[i])
            #     print("F[i]: ", F[i])

        # 进行SDC扫描
        
        # if (n_iter_max_sdc>2):
        #     print("\n\nstart sweep k: ", k)
            # y_sdc, d = sdc_sweep(f_eval, t, y0, y_sdc, F, S, S_p, be_tol=be_tol,do_print=True)
            # print("\nsdc k: ", k)
            # print("y_sdc: ")
            # print(y_sdc)
            # print("")
        # else:
        y_sdc, d = sdc_sweep(f_eval, t, y0, y_sdc, F, S, S_p, be_tol=be_tol)
        # 检查SDC是否收敛
        is_converged = convergence_criteria(d, Y[-1], tol=sdc_tol)

        # 存储校正值
        D.append(d)

        # 更新计数器
        k = k + 1

    Y = np.vstack(Y)
    D = np.vstack(D)

    # 调整数组形状
    Y.resize( (k, n_nodes, m) )
    D.resize( (k, n_nodes, m) )

    # 停止计时
    end = time.time()

    # 打印耗时
    if do_print:
        print("sdc run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    return y_sdc, Y, D, is_converged

def sdc_node(f_eval, t, h, rhs, y_old, be_tol=BE_TOL, do_print=False):

    """
    This function runs SDC on the :math:`i^{th}` time node within the time step.

     .. math::
            y^{[k+1]}_i - hF(y^{[k+1]}_i) = rhs_i

    where

    .. math::
        rhs_i = y^{[k]}_{i-1} + (\\Delta{t}S_i - he_i)\\cdot F(y^{[k]})

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param float t: the time
    :param float h: the time step size
    :param numpy.ndarray rhs: the right hand side of the backward euler system [m x 1]
    :param numpy.ndarray y_old: the approximate solution before the update [number of nodes, size of the problem)
    :param float be_tol: the convergence criteria for the backward Euler solver

    :return: both the improved solutionand the correction at the time :math:`t`
    :rtype: numpy.ndarray, numpy.ndarray
    """
    # 开始计时
    start = time.time()
    # print("sdc node, rhs")
    # print(rhs)
    # print("sdc node, y_old")
    # print(y_old)
    # 在时间节点t处的改进解
    y_new = backward_euler_node(f_eval, t, h, rhs=rhs, x0=y_old, be_tol=be_tol)

    # 延迟校正
    d = y_new - y_old
    # 结束计时
    end = time.time()

    if do_print:
        print("sdc_node run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")
    return y_new, d

def sdc_sweep(f_eval, t, y0, y_old, F, S, S_p, be_tol=BE_TOL, do_print=False):

    """
    This runs SDC on the entire time step interval

    .. math::
        \\begin{cases}
            y^{[k+1]} - \\Delta{t}\\tilde{S}F(y^{[k+1]}) & = y_0 + \\Delta{t}(S- \\tilde{S})F(y^{[k]}) \\\\
            \\delta^{[k]} &=  H(y^{[k]}) = y^{[k+1]} - y^{[k]}
        \\end{cases}

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes over a time step :math:`\\Delta{t}`
    :param numpy.ndarray y0: the initial solutions length, size of the problem
    :param numpy.ndarray y_old: the approximate solution :math:`y^{[k]}` [n_nodes x m]
    :param numpy.ndarray F: the derivative at the approximate solution :math:`F(y^{[k]})`
    :param numpy.ndarray S: the spectral integration (Gaussian quadrature) matrix :math:`S`
    :param numpy.ndarray S_p: the backward Euler integration matrix :math:`\\tilde{S}`
    :param float be_tol: the tolerance for the backward Euler solver

    :return: the approximate solution, the correction
    :rtype: numpy.ndarray (dimensions number of temporal nodes x size of the problem), \
    numpy.ndarray (dimensions number of temporal nodes by size of the problem)
    """
    # 开始计时
    start = time.time()

    # 时间节点的数量
    n_nodes = len(t)

    # 系统的规模
    m = len(y0)

    # 时间步内的解
    y = np.zeros( (n_nodes, m) )

    # 时间步内的校正值
    d = np.zeros(y.shape)

    #
    # 对整个时间步运行SDC
    #
    # 打印SDC扫描运行信息
    # print("sdc sweep run")
    # print("y0: ")
    # print(y0)
    # print("S: ")
    # print(S)
    for i in range(n_nodes):
                
        # print("i: ", i)
        # 预条件积分中的步长
        h = S_p[i, i]
        w = np.zeros(n_nodes)
        w[i] = h

        # 如果时间节点包含时间步的左端点则跳过
        do_skip =  (S_p[i] == 0).all()

        if i == 0 and do_skip:
            y[0,:] = y0[:]
        else:
            # print("y0: ")
            # print(y0)
            # print("S: ")
            # print(S)
            # print("do_skip", do_skip)
            if (i == 0) and (not do_skip):
                rhs = y0 + (S[i,:] -  w).dot(F)
            else:
                rhs =  y[i-1] + (S[i,:] - S[i-1,:] - w).dot(F)
            # print("rhs: ")
            # print(rhs)
            # print("w: ")
            # print(w)
            # print("F: ")
            # print(F)
            # print("val: ")
            # print((S[i,:] - S[i-1,:] - w).dot(F))
            # 在当前时间求解SDC
            y_temp, d_temp = sdc_node(f_eval, t[i], h, rhs, y_old[i], be_tol=be_tol)

            # 存储值
            y[i,:]  = y_temp
            d[i,:]  = d_temp
            
        # 打印当前迭代信息

        # print("y: ")
        # print(y[i,:])
    # 结束计时
    end = time.time()

    if do_print:
        print("sdc_sweep run time: ")
        # print_elapsed_time(start, end)
        print(end - start, " s")

    # print("sdc_sweep y: ")
    # print(y)
    return y, d

def spectral(f_eval, t, y0, S, f_tol=NK_FTOL, y_approx=None, do_print=False, verbose=False):

    """
    This function solves directly the spectral solution (Gauss collocation formulation). That is,
    this function solves

    .. math::
        y - \\Delta{t}SF(y) = y_0

    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes over a time step :math:`\\Delta{t}`
    :param numpy.ndarray y0: the initial solutions length, size of the problem
    :param numpy.ndarray S: the spectral integration matrix :math:`S`
    :param float f_tol: the relative tolerance of the Newton-Krylov solver
    :param numpy.ndarray y_approx: the initial guess for the Newton-Krylov solver
    :param bool verbose: a flag for the built-in Newton-Krylov solver

    :return: the spectral solution
    :rtype: numpy.ndarray
    """

    # 开始计时
    start = time.time()

    # 时间节点的数量
    n_nodes = S.shape[1]

    # 复制初始条件，形状为 [n_nodes x m]
    y0_vec  = np.vstack([y0 for _ in range(n_nodes)])

    # 定义需要求解的根函数
    H       = lambda y: spectral_root_finder(f_eval, t, y, y0_vec, S)

    # 获取解
    if y_approx is None:
        y_approx = y0_vec

    # 使用牛顿-克里洛夫方法求解
    y   = newton_krylov(H, y_approx, f_rtol=f_tol, verbose=verbose)

    # 结束计时
    end = time.time()

    if do_print:
        print_elapsed_time(start, end)

    return y

def spectral_root_finder(f_eval, t, y, y0_vec, S):

    """
    This function calculates the residual in the spectral (Gauss) collocation formulation. That is,

    .. math::
        y - \\Delta{t}SF(t,y) - y_0.

    It is used in the root finding algorithm to solve

    .. math::
        A(y) = y - \\Delta{t}SF(t,y) - y_0 = 0
    
    :param function f_eval: the derivative function :math:`y' = f(t,y)`.
    :param numpy.ndarray t: the time nodes
    :param numpy.ndarray y: an approximate solution [number of temporal nodes x size of problem]
    :param y0_vec: the initial condition vector [number of temporal nodes x size of problem]
    :param numpy.ndarray S: the spectral integration (Gauss quadrature) matrix :math:`S`

    :return: the residual in the spectral collocation formulation
    :rtype:  numpy.ndarray  [ number of temporal nodes x size of problem ]
    """

    # 时间节点的数量，系统的规模
    n_nodes, m = y0_vec.shape

    # 导数函数
    F       = np.zeros(y0_vec.shape)
    temp    = np.zeros(y0_vec.shape)

    # 计算所有节点的导数
    for i in range(n_nodes):
        F[i] = f_eval(t[i], y[i])

    # 谱积分
    for i in range(m):
        temp[:, i] = S.dot(F[:, i])

    # 残差
    res = y - temp - y0_vec

    return res

def step_size_scaler(y, ysteps, n_nodes, k, tol, node_type):

    """
    This function

    :param numpy.ndarray y: an approximation of the ODE system using 1 step size of :math:`h=\\Delta{t}`
    :param numpy.ndarray ysteps: an approximation of the ODE system using :math:`k` steps of size \
    :math:`h=\\frac{\\Delta{t}}{k}`
    :param int n_nodes: the number of temporal nodes within the time step
    :param int k: the amount of mini time steps leading up to :math:`y(\\Delta{t})`
    :param float tol: the absolute error tolerance wanted within the time step
    :param int node_type: the type of temporal nodes

    Given a step size :math:`\\Delta{t}`, calculate :math:`x` the amount the step size should increase \
    or decrease for the proper adaptive time step size :math:`\\Delta{t}_{new}` where

    .. math::
        \\Delta{t}_{new} = x \\Delta{t}

    .. note::
        This function currently runs using only Gauss-Lobatto nodes

    :return:
    """

    # 仅支持高斯-洛巴托节点
    assert node_type == points.GAUSS_LOBATTO

    # aerr是在时间dt处的最大误差
    aerr = LA.norm(y - ysteps, ord=np.inf)
    # print("aerr信息： ")
    # print(aerr)
    # k是达到步长dt所需的子步数
    if node_type == points.GAUSS_LOBATTO:
        scaler = scaler_lobatto(aerr, n_nodes, k, tol)
        # print("scaler信息： ")
        # print(scaler)

    return scaler

def stopping_criteria(t_final, t0):

    """
    This function sends a flag whether or not we have reached the end of the simulation while \
    taking account errors from inexact arithmetic.

    :param float t_final: final time in the ODE simulation
    :param float t0: the start time for the current time step

    :return: a flag indicating whether or not we have reached the end of the simulation
    :rtype: bool
    """

    # 定义最终时间“小”的标准
    T_SMALL = 1e-15

    # 用于补偿不精确算术误差
    EPS = 1e-20

    # 考虑不精确算术误差后的时间
    t = t0 + EPS

    # 最终时间被认为“大”，使用绝对差值
    # 或者最终时间为0，避免除以0
    if (t_final >= T_SMALL) or (t_final == 0):
        do_stop = t >= t_final
    else:
        # 最终时间被认为“小”，使用相对差值
        # 如果t > t_final 或者 abs(t)/t_final < 相应的量级，则停止
        do_stop = (t_final - t) / t_final <= 1e-5

    return do_stop

def update_step_size(dt, t0, t_final):

    """
    This function makes sure that the time step :math:`\\Delta{t}` does not cause the simulation to go \
    past the final time :math:`t_{final}`.

    :param float dt: the current step size :math:`\\Delta{t}`
    :param float t0: the start time for the current time step
    :param float t_final: the final time :math:`t_{final}` of the ODE simulation

    :return: the step size for the next time step
    :rtype: float
    """

    # 假设步长不变
    dt_new = dt

    # 确保步长不会使模拟超过最终时间
    if (t0 + dt) > t_final:
        dt_new = t_final - t0

    return dt_new
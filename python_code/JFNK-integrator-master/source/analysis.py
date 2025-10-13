# This file was written by Dr. Namdi Brandon
# ORCID: 0000-0001-7050-1538
# June 29, 2018

"""
This module contains class :class:`analysis.Result` that help saves the information \
related to the ODE solution. Also, this module contains functions that aid in \
analyzing the results of the various solutions.
"""

# ===============================================
# import
# ===============================================
# 导入NumPy库，用于数值计算
import numpy as np
# 从NumPy库中导入线性代数模块
import numpy.linalg as LA

# 导入copy模块，用于对象的复制操作
import copy

# ===============================================
# class Result
# ===============================================
class Result(object):

    """
    This class saves the information related to the ODE solution and the parameters used in order \
    to generate the solution.

    :param params.Params p: the parameters related to the solver
    :param float t0: the initial time :math:`t_{init}`
    :param float t_final: the final time :math:`t_{final}`
    :param numpy.ndarray t: the time nodes used in the simulation
    :param numpy.ndarray y: the approximation at the time nodes
    :param bool is_adaptive: a flag indicating whether the simulation used adaptive (if True) \
    step sizes or uniform step sizes (if False)
    :param Y: a history of the solution  history for each iteration at each time step
    :param D: a history of the deferred correction for each iteration at each time step
    :param float h: the adaptive time steps used
    :param float dt_init: the initial step size :math:`\\Delta{t}_{init}` used in the adaptive scheme
    :param float be_tol: the convergence criteria for the backward Euler solver
    :param float sdc_tol: the convergence criteria for the SDC solver
    :param float tol: the approximated absolute error at each step for the adaptive solution
    """

    def __init__(self, p, t0, t_final, t, y, is_adaptive, Y=None, D=None, h=None, dt_init=None, be_tol=None,
                 sdc_tol=None, tol=None):


        # the CMAQ case

        # 复制求解器相关参数
        self.p = copy.copy(p)

        # 求解器输入的初始时间
        self.t0 = t0

        # 求解器输入的最终时间
        self.t_final = t_final

        # 时间步
        self.t = t

        # 解，维度为 (节点总数, 问题规模)
        self.y = np.array(y)

        # 指示结果是否来自自适应算法
        self.is_adaptive = is_adaptive

        # 解的历史记录（即，每个时间步，每次迭代的解）
        if Y is not None:
            Y = copy.copy(Y)
        self.Y = Y

        # 校正的历史记录（即，每个时间步，每次迭代的解）
        if D is not None:
            D = copy.copy(D)
        self.D = D

        # 使用的时间步
        if type(h) is np.ndarray:
            h = np.array(h)
        self.h = h

        # 自适应算法的初始步长
        self.dt_init = dt_init

        # 向后欧拉求解器的容差
        self.be_tol = be_tol

        # SDC迭代的容差
        self.sdc_tol = sdc_tol

        # 自适应时间步的容差
        self.tol = tol

        # 每个时间步的SDC迭代次数和总的SDC迭代次数
        n_sdc_per_step, n_sdc = None, None

        if self.Y is not None:
            # 计算每个时间步的SDC迭代次数
            n_sdc_per_step = np.array([len(x) for x in Y])
            # 计算总的SDC迭代次数
            n_sdc = n_sdc_per_step.sum()

        self.n_sdc_per_step = n_sdc_per_step
        self.n_sdc = n_sdc

        return

# ===============================================
# functions
# ===============================================
def analyze_corrections(D, Y):

    """
    This function calculates the following

    #. for each iteration :math:`k`, the magnitude of the correction for each component \
    :math:`i`: :math:`\\| \\delta^{[k]}_i \\|`
    #. for each iteration :math:`k`, the magnitude of the correction for each component \
    :math:`i` on a :math:`log_{10}` scale : :math:`log_{10} \\left( \\| \\delta^{[k]}_i \\| \\right)`
    #. for each iteration :math:`k`, the magnitude of the relative correction for each \
    component :math:`i`: :math:`\\frac{\\| \\delta^{[k]}_i \\| }{ \\| y^{[k]}_i \\| }`
    #. for each iteration :math:`k`, the magnitude of the relative correction for each \
    component :math:`i` on a :math:`log_{10}` scale: \
    :math:`log_{10} \\left( \\frac{\\| \\delta^{[k]}_i \\| }{ \\| y^{[k]}_i \\| }\\right)`

    :param list D: corrections, dimensions (n iterations, p nodes, m problem size)
    :param list Y: approximations, dimensions (n iterations, p nodes, m problem size)

    :return: for each iteration the magnitude of the correction for each component, \
    for each iteration the magnitude of the correction for each component in log \
    base 10, \
    for each iteration the magnitude of the relative correction for each component, \
    for each iteration the magnitude of the relative correction for each component \
    in log base 10
    :rtype: numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
    """
    # 计算每次迭代中每个分量近似解的范数，结果按行堆叠
    # 用来计算相对范数的分母
    y_norm = np.vstack([np.linalg.norm(x, axis=0) for x in Y])

    # 计算每次迭代中每个分量校正值的范数，结果按行堆叠
    # 代表两次迭代校正的大小（相差的大小），逐渐变小说明在收敛
    d_norm = np.vstack([np.linalg.norm(x, axis=0) for x in D])
    # 计算每次迭代中每个分量的相对校正值
    d_norm_rel  = np.divide(d_norm, y_norm)
    # 计算校正值范数的以10为底的对数
    log_d_norm      = np.log10(d_norm)
    # 计算相对校正值范数的以10为底的对数
    log_d_norm_rel  = np.log10(d_norm_rel)

    return d_norm, log_d_norm, d_norm_rel, log_d_norm_rel

# 全局误差分析，提供整体误差的概览
def error_analysis(y_approx, y_solution, threshold=1e-20):

    """
    This function calculates the absolute error or the relative error.

    :param numpy.ndarray y_approx: the approximate solution
    :param numpy.ndarray y_solution: the more accurate solution
    :param float threshold: the threshold to set the components to 0

    :return: the absolute error, the relative error
    :rtype: numpy.ndarray, numpy.ndarray
    """

    # 计算绝对误差：数值解和准确解
    aerr = np.abs(y_approx - y_solution)

    # 分母
    bot = np.abs(y_solution)

    # 计算相对误差
    rerr = aerr / bot

    # 分母为零或接近零的索引
    idx = bot <= threshold
    # 将分母为零或接近零的相对误差设为0
    rerr[idx] = 0

    return aerr, rerr

def get_correction_norms(Y, D):

    """
    This function calculates various measures of the corrections from SDC sweep

    #. for each iteration :math:`k`, the maximum magnitude of the relative correction for each \
    component :math:`i`: :math:`\\max \\frac{\\| \\delta^{[k]}_i \\| }{ \\| y^{[k]}_i \\| }`
    #. for each iteration :math:`k`, the maximum magnitude of the correction for each component \
    :math:`i`: :math:`\\max \\| \\delta^{[k]}_i \\|`
    #. for each iteration :math:`k`, the mean magnitude of the relative correction for each \
    component :math:`i`: :math:`E\\left[\\frac{\\| \\delta^{[k]}_i \\| }{ \\| y^{[k]}_i \\| }\\right]`
    #. for each iteration :math:`k`, the mean magnitude of the correction for each component \
    :math:`i`: :math:`E[\\| \\delta^{[k]}_i \\|]`

    :param list Y: the approximations for the solution at each iteration arrays of dimensions (n nodes, size of problem)
    :param list D: the corrections for each iteration arrays of dimensions (n nodes, size of problem)

    :return: the maximum relative norm, the maximum absolute norm, the mean relative norm, \
    the mean absolute norm
    :rtype: numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
    """

    # 定义计算相对校正值的函数
    f = lambda y, d: np.linalg.norm(d, axis=0) / np.linalg.norm(y, axis=0)
    # 定义计算校正值范数的函数
    g = lambda y: np.linalg.norm(y, axis=0)

    # 过滤掉非有限值（避免除以零导致的NaN或无穷大）
    f_rel = lambda x: x[np.isfinite(x)]

    # 计算每次迭代的最大相对校正值
    rel_max = np.array([f_rel(f(y, d)).max() for y, d in zip(Y, D)])
    # 计算每次迭代的最大校正值
    ab_max = np.array([g(d).max() for d in D])

    # 计算  
    rel_mean = np.array([f_rel(f(y, d)).mean() for y, d in zip(Y, D)])
    # 计算每次迭代的平均校正值
    ab_mean = np.array([g(d).mean() for d in D])

    return rel_max, ab_max, rel_mean, ab_mean

# 部误差分析，提供每个时间节点的误差细节
def get_error_time_nodes(t, y, y_spline, do_relative=True):

    """
    This function calculates the absolute or relative error comparing the approximate \
    solution and the "exact" solution at all of the temporal nodes.

    :param numpy.ndarray t: the temporal nodes on the approximate solution
    :param numpy.ndarray y: the approximate solution
    :param function y_spline: a function that may interpolate the "exact" (more accurate) solution
    :param bool do_relative: a flag indicating whether (if True) to calculate the relative error \
    or not to (if False) to calculate the absolute error

    :return: the absolute or relative errors
    :rtype: numpy.ndarray
    """

    if do_relative:
        # 计算每个时间节点的相对误差
        err = [ relative_norm( y_spline(t[i]) - y[i], y_spline(t[i]) ) for i in range(len(t)) ]
    else:
        # 计算每个时间节点的绝对误差
        err = [ LA.norm( y_spline(t[i]) - y[i], ord=np.inf ) for i in range(len(t)) ]

    err = np.array(err)

    return err

def get_error_time_steps(t, y, y_spline, n_nodes, do_relative=True):

    """
    This function calculates the absolute or relative error comparing the approximate \
    solution and the "exact" solution at the end of each time step.

    :param numpy.ndarray t: the temporal nodes
    :param numpy.ndarray y: the approximate solution
    :param function y_spline: a function that may interpolate the "exact" (more accurate) solution
    :param int n_nodes: the number of nodes per time step
    :param bool do_relative: a flag indicating whether (if True) to calculate the relative error \
    or not to (if False) to calculate the absolute error

    :return: the absolute or relative error
    :rtype: numpy.ndarray
    """
    err = list()

    # 计算时间步数
    n_steps = (len(t) - 1) / (n_nodes - 1)
    n_steps = int(n_steps)

    for j in range(n_steps):

        # 获取时间步结束时的索引
        i = (n_nodes - 1) * j

        if do_relative:
            # 计算相对误差
            x = relative_norm( y_spline(t[i]) - y[i], y_spline(t[i]) )
        else:
            # 计算绝对误差
            x = LA.norm( y_spline(t[i]) - y[i], ord=np.inf )

        err.append(x)

    err = np.array(err)

    return err


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
    if top.ndim == 0:
        top_norm = np.abs(top)
    else:
        top_norm = LA.norm(top, axis=0)

    # 计算分母向量的范数
    if bot.ndim == 0:
        bot_norm = np.abs(bot)
    else:
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


def run_threshold(y, threshold):

    """
    If the value of the solution is below the threshold, set it to zero.

    :param numpy.ndarray y: the approximate solution
    :param float threshold: the threshold

    :return: None
    """
    # 将绝对值小于阈值的解分量设为0
    y[np.abs(y) < threshold] = 0

    return
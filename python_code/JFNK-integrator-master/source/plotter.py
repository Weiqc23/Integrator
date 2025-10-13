# This file was written by Dr. Namdi Brandon
# ORCID: 0000-0001-7050-1538
# June 29, 2018

"""
This module contains information about plotting.
"""

# ===============================================
# import
# ===============================================
import matplotlib.pylab as plt
import numpy as np
import numpy.linalg as LA

import integrator

# ===============================================
# function
# ===============================================
def corrections(d_sdc_norm, d_jfnk_norm, y_spect_norm, do_rerr, labels=None, do_legend=False):
    """
    此函数用于绘制谱延迟校正（SDC）方法和无雅可比牛顿 - 克里洛夫（JFNK）方法在每次迭代中的校正幅度。

    :param numpy.ndarray d_sdc_norm: 谱延迟校正（SDC）方法中近似值的最大范数 :math:`\\| \\delta^{[k]}_{sdc} \\|`
    :param numpy.ndarray d_jfnk_norm: 无雅可比牛顿 - 克里洛夫（JFNK）方法中近似值的最大范数 :math:`\\| \\delta^{[k]}_{jfnk} \\|`
    :param numpy.ndarray y_spect_norm: 谱解的范数 :math:`\\| y_{spect} \\|`，用于缩放校正值以计算相对误差
    :param bool do_rerr: 此标志指示是否绘制相对误差（若为 True）或绝对误差（若为 False）
    :param list labels: 线条的名称列表
    :param bool do_legend: 此标志指示是否显示图例

    :return: None
    """
    #
    # 绘制 SDC 和 JFNK 中的校正值
    #

    # 将 SDC 和 JFNK 的校正范数作为绘图数据
    data = (d_sdc_norm, d_jfnk_norm)

    # 设置 y 轴标签
    ylabel = 'log 10 Norm'

    if do_rerr:
        # 若绘制相对误差，将数据除以谱解的范数
        data = [x / y_spect_norm for x in data]
        # 设置图形主标题为相对校正
        main_title = 'Relative Corrections'
    else:
        # 若绘制绝对误差，设置图形主标题为绝对校正
        main_title = 'Absolute Corrections'

    # 设置子图标题
    titles = ['SDC', 'JFNK']
    # 设置 x 轴标签
    xlabel = 'Iteration'

    # 调用 plot_errors 函数进行绘图
    plot_errors(data, titles, labels, do_legend=do_legend, main_title=main_title, xlabel=xlabel, ylabel=ylabel, ls='-o')

    return

def plot_correction_time_steps(D, Y, do_save=False, fpath=None, do_close=False):

    """
    This function plots the magnitude of the corrections (both absolute and relative) \
    :math:`\\| \\delta \\|` vs iteration for ech time step.

    :param list D: the corrections over the simulation. List of length number of time steps, \
    containing a list of length number of iterations for the given time step, of the corrections \
    for the iteration. The corrections are of dimensions (number of time nodes, size of the problem)
    :param list Y: the corresponding approximations to the given correction.
    :param bool do_save: indicating whether or not to save the data
    :param str fpath: the file path in which to save the data
    :param bool do_close: a flag indicating whether or not to close the plots

    :return: None
    """

    # for each time step
    for k in range( len(D) ):

        # the list of the values per iteration
        d = D[k]
        y = Y[k]

        # for each iteration
        x = np.array( [ LA.norm(x, axis=0).max() for x in d ] )
        u = np.array( [ integrator.relative_norm(dd, yy) for dd, yy in zip(d, y) ] )

        # plot for each time step
        plt.figure(k)
        plt.title('step: ' + str(k))

        # plot the absolute correction
        plt.plot(range(len(x)), np.log10(x), '-o', label='|d|')

        # plot the relative correction
        plt.plot(range(len(u)), np.log10(u), '-o', label='|d|/|y|')

        plt.xlabel('iteration')
        plt.ylabel('log10(|d|)')
        plt.legend(loc='best')

        # save the figure
        if do_save:
            fname = fpath + ('\\step%.3d.pdf' % k)
            plt.savefig(fname, bbox_inches='tight')

        # close the figure
        if do_close:
            plt.close()

    return

def plot_errors(data, titles, labels, do_legend=False, main_title='', xlabel='', ylabel='', ls='-'):
    """
    此函数用于绘制每次迭代中 SDC 方法和 JFNK 方法的校正幅度。

    :param data: 模拟得到的误差数据
    :param list titles: 子图的标题列表
    :param list labels: 每条线的名称列表
    :param bool do_legend: 此标志指示是否显示图例
    :param str main_title: 整个图形的标题
    :param str xlabel: x 轴的标签
    :param str ylabel: y 轴的标签
    :param str ls: 线条样式

    :return: None
    """
    # 计算数据的数量，即子图的数量
    ndata = len(data)

    # 创建一个包含多个子图的画布，子图按列排列，且共享 y 轴
    fig, axes = plt.subplots(ncols=ndata, sharey=True)

    # 设置整个图形的标题
    fig.suptitle(main_title)

    # 用于存储所有线条对象的列表，以便后续添加图例
    lines = []

    # 如果没有提供标签列表，则为每个数据列生成一个 None 标签
    if labels is None:
        labels = [None] * len(data[0].T)

    # 遍历每个子图、对应的标题和数据
    for ax, title, d in zip(axes, titles, data):
        # 设置子图的标题
        ax.set_title(title)

        # 遍历数据的每一列和对应的标签
        for y, label in zip(d.T, labels):
            # 过滤掉非正值，将其替换为一个较小的正数，避免对数运算出错
            y_positive = np.where(y > 0, y, 1e-10)
            # 绘制对数曲线，并保存线条对象
            temp = ax.plot(np.log10(y_positive), ls, label=label)
            # 将线条对象添加到列表中
            lines.append(temp[0])
            # 设置子图的 x 轴标签
            ax.set_xlabel(xlabel)

    # 如果需要显示图例
    if do_legend:
        # 在图形的右上角添加图例
        fig.legend(lines, labels, loc='upper right')

    # 设置第一个子图的 y 轴标签
    ax = axes[0]
    ax.set_ylabel(ylabel)

    return
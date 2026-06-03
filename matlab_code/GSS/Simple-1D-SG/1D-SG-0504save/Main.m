%% 概述
% 使用电势和载流子浓度作为未知变量，使用 SG 格式，计算一维半导体稳态
% 只能求解一维半导体，具体的掺杂情况可以在 InitializeVariables.m 中设置
% 左侧电极偏置电压通过 SetBoundaryVoltage 来指定，右侧电极偏置电压为 0
% 电压爬升的复杂逻辑没实现，只是通过 SetBoundaryVoltage(0.1); NewtonSolve(); 来实现逐步抬升电极电压

%% 全局变量

% 无量纲化基值和物理常数
global L_base V_base E_base N_base J_base G_base t_base q NIE VT
% 器件相关参数
global P_doping N_doping device_length N
% 数值计算相关变量
global n_node dx
% 主要待求变量
global psi n p phi_n phi_p J_n J_p
% solver 相关
global J_mat

% 两个重要数值参数
global lambda2 delta2

%% 主要数据输入

% PN 结算例输入
P_doping = 1e17;    % P 区掺杂浓度
N_doping = 1e17;    % N 区掺杂浓度
n_node = 1e2;       % 有限差分法的总格点数目（包含两个边界点）
device_length = 30; % 器件长度，单位 um


%% 初始化 & Initial Guess

% 初始化各种常量和无量纲化系数
InitializeConstants();

% 初始化各个待求解向量，初始化 doping
InitializeVariables();

% 由 psi n p 计算出全部其他物理量
CalculateVariables();

% 实时绘图
RealTimePlot();


% 设置 MultiPrecision Computing Toolbox 的精度，需要先 addpath
% addpath('/Users/sj/Code/Semiconductor/MathLib/AdvanpixMCT-4.9.3.15018');
% mp.Digits(34);


% 测试将 lambda 赋为 0，从而求解 outer solution
% 注意，当 doping 过低时，好像 outer solution 的解也不太对
% lambda2 = lambda2 * 0;


%% 设置边界条件 & 正式迭代
SetBoundaryVoltage(0.1);
NewtonSolve();
SetBoundaryVoltage(0.2);
NewtonSolve();
SetBoundaryVoltage(0.3);
NewtonSolve();
SetBoundaryVoltage(0.4);
NewtonSolve();
SetBoundaryVoltage(0.5);
NewtonSolve();
SetBoundaryVoltage(0.6);
NewtonSolve();
SetBoundaryVoltage(0.7);
NewtonSolve();
SetBoundaryVoltage(0.8);
NewtonSolve();

% SetBoundaryVoltage(0.9);
% NewtonSolve();
% SetBoundaryVoltage(1.0);
% NewtonSolve();
% SetBoundaryVoltage(1.1);
% NewtonSolve();

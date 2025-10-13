clear all; format short e; clc; close all;
restoredefaultpath; %恢复默认路径
addpath(genpath('math\'));
addpath(genpath('Device\'));
addpath(genpath('other_tool\'));
addpath('problem_dependent\');

warning off;

%%
t0 = 0;
tfinal = 4e-2;
load('Parameters_MMC_Np4.mat');

tic;
sim('ModularMultiLevelConverter.slx');
toc;

fprintf('\n');
fprintf('\n');
get_param('ModularMultiLevelConverter', 'SolverType')       % 返回求解器类型（'Variable-step' 或 'Fixed-step'）
get_param('ModularMultiLevelConverter', 'SolverName')     % 返回实际使用的求解器名称
fprintf('Simulink总步数: %.d \n', length(tout)-1);
fprintf('\n');
fprintf('\n');

% [A,B,C,D,y0,electrical_states,inputs,outputs] = power_analyze('ModularMultiLevelConverter');
% 
%%
[y0,electrical_states,inputs,outputs] = setting_Problem();
% 变化的时间节点和变化后的值
[tseries_g_change, value_g_change] = GetSignalChange(tout, signal_list);

%% Collocation
fprintf('\n');
fprintf('\n');
fprintf('------------------------------------');
fprintf('\n');

h_Collocation = 0.8e-3;
fprintf('Collocation的设定步长为: %.8e 秒\n', h_Collocation);

global nodeType set_etol set_kmax
nodeType = 2; % (Gauss =1,  Radau-II = 2)
n = 4;
set_etol = 1e-5;%当set_etol设置得很大时，即直接在predicator函数里算一步TR，不执行后面的SDC或者GMRES_SDC过程     %1e-4; %1e-5
set_kmax = 100;

initialize_sw_state(y0, signal_list(1,:), n);
global m
[tplot_ss, yplot_ss, zplot_ss] = Collocation(n, m, h_Collocation, t0, y0, tfinal, tseries_g_change, value_g_change);
% build_mex_Collocation;
% [tplot_ss, yplot_ss, zplot_ss] = Collocation_mex(n, m, h_Collocation, t0, y0, tfinal, tseries_g_change, value_g_change);

width = 1200;
height = 600;
screenSize = get(0, 'ScreenSize');  % 获取屏幕分辨率
left = (screenSize(3) - width)/2;   % 水平居中
bottom = (screenSize(4) - height)/2; % 垂直居中
color_blue = [0 0.4470 0.7410];
color_orange = [0.8500 0.3250 0.0980];

set(0, 'DefaultAxesFontName', 'Times New Roman');

figure('Position', [left, bottom, width, height]);
plot(tout,Vout_sim, '-', ...
     'Color', 'k', ...
     'LineWidth', 1.5);
hold on;
scatter(tplot_ss, zplot_ss(:,1), ...
        40, ...                % marker 大小
        color_orange, ...               % marker 边缘颜色
        'x', ...               % marker 样式
        'LineWidth', 1.8);
% plot(tplot_ss, zplot_ss(:,1), '-x', ...
%      'Color', color_orange, ...
%      'LineWidth', 1.8);
grid on;
xlabel('time/(s)','FontName','Times New Roman','FontWeight','light');
ylabel('Vout/(V)','FontName','Times New Roman','FontWeight','light');
% xlim([0.12, 0.16]);
legend('Simulink','Collocation','FontName','Times New Roman','FontWeight','bold','Location','southeast');

%% SDC
fprintf('\n');
fprintf('\n');
fprintf('------------------------------------');
fprintf('\n');

global iformulation use_linear_interpolation
iformulation=0; % the formulation of the problem.
use_linear_interpolation = 1;

h_GMRES = h_Collocation;
fprintf('GMRES-SDC的设定步长为: %.8e 秒\n', h_GMRES);

global use_gmres ExImplicit 
ExImplicit = 1; % (Explicit)=0，(Implicit)=1
use_gmres = 0; % 使用gmres-sdc，取1；仅使用sdc，取0

global check_eig
check_eig = false;

global check_iter_times
check_iter_times = false;

global m
initialize_sw_state(y0, signal_list(1,:), n);
[tplot_gs,yplot_gs,zplot_gs] = GMRES_SDC(n, m, h_GMRES, t0, y0, tfinal, tseries_g_change, value_g_change);

width = 1200;
height = 600;
screenSize = get(0, 'ScreenSize');  % 获取屏幕分辨率
left = (screenSize(3) - width)/2;   % 水平居中
bottom = (screenSize(4) - height)/2; % 垂直居中
color_blue = [0 0.4470 0.7410];
color_orange = [0.8500 0.3250 0.0980];

set(0, 'DefaultAxesFontName', 'Times New Roman');

figure('Position', [left, bottom, width, height]);
plot(tout,Vout_sim, '-', ...
     'Color', 'k', ...
     'LineWidth', 1.5);
hold on;
scatter(tplot_gs,zplot_gs(:,1), ...
        36, ...                % marker 大小
        color_blue, ...               % marker 边缘颜色
        'o', ...               % marker 样式
        'LineWidth', 1.8);
% plot(tplot_TRBDF2,zplot_TRBDF2(:,1), '-o', ...
%      'Color', color_blue , ...
%      'LineWidth', 1.8);
grid on;
xlabel('time/(s)','FontName','Times New Roman','FontWeight','light');
ylabel('Vout/(V)','FontName','Times New Roman','FontWeight','light');
% xlim([0.12, 0.16]);
legend('Simulink','SDC(Implicit)','FontName','Times New Roman','FontWeight','bold','Location','southeast');

%%

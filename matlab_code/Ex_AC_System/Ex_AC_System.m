clear all; format short e; clc; close all;
restoredefaultpath; %恢复默认路径
addpath(genpath('math\'));
addpath(genpath('problem_dependent\'));
% 获取当前脚本所在目录
% current_folder = fileparts(mfilename('fullpath'));
% model_path = fullfile(current_folder, 'TransmissionLine_AC_System.slx');

% 检查文件是否存在
% if exist(model_path, 'file') == 4
%     sim(model_path);
% else
%     error('模型文件不存在: %s', model_path);
% end

warning off;

%%
t0 = 0;
tfinal = 0.1;
dt_sim = 1e-6; % simulink的步长
Line_num = 10; % 传输线段数

tic;
% sim('TransmissionLine_AC_System.slx');
sim('TransmissionLine_AC_System');
% 自动生成时间节点tout
toc;

fprintf('\n');
fprintf('\n');
get_param('TransmissionLine_AC_System', 'SolverType')       % 返回求解器类型（'Variable-step' 或 'Fixed-step'）
get_param('TransmissionLine_AC_System', 'SolverName')     % 返回实际使用的求解器名称
fprintf('Simulink总步数: %.d \n', length(tout)-1);
fprintf('\n');
fprintf('\n');

% [A,B,C,D,y0,electrical_states,inputs,outputs] = power_analyze('TransmissionLine_AC_System');

%%
[y0,electrical_states,inputs,outputs] = setting_Problem();

%% SDC
fprintf('\n');
fprintf('\n');
fprintf('------------------------------------');
fprintf('\n');

global iformulation
iformulation=0; % the formulation of the problem.

% max_h_GMRES = 15000/(-lambda) %稳定域步长限制
h_GMRES = 100e-6; % 0.01; %100e-6;
h = 100e-6;
% h_GMRES = (tfinal-t0)/2;
fprintf('GMRES-SDC的设定步长为: %.8e 秒\n', h_GMRES);

global nodeType use_gmres set_etol set_gtol set_kmax ExImplicit BE_TOL SDC_TOL N_ITER_MAX_NEWTON N_ITER_MAX_SDC
ExImplicit = 1; % S_wave采用前向欧拉格式时，取0；S_wave采用后向欧拉格式时，取1。
use_gmres = 0; % 使用gmres-sdc，取1；仅使用sdc，取0
nodeType = 2; % (用Gauss点，取1；用Radau-II点，取2)
n = 4; % 单个区间内的节点个数
set_etol = 1e-5;%当set_etol设置得很大时，即直接在predicator函数里算一步TR，不执行后面的SDC或者GMRES_SDC过程     %1e-4; %1e-5
set_gtol = 1e-6;
set_kmax = 100;% 最大迭代次数

SDC_TOL = 1e-5;
BE_TOL = 1e-12; % 向后欧拉法（猜测解）精度
N_ITER_MAX_NEWTON = 50; % 最大的JFNK迭代次数
N_ITER_MAX_SDC = 100; % 最大的SDC迭代次数

% ----------------------------------------------------------------1. GMRES_SDC--------------------------------------------------------
global m

[tplot,yplot,zplot] = GMRES_SDC(n, m, h_GMRES, t0, y0, tfinal);

width = 1200;
height = 600;
screenSize = get(0, 'ScreenSize');  % 获取屏幕分辨率
left = (screenSize(3) - width)/2;   % 水平居中
bottom = (screenSize(4) - height)/2; % 垂直居中

set(0, 'DefaultAxesFontName', 'Times New Roman');

figure('Position', [left, bottom, width, height]);
subplot(2,1,1);
plot(tout,u_Load, '-', ...
     'Color', 'k', ...
     'LineWidth', 1.5);
hold on;
plot(tplot,zplot(:,1), '--', ...
     'Color', 'b', ...
     'LineWidth', 1.8);
grid on;
xlabel('time/(s)','FontName','Times New Roman','FontWeight','light');
ylabel('u Load/(V)','FontName','Times New Roman','FontWeight','light');
% xlim([0.15, 0.16]);
legend('Simulink','GMRES_SDC','FontName','Times New Roman','FontWeight','bold','Location','southeast');

subplot(2,1,2);
plot(tout,i_Load, '-', ...
     'Color', 'k', ...
     'LineWidth', 1.5);
hold on;
plot(tplot,zplot(:,2), '--', ...
     'Color', 'b', ...
     'LineWidth', 1.8);
grid on;
xlabel('time/(s)','FontName','Times New Roman','FontWeight','light');
ylabel('i Load/(A)','FontName','Times New Roman','FontWeight','light');
% xlim([0.03+1e-5, 0.03+5e-5]);
% legend('Simulink','Proposed Method','FontName','Times New Roman','FontWeight','bold','Location','southeast');
legend('Simulink','GMRES_SDC','FontName','Times New Roman','FontWeight','bold','Location','southeast');


% ----------------------------------------------------------------2. Collocation--------------------------------------------------------
%% Collocation
fprintf('\n');
fprintf('\n');
fprintf('------------------------------------');
fprintf('\n');

h_Collocation = 100e-6;
fprintf('Collocation的设定步长为: %.8e 秒\n', h_Collocation);
[tplot_ss, yplot_ss, zplot_ss] = Collocation(n, m, h_Collocation, t0, y0, tfinal);

width = 1200;
height = 600;
screenSize = get(0, 'ScreenSize');  % 获取屏幕分辨率
left = (screenSize(3) - width)/2;   % 水平居中
bottom = (screenSize(4) - height)/2; % 垂直居中

set(0, 'DefaultAxesFontName', 'Times New Roman');

figure('Position', [left, bottom, width, height]);
subplot(2,1,1);
plot(tout,u_Load, '-', ...
     'Color', 'k', ...
     'LineWidth', 1.5);
hold on;
plot(tplot_ss,zplot_ss(:,1), '--', ...
     'Color', 'b', ...
     'LineWidth', 1.8);
grid on;
xlabel('time/(s)','FontName','Times New Roman','FontWeight','light');
ylabel('u Load/(V)','FontName','Times New Roman','FontWeight','light');
% xlim([0.03+1e-5, 0.03+5e-5]);
legend('Simulink','Collocation','FontName','Times New Roman','FontWeight','bold','Location','southeast');

subplot(2,1,2);
plot(tout,i_Load, '-', ...
     'Color', 'k', ...
     'LineWidth', 1.5);
hold on;
plot(tplot_ss,zplot_ss(:,2), '--', ...
     'Color', 'b', ...
     'LineWidth', 1.8);
grid on;
xlabel('time/(s)','FontName','Times New Roman','FontWeight','light');
ylabel('i Load/(A)','FontName','Times New Roman','FontWeight','light');
% xlim([0.03+1e-5, 0.03+5e-5]);
legend('Simulink','Collocation','FontName','Times New Roman','FontWeight','bold','Location','southeast');


%%
% ---------------------------------------------------------------3. JFNK--------------------------------------------------------

fprintf('\n');
fprintf('\n');
fprintf('------------------------------------');
fprintf('\n');

p = Params(n,m,nodeType);
[tplot_ssss, yplot_ssss, zplot_ssss] = JFNK_adaptive(@f_eval,t0,tfinal,h,p,y0,N_ITER_MAX_NEWTON,N_ITER_MAX_SDC,BE_TOL,SDC_TOL);

% 1. 图窗宽度 1200 像素
width = 1200;

% 2. 图窗高度 600 像素
height = 600;

% 3. 读取屏幕整屏尺寸 [left, bottom, width, height]
screenSize = get(0, 'ScreenSize');

% 4. 计算窗口左上角 x 坐标，使窗口水平居中
left = (screenSize(3) - width) / 2;

% 5. 计算窗口左上角 y 坐标，使窗口垂直居中
bottom = (screenSize(4) - height) / 2;

% 6. 全局设置：坐标轴文字默认用 Times New Roman 字体
set(0, 'DefaultAxesFontName', 'Times New Roman');

% 7. 新建图窗并指定位置和大小（居中 1200×600）
figure('Position', [left, bottom, width, height]);

% 8. 把图窗分成 2×1 格子，选中第 1 格（上面子图）
subplot(2,1,1);

% 9-11. 用黑色实线画 Simulink 结果：负载电压 u_Load
plot(tout, u_Load, '-', 'Color', 'k', 'LineWidth', 1.5);

% 12. 保持当前坐标轴，后续 plot 不擦除
hold on;

% 13-15. 用蓝色虚线画自己算法结果：负载电压 zplot 第 1 列
plot(tplot_ssss, zplot_ssss(:,1), '--', 'Color', 'b', 'LineWidth', 1.8);

% 16. 打开网格线
grid on;

% 17. 设置横轴标签：time/(s)，字体 Times New Roman，细体
xlabel('time/(s)', 'FontName', 'Times New Roman', 'FontWeight', 'light');

% 18. 设置纵轴标签：u Load/(V)，字体同上
ylabel('u Load/(V)', 'FontName', 'Times New Roman', 'FontWeight', 'light');

% 19. 暂时注释掉：若只想看局部时间可打开 xlim
% xlim([0.15, 0.16]);

% 20. 图例：黑色实线=Simulink，蓝色虚线=Proposed Method，放右下角
legend('Simulink', 'JFNK', 'FontName', 'Times New Roman','FontWeight', 'bold', 'Location', 'southeast');

% 21. 选中第 2 格（下面子图）
subplot(2,1,2);

% 22-24. 用黑色实线画 Simulink 负载电流 i_Load
plot(tout, i_Load, '-', 'Color', 'k', 'LineWidth', 1.5);

% 25. 保持坐标轴
hold on;

% 26-28. 用蓝色虚线画自己算法结果：负载电流 zplot 第 2 列
plot(tplot_ssss, zplot_ssss(:,2), '--', 'Color', 'b', 'LineWidth', 1.8);

% 29. 打开网格
grid on;

% 30-31. 横轴标签同上
xlabel('time/(s)', 'FontName', 'Times New Roman', 'FontWeight', 'light');

% 32-33. 纵轴标签：i Load/(A)
ylabel('i Load/(A)', 'FontName', 'Times New Roman', 'FontWeight', 'light');

% 34. 局部放大注释掉
% xlim([0.03+1e-5, 0.03+5e-5]);

% 35. 图例同上
legend('Simulink','JFNK','FontName','Times New Roman','FontWeight','bold','Location','southeast');


%% 比较 yplot_ssss 和 yplot 的误差
fprintf('\n====================================\n');
% [n_rows, n_cols] = size(yplot);
% for (i = 1:n_rows)
%     fprintf("i: %d\n",i);
%     disp("sdc_gemras ans:");
%     disp(num2str(yplot(i,:), '%.12f '));
%     disp("sdc_gemras_sdc ans:");
%     disp(num2str(yplot_sdc(i,:), '%.12f '));
%     errest = max(abs(yplot(i,:)-yplot_sdc(i,:)));
%     fprintf("sdc_gemras方法和sdc_gemras_sdc方法最大误差: %.12f\n",errest);
%     disp("jfnk ans:");
%     disp(num2str(yplot_ssss(i,:), '%.12f '));
%     errest = max(abs(yplot_ssss(i,:)-yplot(i,:)));
%     fprintf("sdc_gemras方法和jfnk方法最大误差: %.12f\n",errest);
% end


disp('sdc_gemras最后一行');
disp(num2str(yplot(end,:), '%.12f '));
disp('sdc_gemras倒数第二行');
disp(num2str(yplot(end-1,:), '%.12f '));

disp('jfnk最后一行');
disp(num2str(yplot_ssss(end,:), '%.12f '));
disp('jfnk倒数第二行');
disp(num2str(yplot_ssss(end-1,:), '%.12f '));

err_last = yplot(end,:) - yplot_ssss(end,:);
err_second_last = yplot(end-1,:) - yplot_ssss(end-1,:);

disp('误差end');
disp(num2str(err_last, '%.12e '));   % 用科学计数法更清晰

disp('误差end-1');
disp(num2str(err_second_last, '%.12e '));


% [nRows, nCols] = size(yplot);
% last_row = yplot(end, :);
% second_last_row = yplot(end-1, :);
% output_cell = cell(7, nCols + 1);
% output_cell{1, 1} = sprintf('yplot形状: %d x %d', nRows, nCols);
% output_cell(2, 2:end) = num2cell(last_row);
% output_cell(3, 2:end) = num2cell(second_last_row);
% output_cell{2, 1} = 'yplot(end)';
% output_cell{3, 1} = 'yplot(end-1)';
% 
% last_row_ssss = yplot_ssss(end, :);
% second_last_row_ssss = yplot_ssss(end-1, :);
% 
% output_cell{5, 1} = sprintf('yplot形状: %d x %d', nRows, nCols);
% output_cell(6, 2:end) = num2cell(last_row_ssss);
% output_cell(7, 2:end) = num2cell(second_last_row_ssss);
% output_cell{6, 1} = 'yplot_ssss(end)';
% output_cell{7, 1} = 'yplot_ssss(end-1)';
% writecell(output_cell, 'gemras_sdc和jfnk对比.xlsx');

% abs_error = abs(yplot_ssss' - yplot');
% % rel_error = abs_error ./ abs(yplot_ssss);
% fprintf("sdc_gemras方法和jfnk方法绝对误差: %.12f\n",max(abs_error));
% 
% fprintf("\n\n")
% fprintf('sdc_gemras方法和jfnk方法相对误差: %.12f\n', max(max(rel_error)));


%% run_setting.m
% 运行 setting_Problem 函数，补全所需变量，并输出系统矩阵

clc; clear; close all;

% ========== 补全模型需要的变量 ==========

Line_num = 10;
t0 = 0;
tfinal = 0.1;
dt_sim = 1e-6; 
Line_num = 10;
 

[y0, electrical_states, inputs, outputs] = setting_Problem();

disp('初始状态向量 y0 = ');
disp(y0);

disp('电气状态变量:');
disp(electrical_states);

disp('输入变量:');
disp(inputs);

disp('输出变量:');
disp(outputs);

fprintf('系统矩阵 A, B, C, D 已导出到 CSV 文件\n');

% 已在 setting_Problem.m 里保存 A.csv, B.csv, C.csv, D.csv, y0.csv, inputs.txt

function [y0,electrical_states,inputs,outputs] = setting_Problem()
%% 问题描述
% 电路的状态空间方程可以抽象为如下的ODE系统:
% dy/dt = Ay + Bu;
% z = Cx + Du;
% y 为系统中独立的状态变量(electrical_states)；z为非状态变量
% u=u(t)，在线性系统下为已知函数

%%
global A B C D m mi mo
% 调用simulink的Power_analyze函数，可以自动获得方程矩阵，以及变量y，u的排列顺序
% [A,B,C,D,y0,electrical_states,inputs,outputs] = power_analyze('TransmissionLine_AC_System');
addpath('C:\Users\wei13\Desktop\Basic_Code\Basic_Code\Ex_AC_System');
fprintf('start power analyze\n');
try
    [A,B,C,D,y0,electrical_states,inputs,outputs] = power_analyze('TransmissionLine_AC_System');
    % electrical_states：电气状态信息（状态变量名称、位置、类型（电容电压&电感电流）、数值）
    % inputs：输入变量信息（输入信号名称，类型（电压源&电流源），位置，波形）
    % outputs：输出变量信息（输出信号名称，测量点，物理量（电压、电流、功率））
catch ME
    fprintf('=== Exception message ===\n%s\n\n', ME.message);
    fprintf('=== getReport(ME,''extended'') ===\n');
    disp(getReport(ME,'extended'));   % 完整错误 + 堆栈
    % - 错误类型
    fprintf('\n=== ME.stack items ===\n');
    for k = 1:numel(ME.stack)
        fprintf('  %d: %s  (line %d)  in %s\n', k, ME.stack(k).name, ME.stack(k).line, ME.stack(k).file);
        % - k: 堆栈层级（从最内层开始）
        % - ME.stack(k).name: 函数名
        % - ME.stack(k).line: 出错行号
        % - ME.stack(k).file: 文件名
    end
end
m = length(A);
% 维度
mi = size(D,2);
% D的列数：输入变量数
mo = size(D,1);
% D的行数：输出变量数

y0 = y0';

% fprintf("output to csv")
% writematrix(A, "A.csv");
% writematrix(B, "B.csv");
% writematrix(y0(:), "y0.csv");

% fid = fopen("inputs.txt","w");
% for i = 1:length(inputs)
%     fprintf(fid, "%s\n", inputs{i});
% end
% fclose(fid);
% 
% writematrix(C, "C.csv");
% writematrix(D, "D.csv");

end
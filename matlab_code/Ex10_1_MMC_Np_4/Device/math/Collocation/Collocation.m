function [tplot, yplot, zplot] = Collocation(n, m, h0, t0, y0, tfinal, tseries_g_change, value_g_change)
global nodeType
if (nodeType==1)% Gauss
    [A,B,tc]=chebnodes(n);  %set up the integration matrix, and construct tc.
else %Radau-II
    [A,tc]=tlgr(n);  %set up the integration matrix, and construct tc.
    B = [];
end

global set_etol
etol = set_etol;

global set_kmax
kmax = set_kmax;

%%
% 变换次数
num_g_change = size(value_g_change,1);
num_pre_guess = max(num_g_change, round((tfinal-t0)/h0)+1)*(n+1);
yplot = zeros(num_pre_guess,m);
global mo
zplot = zeros(num_pre_guess,mo);
tplot = zeros(num_pre_guess,1);

tplot(1) = t0;
yplot(1, :) = y0;
zplot(1, :) = GetOutputs_SingleStep(t0, y0);
plot_count = 1;
nSteps = 0;
%%
% 初始化全局计时器变量
% global time_rhs time_dfdy time_solve time_switch
% time_rhs = 0;
% time_dfdy = 0;
% time_solve = 0;
% time_switch = 0;

time_Collocation_onestep = 0;
time_save_results = 0;

ynow = y0;

% 改为更稳的计时方式
startTime = datetime('now');
%%

% get_signal = true
global get_signal
for i = 1 : num_g_change
    
    tic;

    tnow = tseries_g_change(i);
    tfinal = tseries_g_change(i+1);
    dt = min(h0, tfinal-tnow);

    if get_signal 
        signal_list = value_g_change(i,:);
        GetSwitchStatusFromSignal(tnow, ynow, signal_list);
    end

    time_Collocation_onestep = time_Collocation_onestep + toc;

    while tnow < tfinal

        nSteps = nSteps+1;

        if dt > tfinal - tnow
            dt = tfinal - tnow;
        end


        tic;
        [tnow,ynow, tnode,ynode,znode] = Collocation_onestep(m, n, tnow, ynow, dt, tc, A, B, kmax, etol); %m,y0,dt,n,A
        time_Collocation_onestep = time_Collocation_onestep + toc;

        tic;
            numNew = length(tnode);
            % 填充预分配的数组
            itv = plot_count+1:plot_count+numNew;
            tplot(itv) = tnode;
            yplot(itv, :) = ynode;
            zplot(itv, :) = znode;

            plot_count = plot_count + numNew;
            
        time_save_results = time_save_results + toc;

    end

end

%%
num_step = find(tplot==0,2);
if (length(num_step)>1)
    plot_interval = (num_step(1):num_step(2)-1)';
    tplot = tplot(plot_interval);
    yplot = yplot(plot_interval,:);
    zplot = zplot(plot_interval,:);
end

endTime = datetime('now');
elapsedTime = seconds(endTime - startTime);

% 打印结果
fprintf('Collocation() Total Execution Time: %.4f s\n', elapsedTime);
fprintf('Collocation() Total Execution Steps: %d \n', int32(nSteps));

% fprintf('--- 各部分累计时间 ---\n');
% fprintf('  - rhs 部分：    %.4f 秒\n', time_rhs);
% fprintf('  - dfdy 部分：   %.4f 秒\n', time_dfdy);
% fprintf('  - 求解部分：    %.4f 秒\n', time_solve);
% fprintf('  - switch 部分： %.4f 秒\n', time_switch);
% 
% sumParts = time_rhs + time_dfdy + time_solve + time_switch;
% fprintf('各部分耗时总和： %.4f 秒\n', sumParts);
% fprintf('差值（总时间 - 各部分之和）： %.6f 秒\n', elapsedTime - sumParts);

fprintf('--- The cumulative time of each part ---\n');
fprintf('  - Collocation_onestep():  %.4f s\n', time_Collocation_onestep);
fprintf('  - save_results():  %.4f s\n', time_save_results);


end

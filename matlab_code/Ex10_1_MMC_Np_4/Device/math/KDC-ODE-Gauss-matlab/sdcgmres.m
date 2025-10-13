function [tplot,yplot,zplot]=sdcgmres(iprob,m,t0,tfinal,y0,h0,n,kmax,gtol,etol,k0, tseries_g_change, value_g_change)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INTENT (out)
%   ysolfin (row vector) solution y at final time
%   res (row vector) residues outputed by gmres
%   indres (row vector) index for res : res(indres(i-1)+1:indres(i)), 
%          residues outputs from gmres' during the i_th time step
%        
%   err (row vector) error after each gmres
%   inderr (row vector) index for err: err(inderr(i-1)+1:inderr(i)),
%          errors outputs from gmres' during the i_th time step
%   iter (row vector) iteration numbers outputed by each gmres.
% 
% INTENT (in) 
%   Most are scalars.
%   m size of the problem
%   t0 initial time
%   tfinal fintal time
%   y0 (row vector) initial value for y
%   h0 step size
%   n number of grid points used for each time step
%   kmax maximal number that the gmres is done
%   gtol tolerance for gmres
%   etol tolerance for err/res
%   k0 gmresk selecting parameter
%
% Last change: Jingfang Huang, 03/10/2005.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initiallization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global nodeType
if (nodeType==1)% Gauss
    [As,Bs,tcs]=chebnodes(n);  %set up the integration matrix, and construct tc.
    A = zeros(n+1); A(2:end, 2:end) = As;
    B = zeros(1,n+1); B(2:end) = Bs;
    tc = zeros(1,n+1); tc(2:end) = tcs;
else %Radau-II
    [As,tcs]=tlgr(n);  %set up the integration matrix, and construct tc.
    A = zeros(n+1); A(2:end, 2:end) = As;
    tc = zeros(1,n+1); tc(2:end) = tcs;
    B = [];
end


tnow=t0; dt=h0; ynow=y0; % Now set up the initial values for iteration.
znow = GetOutputs_SingleStep(tnow, ynow);
count=0; %monitor values.

% 开关信号改变的次数
num_g_change = size(value_g_change,1);
guess_coeff = 2; %被动开关会导致步数增加，留的裕量
num_pre_guess = guess_coeff*round( max(num_g_change, (tfinal-t0)/dt+1)*(n+1)) ;
yplot = zeros(num_pre_guess,m);
global mo
zplot = zeros(num_pre_guess,mo);
tplot = zeros(num_pre_guess,1);

tplot(1) = tnow;
yplot(1, :) = ynow;
zplot(1, :) = znow;
plot_count = 1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   the main marching scheme. No adaptive steps yet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_GMRES_onestep = 0;
time_save_results = 0;
global count_atv count_updated
count_atv = 0;
count_updated = 0;

% 改为更稳的计时方式
startTime = datetime('now');

% global if_save_state

global check_iter_times
if (check_iter_times)
    tplot_whole = t0;
    count_updated_plot =[];
end

global check_eig
max_eig_save = 0;
if (check_eig)
    dfdy=jeval(m,ynow,tnow);
    L = (dfdy(:,:,1))';
    max_eig_save = -max(abs(real(eig(L))));
end

%%
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

    time_GMRES_onestep = time_GMRES_onestep + toc;

    while tnow<tfinal
        count=count+1;

        if dt>tfinal-tnow
            dt=tfinal-tnow;  % find the right time step. This is important for last step.
        end

        if (dt<1e-10)
            break;
        end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % march one-step. main code.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        [tnow, ynow, znow, tnode,ynode, znode, iters]=...
            onestep(m,tnow,ynow,dt,n,tc,kmax,gtol,etol,k0,A,B);
        time_GMRES_onestep = time_GMRES_onestep + toc;

        if (check_eig)
            if (nodeType==1)% Gauss
                vec_t = [tnode(2:end), tnow];
                vec_y = [ynode(2:end,:); ynow];
            else % Radau-II
                vec_t = tnode(2:end);
                vec_y = ynode(2:end,:);
            end
            dfdy=jeval(m,vec_y,vec_t);
            for j = 1 : length(vec_t)
                L = (dfdy(:,:,j))';
                max_eig_save = [max_eig_save; -max(abs(real(eig(L))))];
            end
        end

        if (check_iter_times)
            tplot_whole = [tplot_whole; tnow];
            count_updated_plot = [count_updated_plot; iters];
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  monitor values.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if (ierrmsg==1)
        %     sprintf('Bad Results, no convergence, stopped');
        %     break
        % end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % record errors for each iteration.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;

        % 换用一种更高效的方式:
        if (nodeType==1) %Gauss
            % 当前新增的数据量
            numNew = length(tnode);
            % 填充预分配的数组
            itv = plot_count+1:plot_count+numNew;
            tplot(itv) = [tnode(2:end), tnow]';
            yplot(itv, :) = [ynode(2:end,:); ynow];
            zplot(itv, :) = [znode(2:end,:); znow];
        else %Radau-II
            % 当前新增的数据量
            numNew = length(tnode) - 1;
            % 填充预分配的数组
            itv = plot_count+1:plot_count+numNew;
            tplot(itv) = tnode(2:end)';
            yplot(itv, :) = ynode(2:end, :);
            zplot(itv, :) = znode(2:end, :);
        end
        plot_count = plot_count + numNew;

        time_save_results = time_save_results + toc;

       % if (if_save_state)
       %     Save_State(tnode);
       % end


    end

end

%%
% tplot(plot_count+1) = tfinal;
% yplot(plot_count+1,:) = ynow;
% zplot(plot_count+1,:) = znow;

num_step = find(tplot==0,2);
if (length(num_step)>1)
    plot_interval = (num_step(1):num_step(2)-1)';
    tplot = tplot(plot_interval);
    yplot = yplot(plot_interval,:);
    zplot = zplot(plot_interval,:);
end

endTime = datetime('now');
elapsedTime = seconds(endTime - startTime);

fprintf('GMRES-SDC Total Execution Time: %.4f s\n', elapsedTime);

fprintf('--- The cumulative time of each part ---\n');
fprintf('  - onestep():  %.4f s\n', time_GMRES_onestep);
fprintf('  - save_result():  %.4f s\n', time_save_results);

fprintf('GMRES-SDC Total Execution Steps %d (whole step)\n', int32(count));
fprintf('  - updated() Num of Execution:  %.0f \n', count_updated);
fprintf('  - atv() Num of Execution (or "(I-dt*S)*" Num of Execution):  %.0f \n', count_atv);
% fprintf('  - 牛顿迭代的总执行次数）：    %.0f 次\n', length(iter));

if (check_eig) % 所有时刻的特征值画图
    % 0. 使用 semilogy 绘制实部随时间的变化. 核心步骤：绘制绝对值并修改坐标轴标签 ---
    figure;
    subplot(2,1,1);
    % 步骤 0a: 同样使用 semilogy 绘制绝对值
    semilogy(tplot, abs(real(max_eig_save)), 'b-', 'LineWidth', 2);
    % 步骤 0b: 获取坐标轴句柄并设置负数标签 (同之前的方法)
    ax = gca;
    yticks_positive = ax.YTick;
    yticklabels_negative = arrayfun(@(x) sprintf('-10^{%g}', log10(x)), yticks_positive, 'UniformOutput', false);
    ax.YTickLabel = yticklabels_negative;
    % -----------------------------------------------------------------
    % 步骤 0c: 【关键步骤】反转 Y 轴的方向
    set(gca, 'YDir', 'reverse');
    % -----------------------------------------------------------------
    % 步骤 0d: 添加图形注释
    title('max\_eig\_save');
    xlabel('time (t)');
    ylabel('RealParts (log)');
    grid on;

    % % 1. 绘制所有不重复的点
    % % 使用 unique 函数获取向量中的唯一值
    % unique_points = unique(max_eig_save);
    % 
    % figure; % 创建一个新的图形窗口
    % plot(real(unique_points), imag(unique_points), 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    % title('max\_eig\_save');
    % xlabel('Real Part');
    % ylabel('Imaginary Part');
    % grid on;
    % axis equal; % 使x轴和y轴的刻度长度相等，避免图形畸变

    % 2. 找到实部绝对值最大的数并打印
    % 计算每个元素的实部的绝对值
    abs_real_parts = abs(real(max_eig_save));
    [~, index_max] = max(abs_real_parts);
    max_abs_real_num = max_eig_save(index_max); % 获取该索引对应的原始复数
    fprintf('The largest real part of the eigenvalue is: %.4e + (%.4e)i\n', real(max_abs_real_num), imag(max_abs_real_num));

    % 3. 找到实部绝对值最小的数并打印
    % 计算每个元素的实部的绝对值
    [~, index_min] = min(abs_real_parts);
    min_abs_real_num = max_eig_save(index_min); % 获取该索引对应的原始复数
    fprintf('The smallest real part of the eigenvalue is: %.4e + (%.4e)i\n', real(min_abs_real_num), imag(min_abs_real_num));

end

if (check_iter_times)

subplot(2,1,2);
plot(tplot_whole(2:end), count_updated_plot, '-o', ...
     'LineWidth', 0.5, ...    % 线条粗细
     'MarkerSize', 1, ...     % 圆圈大小
     'MarkerFaceColor', 'b'); % 圆圈内部填充颜色（可选）
grid on;

end


return
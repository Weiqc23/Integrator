function [tplot,yplot,zplot]=sdcgmres(m,t0,tfinal,y0,h0,n,kmax,gtol,etol,k0)
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
    [A,B,tc]=chebnodes(n);  %set up the integration matrix, and construct tc.
else %Radau-II
    [A,tc]=tlgr(n);  %set up the integration matrix, and construct tc.
    B = zeros(n);
end


tnow=t0; dt=h0; ynow=y0; % Now set up the initial values for iteration.
znow = GetOutputs_SingleStep(tnow, ynow);
count=0; %monitor values.


yplot = ynow;% y0
zplot = znow;
tplot = tnow;% t0

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

global check_iter_times
if (check_iter_times)% 记录每一时步SDC（牛顿迭代）次数
    tplot_whole = t0;
    count_updated_plot =[];
end


%%
    while tnow<tfinal
        count=count+1;
        v0 = ynow;
        if dt>tfinal-tnow
            dt=tfinal-tnow;  % find the right time step. This is important for last step.
        end
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % march one-step. main code.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        % tnow ~ tnow + dt sdc求解
        [tnow, ynow, znow, tnode, ynode, znode, iters]=...
            onestep(m,tnow,ynow,dt,n,tc,kmax,gtol,etol,k0,A,B);
            
        time_GMRES_onestep = time_GMRES_onestep + toc;

        if (check_iter_times) % 记录每一时步SDC（牛顿迭代）次数
            tplot_whole = [tplot_whole; tnow];
            count_updated_plot = [count_updated_plot; iters];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save solution.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;

        if (nodeType==1) %Gauss
            % 当前新增的数据量
            vec_t = [tnode(2:end), tnow]';
            vec_y = [ynode(2:end,:); ynow];
            vec_z = [znode(2:end,:); znow];
        else %Radau-II
            % 当前新增的数据量
            vec_t = tnode(:)';
            vec_y = ynode(:, :);
            vec_z = znode(:, :);
        end

        tplot = [tplot, vec_t];
        yplot = [yplot; vec_y];
        zplot = [zplot; vec_z];

        time_save_results = time_save_results + toc;

    end


%%

endTime = datetime('now');
elapsedTime = seconds(endTime - startTime);

fprintf('GMRES-SDC Total Execution Time: %.4f s\n', elapsedTime);

fprintf('--- The cumulative time of each part ---\n');
fprintf('  - onestep():  %.4f s\n', time_GMRES_onestep);
fprintf('  - save_result():  %.4f s\n', time_save_results);

fprintf('GMRES-SDC Total Execution Steps %d (whole step)\n', int32(count));
fprintf('  - updated() Num of Execution:  %.0f \n', count_updated);
fprintf('  - atv() Num of Execution (or "(I-dt*S)*" Num of Execution):  %.0f \n', count_atv);


if (check_iter_times)

    subplot(2,1,2);
    plot(tplot_whole(2:end), count_updated_plot, '-o', ...
         'LineWidth', 0.5, ...    % 线条粗细
         'MarkerSize', 1, ...     % 圆圈大小
         'MarkerFaceColor', 'b'); % 圆圈内部填充颜色（可选）
    grid on;

end


return
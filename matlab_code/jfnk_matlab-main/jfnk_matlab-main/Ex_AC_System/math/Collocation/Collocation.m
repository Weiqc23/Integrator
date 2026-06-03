function [tplot, yplot, zplot] = Collocation(n, m, dt, t0, y0, tfinal)
global nodeType
if (nodeType==1)% Gauss
    [A,B,tc]=chebnodes(n);  %set up the integration matrix, and construct tc.
else %Radau-II
    [A,tc]=tlgr(n);  %set up the integration matrix, and construct tc.
    B = zeros(n);
end
%%
tnow = t0;
ynow = y0;
znow = GetOutputs_SingleStep(tnow, ynow);

yplot = ynow;
zplot = znow;
tplot = tnow;

nSteps = 0;
%%
% 初始化全局计时器变量
time_Collocation_onestep = 0;
time_save_results = 0;

% 改为更稳的计时方式
startTime = datetime('now');
%%

    while tnow < tfinal

        nSteps = nSteps+1;

        if dt > tfinal - tnow
            dt = tfinal - tnow;
        end


        tic;
        [tnow,ynow,znow, tnode,ynode,znode] = Collocation_onestep(m, tnow, ynow, dt, n, tc, A); %m,y0,dt,n,A
        time_Collocation_onestep = time_Collocation_onestep + toc;

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
            vec_t = tnode(1:end)';
            vec_y = ynode(1:end, :);
            vec_z = znode(1:end, :);
        end

        tplot = [tplot; vec_t];
        yplot = [yplot; vec_y];
        zplot = [zplot; vec_z];

        time_save_results = time_save_results + toc;

    end


%%
endTime = datetime('now');
elapsedTime = seconds(endTime - startTime);

% 打印结果
fprintf('Collocation() Total Execution Time: %.4f s\n', elapsedTime);
fprintf('Collocation() Total Execution Steps: %d \n', int32(nSteps));
fprintf('--- The cumulative time of each part ---\n');
fprintf('  - Collocation_onestep():  %.4f s\n', time_Collocation_onestep);
fprintf('  - save_results():  %.4f s\n', time_save_results);


end

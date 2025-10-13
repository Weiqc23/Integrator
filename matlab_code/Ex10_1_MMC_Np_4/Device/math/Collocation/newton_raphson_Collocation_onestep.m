function ynode = newton_raphson_Collocation_onestep(m, y0, tnode, dt, p, S, kmax, etol)
% 行向量风格的一步 Collocation + Newton 求解器
% Inputs:
%   m  - 状态维数
%   y0 - 1 x m 行向量，步起始值
%   t  - p x 1 节点位置（相对位置，举例在 [0,1]）
%   dt - 时间步长
%   p  - 节点数
%   S  - p x p 谱积分矩阵
% Output:
%   y_solution - p x m 矩阵（第 i 行为节点 i 的解）

    iter = 0;
    converged = false;

    ynode = repmat(y0, p, 1);    % p x m, 初始猜测值

    %% Newton 迭代主循环
    while ~converged && iter < kmax

        % ------- 1) 在每个节点评估 f-------
        f = fright(m,ynode,tnode);


        % ------- 2) 计算残差 R = Y - y0 - dt * S * F -------
        R = ynode - repmat(y0, p, 1) - dt * (S * f); % p x m
        res_norm = max(abs(R(:)));
        if res_norm < etol
            converged = true;
            break;
        end


        % ------- 3) 在每个节点评估f 并组装大雅可比矩阵 Jbig （pm x pm） -------
        %  J = A + BSW * PvPxtrue
        dfdy=jeval(m,ynode,tnode);  % evaluate the jacobian matrix.

        pm = p * m;
        Jbig = zeros(pm, pm);  % 若 p*m 很大，建议用 sparse

        for i = 1:p
            row_idx = (i-1)*m + (1:m);    % 大矩阵中第 i 块行的索引
            for j = 1:p
                col_idx = (j-1)*m + (1:m); % 第 j 块列的索引

                % 注意这里要用 dfdy 的转置：见后面解释
                % 这里转置是因为jeval返回的是已经转置一次的
                blk = - dt * S(i, j) * (dfdy(:,:,j)'); % m x m

                % if i == j
                %     blk = blk + eye(m);
                % end
                Jbig(row_idx, col_idx) = blk;
            end
        end

        Jbig = Jbig+eye(pm);

        % ------- 4) 求解线性增量方程 Jbig * dYvec = -Rvec -------
        Rvec = reshape(R', pm, 1);       % 按行堆叠成列向量（vec of R'）
        dYvec = - (Jbig \ Rvec);         % 求解（直接解法，视规模可换 GMRES）
        dY = reshape(dYvec, m, p)';      % 恢复为 p x m

        % ------- 5) 更新与收敛判断 -------
        ynode = ynode + dY;
        iter = iter + 1;
     
    end

    if ~converged
        fprintf('Warning: Newton did not converge in %d iters, residual=%g\n', kmax, res_norm);
    end

end

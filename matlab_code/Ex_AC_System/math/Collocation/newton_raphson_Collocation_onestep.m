function y_solution = newton_raphson_Collocation_onestep(m, y0, t, dt, p, S)
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

global set_etol set_kmax

    %% 参数与初值
    tol = set_etol;
    max_iter = set_kmax;

    % 初始猜测：把每个节点都设为起始值 y0（行向量）
    Y = repmat(y0, p, 1);    % p x m

    iter = 0;
    converged = false;

    %% Newton 迭代主循环
    while ~converged && iter < max_iter

        % ------- 1) 在每个节点评估 f 和雅可比 dfdy -------
        [f, dfdy] = Get_fright_and_jeval(m, Y, t); % fi: p x m, dfdy: m x m x p


        % ------- 2) 计算残差 R = Y - y0 - dt * S * F -------
        R = Y - repmat(y0, p, 1) - dt * (S * f); % p x m
        res_norm = max(abs(R(:)));
        if res_norm < tol
            converged = true;
            break;
        end

        % ------- 3) 组装大雅可比矩阵 Jbig （pm x pm） -------
        pm = p * m;
        Jbig = zeros(pm, pm);  % 若 p*m 很大，建议用 sparse

        for i = 1:p
            row_idx = (i-1)*m + (1:m);    % 大矩阵中第 i 块行的索引
            for j = 1:p
                col_idx = (j-1)*m + (1:m); % 第 j 块列的索引

                % 注意这里要用 dfdy 的转置：见后面解释
                blk = - dt * S(i, j) * (dfdy(:,:,j)'); % m x m

                if i == j
                    blk = blk + eye(m);
                end
                Jbig(row_idx, col_idx) = blk;
            end
        end

        % ------- 4) 求解线性增量方程 Jbig * dYvec = -Rvec -------
        Rvec = reshape(R', pm, 1);       % 按行堆叠成列向量（vec of R'）
        dYvec = - (Jbig \ Rvec);         % 求解（直接解法，视规模可换 GMRES）
        dY = reshape(dYvec, m, p)';      % 恢复为 p x m

        % ------- 5) 更新与收敛判断 -------
        Y = Y + dY;

        if max(abs(dY(:))) < tol
            converged = true;
            break;
        end

        iter = iter + 1;
    end

    if ~converged
        fprintf('Warning: Newton did not converge in %d iters, residual=%g\n', max_iter, res_norm);
    end

    y_solution = Y;  % p x m
end

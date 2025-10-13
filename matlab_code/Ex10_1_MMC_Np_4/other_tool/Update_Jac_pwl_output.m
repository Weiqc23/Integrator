function [] = Update_Jac_pwl_output(p,m,dt,S)

global A Bsw PuPx_sw
dfdy = A + Bsw*PuPx_sw';

global Jac_pwl_output

if_Jac_expand = ~isempty(S);
if (if_Jac_expand)
    pm = p * m;
    Jac_pwl_output = zeros(pm, pm);  % 若 p*m 很大，建议用 sparse

    for i = 1:p
        row_idx = (i-1)*m + (1:m);    % 大矩阵中第 i 块行的索引
        for j = 1:p
            col_idx = (j-1)*m + (1:m); % 第 j 块列的索引
            % 注意这里要用 dfdy 的转置：见后面解释
            blk = - dt * S(i, j) * dfdy; % m x m
            Jac_pwl_output(row_idx, col_idx) = blk;
        end
    end

else
    Jac_pwl_output = dfdy';
end

end
function ynew=predictor(y0,tnode,m,n, kmax, etol)

ynew=ones(n+1,1)*y0;

global ExImplicit
if (ExImplicit==0)
    return;
else
    y = y0;
    ynew(1,:) = y;
    for i = 1:n
    t0 = tnode(i);
    t = tnode(i+1);
    % y = newton_raphson_TR_onestep(y, t0, t, m, kmax, etol);
    y = newton_raphson_BE_onestep(y, t0, t, m, kmax, etol);
    % y = newton_raphson_FE_onestep(y, t0, t, m);
    ynew(i+1,:) = y;
    end
end

end

function y = newton_raphson_BE_onestep(y0, t0, t, m, kmax, etol)
% dy/dt = f(y,t) 用TR
% F(y) = h*f(y,t) - (y-y0) = 0;

    y = y0;
    h = t-t0;

    iter = 0;         % 迭代计数器

    while iter < kmax
        % 计算残差 (方程左边 - iLight)
        f = fright(m, y, t);
        residual = h*f-(y-y0);

        % 检查收敛条件
        if max(abs(residual)) < etol
            break;
        end

        % 计算雅可比矩阵 (导数)
        dfdy = jeval(m,y,t);
        J = h*dfdy-eye(m);

        % 牛顿迭代更新: Vd_new = Vd - residual / J
        y = y- residual/J;
        % y = y -  J\residual;

        % % 打印当前迭代信息
        % fprintf('%d\t%.6f\t%.6e\n', iter, Vd, residual);

        iter = iter + 1;
    end

    if iter >= kmax
        fprintf('Warning! Reach the maximum number of iterations (%d)\n', int32(kmax));
    end

end
function y = newton_raphson_BE_onestep(y0, t0, t, m)
% dy/dt = f(y,t) 用TR
% F(y) = h*f(y,t) - (y-y0) = 0;

global set_kmax set_etol
    be_tol = set_etol;
    k_max = set_kmax;
    y = y0;
    h = t-t0;

    iter = 0;         % 迭代计数器

    while iter < k_max
        % 计算残差 (方程左边 - iLight)
        f = fright(m, y, t);
        residual = h*f-(y-y0);

        % 检查收敛条件
        if max(abs(residual)) < be_tol
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
        % disp("newton_raphson_BE_onestep:");
        % disp("f");
        % disp(num2str(f, '%.12f '));
        % disp("y0");
        % disp(num2str(y0, '%.12f '));
        % disp("y");
        % disp(num2str(y, '%.12f '));
        iter = iter + 1;
    end

    if iter >= k_max
        fprintf('Warning! Reach the maximum number of iterations (%d)\n', int32(k_max));
    end

end
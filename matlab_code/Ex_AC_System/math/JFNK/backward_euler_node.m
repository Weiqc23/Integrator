% function [y] = backward_euler_node(f_eval, t, h, rhs, y0, be_tol) 
% 
%     y = zeros(size(rhs));
%     if (h==0)
%         y = rhs;
%     else
%         % 调用求解器计算y
%         F = @(x) x - h * f_eval(t, x) - rhs;
%         options = optimoptions('fsolve', 'Display', 'off', 'TolFun', be_tol);
%         y = fsolve(F, y0, options);
% 
%         % options = optimoptions('fsolve', ...
%                                'TolFun', be_tol, ...
%                                'Algorithm', 'trust-region-dogleg', ...
%                                'Display', 'off');
%         % tic;
%         % [y,~,exitflag] = fsolve(F, y0, options);
%         % elapsed = toc;
%         % fprintf('一次fsolve耗时 %.6f 秒\n', elapsed);
% 
%         % if exitflag <= 0
%         %     error('backward_euler_node:fsolve_fail', ...
%         %           'fsolve did not converge at t=%g, h=%g', t, h);
%         % end
%     end
% 
% end

% 牛顿迭代求解y - h*f(t(i),y) = rhs
function y = backward_euler_node(f_eval, t, h, rhs, y0, be_tol, do_print)
    global A B C D m mi mo 
    global set_kmax set_etol
    be_tol = set_etol;
    k_max = set_kmax;
    y = y0(:)';

    iter = 0;  

    while iter < k_max
    % while (1)
        % 计算残差 (方程左边 - iLight)
        f = f_eval(t,y)';
        residual = h*f-(y-rhs');

        % 检查收敛条件
        if max(abs(residual)) < be_tol
            break;
        end

        % 计算雅可比矩阵 (导数)
        dfdy = jeval(m,y,t);
        J = h*dfdy-eye(m);

        % 牛顿迭代更新: Vd_new = Vd - residual / J
        delta_y = residual/J;
        % disp("J");
        % disp(num2str(J,'%.12f '));
        % disp("residual");
        % disp(num2str(residual,'%.12f '));
        % disp("delta_y");
        % disp(num2str(delta_y, '%.12f '));
        % disp("y");
        % disp(num2str(y, '%.12f '));
        y = y- residual/J;
        % y = y -  J\residual;

        % % 打印当前迭代信息
        % fprintf('%d\t%.6f\t%.6e\n', iter, Vd, residual);
        % disp("back_euler_node:");
        % disp("f");
        % disp(num2str(f, '%.12f '));   
        % disp("y0");
        % disp(num2str(y0', '%.12f '));
        % disp("y");
        % disp(num2str(y, '%.12f '));
        iter = iter + 1;
    end

    y = y(:);
    if iter >= k_max
        fprintf('Warning! Reach the maximum number of iterations (%d)\n', int32(k_max));
    end
end


% 牛顿迭代求解y - h*f(t(i),y) = rhs(但只迭代一次)
% function y = backward_euler_node(f_eval, t, h, rhs, y0, be_tol, do_print)
%     global A B C D m mi mo 
%     global set_kmax set_etol
%     be_tol = set_etol;
%     k_max = set_kmax;
%     y = y0(:)';
% 
%     iter = 0;  
% 
%     while iter < 1
%         % 计算残差 (方程左边 - iLight)
%         f = f_eval(t,y)';
%         residual = h*f-(y-rhs');
% 
%         % 检查收敛条件
%         if max(abs(residual)) < be_tol
%             break;
%         end
% 
%         % 计算雅可比矩阵 (导数)
%         dfdy = jeval(m,y,t);
%         J = h*dfdy-eye(m);
% 
%         % 牛顿迭代更新: Vd_new = Vd - residual / J
%         delta_y = residual/J;
%         % disp("J");
%         % disp(num2str(J,'%.12f '));
%         % disp("residual");
%         % disp(num2str(residual,'%.12f '));
%         % disp("delta_y");
%         % disp(num2str(delta_y, '%.12f '));
%         % disp("y");
%         % disp(num2str(y, '%.12f '));
%         y = y- residual/J;
%         % y = y -  J\residual;
% 
%         % % 打印当前迭代信息
%         % fprintf('%d\t%.6f\t%.6e\n', iter, Vd, residual);
%         % disp("back_euler_node:");
%         % disp("f");
%         % disp(num2str(f, '%.12f '));   
%         % disp("y0");
%         % disp(num2str(y0', '%.12f '));
%         % disp("y");
%         % disp(num2str(y, '%.12f '));
%         iter = iter + 1;
%     end
% 
%     y = y(:);
%     if iter >= k_max
%         fprintf('Warning! Reach the maximum number of iterations (%d)\n', int32(k_max));
%     end
% end



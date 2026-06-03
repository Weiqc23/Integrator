% 线性版本sdc求解y - h*f(t(i),y) = rhs
% y - h*(Ay + Bu) = rhs
function y = backward_euler_node_sdc(f_eval, t, h, rhs, y0, be_tol, do_print)
    global A B C D m mi mo 
    % disp(['size(A)=', mat2str(size(A)), ', size(B)=', mat2str(size(B)), ', size(y)=', mat2str(size(y))]);
    y0 = y0(:);

    m = numel(rhs);
    I = speye(m);
    ut = Get_u_From_x(y0,t);

    % 线性系统 (I - h*A) * y = rhs + h*B*u(t)
    y = (I - h*A) \ (rhs + h*B*ut);
end



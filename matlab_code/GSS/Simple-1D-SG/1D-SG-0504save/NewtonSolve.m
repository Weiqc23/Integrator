function [] = NewtonSolve()
    global n_node N_base
    global J_mat delta_x fx
    global psi n p N
    % 使用牛顿法计算

    fx = Ftotal();
    fprintf('Initial norm = %.6f\n', norm(fx, 2));
    J_mat = zeros(3 * (n_node - 2), 3 * (n_node - 2));


    iter = 1;
    step_size = 0.1;
    err_tol = 1e-6;
    
    while 1
        % 计算牛顿法的修正
        
        UpdateJ();

        J_mat = sparse(J_mat);

        delta_x = -J_mat \ fx;

        % 高精度求解
        % delta_x = -mp(J_mat) \ mp(fx);


        % 通过 damping 计算步长
        step_size = BankRoseDamping();
        if step_size < 1e-10
            error('Bank-Rose damping fail, cannot converge\n');
        end


        % 修正现在的解
        psi(2:end-1) = psi(2:end-1) + step_size * delta_x(1:n_node-2);
        n(2:end-1) = n(2:end-1) + step_size * delta_x((n_node-2)+1:2*(n_node-2));
        p(2:end-1) = p(2:end-1) + step_size * delta_x(2*(n_node-2)+1:3*(n_node-2));
    
        % 修一下负浓度
        for idx = 1:n_node
            if n(idx) < 0
                n(idx) = 1e-100;
            end
            if p(idx) < 0
                p(idx) = 1e-100;
            end
        end


        fx = Ftotal();

        CalculateVariables();
        RealTimePlot(); % 实时绘制当前求解的状态
        subplot(4, 4, 16)
        plot(delta_x)
        title("dx")
    
        % condition number
        condition = condest(J_mat);

        fprintf('iter %d, condition(J) = %.3e, step = %.2e, step*|dx| = %.5e, |fx| = %.5e\n', iter, condition, step_size, step_size * norm(delta_x, 2), norm(fx, 2));


        iter = iter + 1;
        if norm(fx, 2) < err_tol
            break;
        end
    end
    fprintf('Newton solver converge\n');
end


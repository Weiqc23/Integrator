% function[y,Y,D,is_converged,is_stiff,ratios] = JFNK(f_eval,t0,t,t_m,y0,y_approx,y_approx_m,S,S_m,S_p,S_p_m,spectral_radius,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)
function[y,Y,D,is_converged,is_stiff,ratios] = JFNK(f_eval, t, y0, y_approx, S, S_p, spectral_radius,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)

    % [y_init,Y_init,D_init,is_converged,is_stiff,ratios] = JFNK_initial(@f_eval,t0,t,t_m,y0,y_approx,y_approx_m,S,S_m,S_p,S_p_m,spectral_radius, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
    [y_init,Y_init,D_init,is_converged,is_stiff,ratios] = JFNK_initial(@f_eval,t,y0,y_approx,S,S_p,spectral_radius, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
    % disp('y_init');
    % disp(num2str(y_init, '%.12f '));
    % fprintf("\n\n");
    y = y_init;
    Y = Y_init;
    D = D_init;

    % Y = {};
    % D = {};
    
    % disp("y_init");
    % disp(y_init);

    if (is_converged)
        y = y_init;
        return;
    end
    
    if (~is_converged && is_stiff)
        [y_new,Y_new,D_new,is_converged] = JFNK_iterations(@f_eval,t,y0,y_init,S,S_p, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
        Y = [Y_init,Y_new];
        D = [D_init,D_new];
        y = y_new;
        % disp('y');
        % disp(num2str(y, '%.12f '));
        % fprintf("\n\n");
    end
    

end
function[y,Y,D,is_converged,is_stiff,ratios] = JFNK_step(f_eval, p, dt, t, y0,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)
    
    S = dt*p.S;
    S_p = dt*p.S_p;

    [y_be] = backward_euler(@f_eval,t,y0,S_p,be_tol);
    % disp("predictor:");
    % disp(num2str(y_be, '%.12f '));
    [y,Y,D,is_converged,is_stiff,ratios] = JFNK(@f_eval, t, y0, y_be, S, S_p, p.spectral_radius ,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
    
end
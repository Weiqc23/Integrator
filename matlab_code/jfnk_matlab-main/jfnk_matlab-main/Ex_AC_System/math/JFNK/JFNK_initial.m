% function[y_sdc_m,Y_m,D_m,is_converged,is_stiff,ratios_m] = JFNK_initial(f_eval, t0, t, t_m, y0, y_approx, y_approx_m,S,S_m,S_p,S_p_m, spectral_radius, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)
    function[y_sdc,Y,D,is_converged,is_stiff,ratios] = JFNK_initial(f_eval, t, y0, y_approx, S, S_p, spectral_radius, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)

    is_stiff = false;
    
    i = 0;
    % ratios = {};
    ratios = {};
    % disp(['JFNK_initial:size(y_approx)=', mat2str(size(y_approx)), ', numel(y_approx)=', num2str(numel(y_approx))]);
    % [y_sdc,y_sdc_m,Y_sdc,D_sdc,D_sdc_m,is_converged] = SDC(@f_eval, t0, t, t_m, y0, y_approx, y_approx_m, S, S_m, S_p, S_p_m, 2, be_tol, sdc_tol);
    % [y_sdc,Y_sdc,D_sdc,is_converged] = SDC(@f_eval, t0, t, y0, y_approx, S, S_p, 2, be_tol, sdc_tol);
    [y_sdc,Y_sdc,D_sdc,is_converged] = SDC(@f_eval, t, y0, y_approx, S, S_p, 2, be_tol, sdc_tol);

    % Y=Y_sdc;
    Y=Y_sdc;
    % D=D_sdc;
    D=D_sdc;

    if (is_converged)
        % ratios = {};
        ratios = {};
    end
    
    while ((~is_converged) && (~is_stiff) && (i<n_iter_max_sdc))
        
        ratio = norm(D{end}, 'fro') / norm(D{end-1}, 'fro');
        
        ratios{end+1} = ratio;
        
        % is_stiff = ((ratio/spectral_radius)>0.1);
        % 强制非线性
        is_stiff = false;

        if (~is_stiff)
            % disp(['JFNK_initial:size(y_sdc)=', mat2str(size(y_sdc)), ', numel(y_old)=', num2str(numel(y_sdc))]);
            [y_sdc,Y_sdc,D_sdc,is_converged] = SDC(@f_eval,t,y0,y_sdc,S,S_p,1,be_tol,sdc_tol);
            % [y_sdc_m,Y_sdc_m,D_sdc_m,is_converged] = SDC(@f_eval,t0,t_m,y0,y_sdc_m,S_m,S_p_m,1,be_tol,sdc_tol);
            Y = {Y{:}, Y_sdc{1}};
            % Y_m = {Y_m{:}, Y_sdc_m{1}};
            D = {D{:}, D_sdc{1}};
            % D_m = {D_m{:}, D_sdc_m{1}};
        end
        i=i+1;
    end

    if (i>=n_iter_max_sdc) 
        % fprintf('initial_error\n');
    end

    % disp(['jfnk_initial']);
    % disp('y0');
    % disp(num2str(y0', '%.12f '));
    % disp('y_approx');
    % disp(num2str(y_approx, '%.12f '));
    % disp('y_sdc');
    % disp(num2str(y_sdc, '%.12f '));
    % fprintf("\n\n");

end


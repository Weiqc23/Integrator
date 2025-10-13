

function [t, y, Y, D] = JFNK_uniform(f_eval, t_init, t_final, p, y0, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)

    dt = t_final-t_init;
    
    if (p.t(end) ~= 1) 
        t_step = [p.t,1];
    else
        t_step = p.t;
    end

    nodes = dt*t_step + t_init;
    nodes = nodes(:);

    % t_init~t_final中真实的节点
    is_stiff = false;
    is_converged = false;
    ratios = {};
    y0 = y0(:);
    % y_all = y0';
    % t_all = t_init;
    % Y_all = {};
    % D_all = {};
    
    % [y, Y, D, is_converged, is_stiff, ratios] = JFNK_step(@f_eval, p,dt,t_init,nodes,nodes_m,y0,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
    [y, Y, D, is_converged, is_stiff, ratios] = JFNK_step(@f_eval, p,dt,nodes,y0,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
    t = nodes;
    % y_all = [y_all;y(1:end,:)];
    % t_all = [t_all;nodes(2:end)];
    % Y_all{end+1} = Y;
    % D_all{end+1} = D;

end
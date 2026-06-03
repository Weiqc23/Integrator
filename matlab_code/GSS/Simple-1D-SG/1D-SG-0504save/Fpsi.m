function [f_psi_] = Fpsi()
    % 计算泊松方程的残差并返回
    global n_node dx
    global vt nie N_base V_base
    global psi phi_n phi_p N n p
    global lambda2

    f_psi_ = zeros(n_node - 2, 1);

    for idx = 2:(n_node-1)
        % 仅遍历非边界格点

        psi_w    = psi(idx-1);
        psi_o    = psi(idx);
        psi_e    = psi(idx+1);

        n_o  = n(idx);
        p_o  = p(idx);
        

        f_psi_(idx-1) = lambda2 * (psi_w + psi_e - 2 * psi_o) / dx ^ 2 ...
        + p_o - n_o + N(idx);
        
    end
    
end


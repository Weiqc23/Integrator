function [f_n_] = Fn()
    % 计算电子连续性方程的残差并返回
    global n_node dx
    global vt nie
    global psi phi_n phi_p mu_n tau_n tau_p n p

    f_n_ = zeros(n_node - 2, 1);

    for idx = 2:(n_node-1)
        % 仅遍历非边界格点

        psi_w    = psi(idx-1);
        psi_o    = psi(idx);
        psi_e    = psi(idx+1);

        n_w  = n(idx-1);
        n_o  = n(idx);
        n_e  = n(idx+1);

        p_o  = p(idx);

        R = (n_o * p_o - nie ^ 2) / (tau_n * (n_o + nie) + (tau_p * (p_o + nie)));
        % R = 0;

        % 标准 SG 差分格式(使用 bern)
        delta_w = (psi_o - psi_w) / vt;
        delta_e = (psi_e - psi_o) / vt;

        % 默认方向维 dpsi 右减左，m_dpsi 为左减右
        bern_dpsi_e = Bern(delta_e);
        bern_m_dpsi_e = Bern(-delta_e);
        bern_dpsi_w = Bern(delta_w);
        bern_m_dpsi_w = Bern(-delta_w);

        J_n_e = mu_n * vt / dx * (bern_dpsi_e * n_e - bern_m_dpsi_e * n_o);
        J_n_w = mu_n * vt / dx * (bern_dpsi_w * n_o - bern_m_dpsi_w * n_w);

        f_n_(idx-1) = (J_n_e - J_n_w) / dx - R;
    end
end


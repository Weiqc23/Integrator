function [] = CalculateVariables()
    global psi phi_n phi_p n p J_n J_p J R E
    global  mu_n mu_p vt tau_n tau_p nie
    global V_base N_base J_base L_base G_base
    global n_node dx

    % 基于当前计算得到的 psi n p，还原 n p J 等物理量的值

    % 均使用无量纲化后的量计算，得到的也是无量纲化后的值

    % n = nie * exp((psi - phi_n) / vt);
    % p = nie * exp((phi_p - psi) / vt);

    phi_n = psi - vt * log(n / nie);
    phi_p = psi + vt * log(p / nie);
     
    % R = (n .* p - nie ^ 2) ./ (tau_n * (n + nie) + tau_p * (p + nie));
    R = 0;

    E = (psi(1:end-1) - psi(2:end)) / dx;

    %% 简单差分
    for idx = 1:n_node-1
        % 使用无量纲的量计算

        n_l = n(idx);
        n_r = n(idx+1);
        p_l = p(idx);
        p_r = p(idx+1);

        dpsi = (psi(idx+1) - psi(idx)) / vt;
        bern_dpsi = Bern(dpsi);
        bern_m_dpsi = Bern(-dpsi);

        J_n(idx) = mu_n * vt / dx * (bern_dpsi * n_r - bern_m_dpsi * n_l);
        J_p(idx) = mu_p * vt / dx * (bern_dpsi * p_l - bern_m_dpsi * p_r);

    end

    J = J_n + J_p;
end


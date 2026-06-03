function [] = UpdateJ()
    global nie vt mu_n mu_p tau_n tau_p
    global n_node dx
    global psi n p
    global J_mat
    global lambda2

    % 计算并组装雅可比矩阵
    % 使用 SG 格式
    % 伯努利函数及其导数使用 Bern.m 和 Pd1Bern.m 计算



    for position = 2:(n_node-1)
        % idx 为矩阵下标
        idx = position - 1;

        %% 计算需要的变量
        psi_w    = psi(position-1);
        psi_o    = psi(position);
        psi_e    = psi(position+1);

        n_w  = n(position-1);
        n_o  = n(position);
        n_e  = n(position+1);
        
        p_w  = p(position-1);
        p_o  = p(position);
        p_e  = p(position+1);

        % 标准 SG 差分格式(使用 bern)
        delta_w = (psi_o - psi_w) / vt;
        delta_e = (psi_e - psi_o) / vt;

        % 伯努利函数
        % 默认方向维 dpsi 右减左，m_dpsi 为左减右
        bern_dpsi_e = Bern(delta_e);
        bern_m_dpsi_e = Bern(-delta_e);
        bern_dpsi_w = Bern(delta_w);
        bern_m_dpsi_w = Bern(-delta_w);
        
        % 伯努利函数导数
        % 默认方向维 dpsi 右减左，m_dpsi 为左减右
        d_berne_dpsie = Pd1Bern(delta_e) / vt;
        d_berne_dpsio = -d_berne_dpsie;
        d_bernme_dpsie = -Pd1Bern(-delta_e) / vt;
        d_bernme_dpsio = -d_bernme_dpsie;

        d_bernw_dpsio = Pd1Bern(delta_w) / vt;
        d_bernw_dpsiw = -d_bernw_dpsio;
        d_bernmw_dpsio = -Pd1Bern(-delta_w) / vt;
        d_bernmw_dpsiw = -d_bernmw_dpsio;

        
        R = (n_o * p_o - nie ^ 2) / (tau_n * (n_o + nie) + (tau_p * (p_o + nie)));
        dRdn = -(tau_n * (n_o * p_o - nie^2)) / (tau_n * (n_o + nie) + tau_p * (p_o + nie))^2 + p_o / (tau_n * (n_o + nie) + tau_p * (p_o + nie)); % R对n的导数
        dRdp = -(tau_p * (p_o * n_o - nie^2)) / (tau_p * (p_o + nie) + tau_n * (n_o + nie))^2 + n_o / (tau_p * (p_o + nie) + tau_n * (n_o + nie)); % R对p的导数
        % R = 0;
        % dRdn = 0;
        % dRdp = 0;

        %% 泊松方程
        if idx == 1
            % 若是 J 的第 1 行，即不需要组装 w 格点的位置
            % d f_psi / d psi_o
            J_mat(idx, idx) = -lambda2 * 2 / dx^2;
            % d f_psi / d psi_e
            J_mat(idx, idx+1) = lambda2 / dx^2;
            % d f_psi / d n_o
            J_mat(idx, idx+(n_node-2)) = -1.0;
            % d f_psi / d p_o
            J_mat(idx, idx+(n_node-2)*2) = 1.0;
        elseif idx == (n_node - 2)
            % 若是 J 的第 n_node - 2 行，即不需要组装 e 格点的位置
            % d f_psi / d psi_w
            J_mat(idx, idx-1) = lambda2 / dx^2;
            % d f_psi / d psi_o
            J_mat(idx, idx) = -lambda2 * 2 / dx^2;
            % d f_psi / d n_o
            J_mat(idx, idx+(n_node-2)) = -1.0;
            % d f_psi / d p_o
            J_mat(idx, idx+(n_node-2)*2) = 1.0;
        else
            % d f_psi / d psi_w
            J_mat(idx, idx-1) = lambda2 / dx^2;
            % d f_psi / d psi_o
            J_mat(idx, idx) = -lambda2 * 2 / dx^2;
            % d f_psi / d psi_e
            J_mat(idx, idx+1) = lambda2 / dx^2;
            % d f_psi / d n_o
            J_mat(idx, idx+(n_node-2)) = -1.0;
            % d f_psi / d p_o
            J_mat(idx, idx+(n_node-2)*2) = 1.0;
        end
        %% 电子连续性方程
        if idx == 1
            % 若是第 1 行，即不需要组装 w 格点的位置
            % d f_n / d psi_o
            J_mat(idx+(n_node-2), idx) = mu_n * vt / dx^2 * ( ...
                (d_berne_dpsio * n_e - d_bernme_dpsio * n_o) ...
                -(d_bernw_dpsio * n_o - d_bernmw_dpsio * n_w) ...
            );
            % d f_n / d psi_e
            J_mat(idx+(n_node-2), idx+1) = mu_n * vt / dx^2 * (d_berne_dpsie * n_e - d_bernme_dpsie * n_o);
            % d f_n / d n_o
            J_mat(idx+(n_node-2), idx+(n_node-2)) = mu_n * vt / dx^2 * (-bern_m_dpsi_e - bern_dpsi_w) - dRdn;
            % d f_n / d n_e
            J_mat(idx+(n_node-2), idx+(n_node-2)+1) = mu_n * vt / dx^2 * (bern_dpsi_e);
            % d f_n / d p_o
            J_mat(idx+(n_node-2), idx+2*(n_node-2)) = -dRdp;
        elseif idx == (n_node - 2)
            % 若是最后一行，即不需要组装 e 格点的位置
            % d f_n / d psi_w
            J_mat(idx+(n_node-2), idx-1) = -mu_n * vt / dx^2 * (d_bernw_dpsiw * n_o - d_bernmw_dpsiw * n_w);
            % d f_n / d psi_o
            J_mat(idx+(n_node-2), idx) = mu_n * vt / dx^2 * ( ...
                (d_berne_dpsio * n_e - d_bernme_dpsio * n_o) ...
                -(d_bernw_dpsio * n_o - d_bernmw_dpsio * n_w) ...
            );
            % d f_n / d n_w
            J_mat(idx+(n_node-2), idx+(n_node-2)-1) = -mu_n * vt / dx^2 * (-bern_m_dpsi_w);
            % d f_n / d n_o
            J_mat(idx+(n_node-2), idx+(n_node-2)) = mu_n * vt / dx^2 * (-bern_m_dpsi_e - bern_dpsi_w) - dRdn;
            % d f_n / d p_o
            J_mat(idx+(n_node-2), idx+2*(n_node-2)) = -dRdp;
        else
            % d f_n / d psi_w
            J_mat(idx+(n_node-2), idx-1) = -mu_n * vt / dx^2 * (d_bernw_dpsiw * n_o - d_bernmw_dpsiw * n_w);
            % d f_n / d psi_o
            J_mat(idx+(n_node-2), idx) = mu_n * vt / dx^2 * ( ...
                (d_berne_dpsio * n_e - d_bernme_dpsio * n_o) ...
                -(d_bernw_dpsio * n_o - d_bernmw_dpsio * n_w) ...
            );
            % d f_n / d psi_e
            J_mat(idx+(n_node-2), idx+1) = mu_n * vt / dx^2 * (d_berne_dpsie * n_e - d_bernme_dpsie * n_o);
            % d f_n / d n_w
            J_mat(idx+(n_node-2), idx+(n_node-2)-1) = -mu_n * vt / dx^2 * (-bern_m_dpsi_w);
            % d f_n / d n_o
            J_mat(idx+(n_node-2), idx+(n_node-2)) = mu_n * vt / dx^2 * (-bern_m_dpsi_e - bern_dpsi_w) - dRdn;
            % d f_n / d n_e
            J_mat(idx+(n_node-2), idx+(n_node-2)+1) = mu_n * vt / dx^2 * (bern_dpsi_e);
            % d f_n / d p_o
            J_mat(idx+(n_node-2), idx+2*(n_node-2)) = -dRdp;
        end
        %% 空穴连续性方程
        if idx == 1
            % 若是第 1 行，即不需要组装 w 格点的位置
            % d f_p / d psi_o
            J_mat(idx+2*(n_node-2), idx) = -mu_p * vt / dx^2 * ( ...
                (d_berne_dpsio * p_o - d_bernme_dpsio * p_e) ...
                -(d_bernw_dpsio * p_w - d_bernmw_dpsio * p_o) ...
            );
            % d f_p / d psi_e
            J_mat(idx+2*(n_node-2), idx+1) = -mu_p * vt / dx^2 * (d_berne_dpsie * p_o - d_bernme_dpsie * p_e);
            % d f_p / d n_o
            J_mat(idx+2*(n_node-2), idx+(n_node-2)) = -dRdn;
            % d f_p / d p_o
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)) = -mu_p * vt / dx^2 * ((bern_dpsi_e) - (-bern_m_dpsi_w)) - dRdp;
            % d f_p / d p_e
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)+1) = -mu_p * vt / dx^2 * (-bern_m_dpsi_e);
        elseif idx == (n_node - 2)
            % 若是最后一行，即不需要组装 e 格点的位置
            % d f_p / d psi_w
            J_mat(idx+2*(n_node-2), idx-1) = mu_p * vt / dx^2 * (d_bernw_dpsiw * p_w - d_bernmw_dpsiw * p_o);
            % d f_p / d psi_o
            J_mat(idx+2*(n_node-2), idx) = -mu_p * vt / dx^2 * ( ...
                (d_berne_dpsio * p_o - d_bernme_dpsio * p_e) ...
                -(d_bernw_dpsio * p_w - d_bernmw_dpsio * p_o) ...
            );
            % d f_p / d n_o
            J_mat(idx+2*(n_node-2), idx+(n_node-2)) = -dRdn;
            % d f_p / d p_w
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)-1) = mu_p * vt / dx^2 * (bern_dpsi_w);
            % d f_p / d p_o
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)) = -mu_p * vt / dx^2 * ((bern_dpsi_e) - (-bern_m_dpsi_w)) - dRdp;
        else
            % d f_p / d psi_w
            J_mat(idx+2*(n_node-2), idx-1) = mu_p * vt / dx^2 * (d_bernw_dpsiw * p_w - d_bernmw_dpsiw * p_o);
            % d f_p / d psi_o
            J_mat(idx+2*(n_node-2), idx) = -mu_p * vt / dx^2 * ( ...
                (d_berne_dpsio * p_o - d_bernme_dpsio * p_e) ...
                -(d_bernw_dpsio * p_w - d_bernmw_dpsio * p_o) ...
            );
            % d f_p / d psi_e
            J_mat(idx+2*(n_node-2), idx+1) = -mu_p * vt / dx^2 * (d_berne_dpsie * p_o - d_bernme_dpsie * p_e);
            % d f_p / d n_o
            J_mat(idx+2*(n_node-2), idx+(n_node-2)) = -dRdn;
            % d f_p / d p_w
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)-1) = mu_p * vt / dx^2 * (bern_dpsi_w);
            % d f_p / d p_o
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)) = -mu_p * vt / dx^2 * ((bern_dpsi_e) - (-bern_m_dpsi_w)) - dRdp;
            % d f_p / d p_e
            J_mat(idx+2*(n_node-2), idx+2*(n_node-2)+1) = -mu_p * vt / dx^2 * (-bern_m_dpsi_e);
        end
    end
end


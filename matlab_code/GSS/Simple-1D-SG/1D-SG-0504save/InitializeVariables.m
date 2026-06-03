function [] = InitializeVariables()
    global L_base N_base t_base mu_base V_base TAU_N TAU_P MU_N MU_P NIE VT
    global P_doping N_doping device_length
    global n_node dx
    global n p N psi phi_n phi_p mu_n mu_p tau_n tau_p E J_n J_p J R nie vt

    
    %% 主要变量声明

    % 定义在格点上的变量
    n = zeros(n_node, 1);
    p = zeros(n_node, 1);
    N = zeros(n_node, 1);
    psi = zeros(n_node, 1);
    phi_n = zeros(n_node, 1);
    phi_p = zeros(n_node, 1);
    R = zeros(n_node, 1);

    % 定义在格点中间的变量
    E = zeros(n_node - 1, 1);
    J_n = zeros(n_node - 1, 1);
    J_p = zeros(n_node - 1, 1);
    J = zeros(n_node - 1, 1);
    
    % 直接赋归一化后的值
    dx = device_length / (n_node - 1) * 1e-4 / L_base; % 归一化后的格点间距。length 单位 um，先换成厘米再归一化
    tau_n = TAU_N / t_base;
    tau_p = TAU_P / t_base;
    mu_n = MU_N / mu_base;
    mu_p = MU_P / mu_base;
    nie = NIE / N_base;
    vt = VT / V_base;

    %% Initial guess: Constant Abrupt PN junction

    for index = 1:floor((n_node - 1) / 2)
        % P 区
        % 直接赋归一化后的值
        N(index) = -P_doping / N_base; % P 区净掺杂
        n(index) = (NIE ^ 2 / P_doping) / N_base;
        p(index) = P_doping / N_base;
        phi_n(index) = 0;
        phi_p(index) = 0;
        psi(index) = (-vt * log(P_doping / NIE));
    end
    for index = floor((n_node - 1) / 2)+1:n_node
        % N 区
        % 直接赋归一化后的值
        N(index) = N_doping / N_base; % N 区净掺杂
        n(index) = N_doping / N_base;
        p(index) = (NIE ^ 2 / N_doping) / N_base;
        phi_n(index) = 0;
        phi_p(index) = 0;
        psi(index) = (vt * log(N_doping / NIE));
    end


    %% Initial guess: Pure P Doping

%     for index = 1:length(n)
%         % P 区
%         % 直接赋归一化后的值
%         N(index) = -P_doping / N_base; % P 区净掺杂
%         n(index) = (NIE ^ 2 / P_doping) / N_base;
%         p(index) = P_doping / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (-vt * log(P_doping / NIE));
%     end


    %% Initial guess: PNPN

%     for index = 1:floor((n_node - 1) / 4)
%         % P 区
%         % 直接赋归一化后的值
%         N(index) = -P_doping / N_base; % P 区净掺杂
%         n(index) = (NIE ^ 2 / P_doping) / N_base;
%         p(index) = P_doping / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (-vt * log(P_doping / NIE));
%     end
%     for index = floor((n_node - 1) / 4)+1:floor((n_node - 1) * 2 / 4)
%         % N 区
%         % 直接赋归一化后的值
%         N(index) = N_doping / N_base; % N 区净掺杂
%         n(index) = N_doping / N_base;
%         p(index) = (NIE ^ 2 / N_doping) / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (vt * log(N_doping / NIE));
%     end
%     for index = floor((n_node - 1) * 2 / 4)+1:floor((n_node - 1) * 3 / 4)
%         % P 区
%         % 直接赋归一化后的值
%         N(index) = -P_doping / N_base; % P 区净掺杂
%         n(index) = (NIE ^ 2 / P_doping) / N_base;
%         p(index) = P_doping / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (-vt * log(P_doping / NIE));
%     end
%     for index = floor((n_node - 1) * 3 / 4) + 1:n_node
%         % N 区
%         % 直接赋归一化后的值
%         N(index) = N_doping / N_base; % N 区净掺杂
%         n(index) = N_doping / N_base;
%         p(index) = (NIE ^ 2 / N_doping) / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (vt * log(N_doping / NIE));
%     end





    %% Initial guess: PNP

%     for index = 1:floor((n_node - 1) / 3)
%         % P 区
%         % 直接赋归一化后的值
%         N(index) = -P_doping / N_base; % P 区净掺杂
%         n(index) = (NIE ^ 2 / P_doping) / N_base;
%         p(index) = P_doping / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (-vt * log(P_doping / NIE));
%     end
%     for index = floor((n_node - 1) / 3)+1:floor((n_node - 1) * 2 / 3)
%         % N 区
%         % 直接赋归一化后的值
%         
%         relax = 1e-3;
% 
%         N(index) = N_doping / N_base * relax; % N 区净掺杂
%         n(index) = N_doping / N_base * relax;
%         p(index) = (NIE ^ 2 / (N_doping * relax)) / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (vt * log(N_doping * relax / NIE));
%     end
%     for index = floor((n_node - 1) * 2 / 3)+1:n_node
%         % P 区
%         % 直接赋归一化后的值
%         N(index) = -P_doping / N_base; % P 区净掺杂
%         n(index) = (NIE ^ 2 / P_doping) / N_base;
%         p(index) = P_doping / N_base;
%         phi_n(index) = 0;
%         phi_p(index) = 0;
%         psi(index) = (-vt * log(P_doping / NIE));
%     end


end

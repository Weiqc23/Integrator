function [] = InitializeConstants()
    global L_base V_base E_base N_base J_base G_base t_base mu_base
    global q VT NIE
    global TAU_N TAU_P MU_N MU_P
    global lambda2 delta2

    % 物理常数
    EPSILON_VACUUM = 8.8541878128e-14;
    EPSILON_SI = EPSILON_VACUUM * 11.7;
    q = 1.602176565e-19;
    N_C = 2.8e19;
    N_V = 1.04e19;
    E_G = 1.08;
    BOLTZMANN_CONSTANT = 1.38064852e-23;
    VT = BOLTZMANN_CONSTANT * 300 / q;  % kT/q
    NIE = sqrt(N_C * N_V) * exp(-E_G / (2 * VT));

    % 使用常值模型
    TAU_N = 1e-7;
    TAU_P = 1e-7;
    MU_N = 1000;
    MU_P = 500;


    % 归一化常数
    % 考虑四个独立的常数
    V_base = VT;
    N_base = 1e17;
    D_base = 27;
    L_base = 30 * 1e-4; % 20 um

    t_base = L_base ^ 2 / D_base;
    E_base = V_base / L_base;
    J_base = q * D_base * N_base / L_base;
    mu_base = D_base / V_base;
    G_base = D_base * N_base / L_base ^ 2;


    fprintf('V_base = %.8e\n', V_base);
    fprintf('N_base = %.8e\n', N_base);
    fprintf('D_base = %.8e\n', D_base);
    fprintf('L_base = %.8e\n', L_base);
    fprintf('t_base = %.8e\n', t_base);
    fprintf('E_base = %.8e\n', E_base);
    fprintf('J_base = %.8e\n', J_base);
    fprintf('mu_base= %.8e\n', mu_base);
    fprintf('G_base = %.8e\n', G_base);

    lambda2 = EPSILON_SI * V_base / (q * N_base * L_base * L_base);
    delta2 = (NIE / N_base);
    fprintf('lambda^2 = %.8e\n', lambda2);
    fprintf('delta^2 = %.8e\n', delta2);
end


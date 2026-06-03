function [] = SetBoundaryVoltage(V_Anode)
    % 无量纲化基值和物理常数
    global V_base N_base NIE
    % 器件相关参数
    global P_doping N_doping N
    global nie vt
    global psi n p phi_n phi_p


    % 根据欧姆接触设置边界条件
    % 左侧电极偏置电压为 V_Anode，右侧电极偏置电压为 0
    % 直接赋无量纲化之后的值

    %% 左边界
    phi_n(1) = V_Anode / V_base;
    phi_p(1) = V_Anode / V_base;
    if N(1) < 0
        % 左边界为 P 区
        doping = -N(1) * N_base;
        p(1) = (doping + sqrt(doping ^ 2 + 4 * NIE ^ 2)) / 2 / N_base;
        n(1) = (NIE ^ 2) / p(1) / N_base ^ 2;
        psi(1) = phi_p(1) - vt * log(p(1)/(nie)); % 使用多子浓度
    else
        % 左边界为 N 区
        doping = N(1) * N_base;
        n(1) = (doping + sqrt(doping ^ 2 + 4 * NIE ^ 2)) / 2 / N_base;
        p(1) = (NIE ^ 2) / n(1) / N_base ^ 2;
        psi(1) = phi_n(1) + vt * log(n(1)/(nie)); % 使用多子浓度
    end


    %% 右边界
    phi_n(end) = 0.0;
    phi_p(end) = 0.0;
    if N(end) < 0
        % 右边界为 P 区
        doping = -N(end) * N_base;
        p(end) = (doping + sqrt(doping ^ 2 + 4 * NIE ^ 2)) / 2 / N_base;
        n(end) = (NIE ^ 2) / p(end) / N_base ^ 2;
        psi(end) = phi_p(end) - vt * log(p(end)/(nie)); % 使用多子浓度
    else
        % 右边界为 N 区
        doping = N(end) * N_base;
        n(end) = (doping + sqrt(doping ^ 2 + 4 * NIE ^ 2)) / 2 / N_base;
        p(end) = (NIE ^ 2) / n(end) / N_base ^ 2;
        psi(end) = phi_n(end) + vt * log(n(end)/(nie)); % 使用多子浓度
    end


    fprintf('\nSet Anode voltage to %f\n', V_Anode);
end


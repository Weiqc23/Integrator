function [] = RealTimePlot() 
    global n p psi R  J_p J_n E phi_n phi_p  slotboom_n slotboom_p
    global dx n_node
    global L_base V_base E_base N_base J_base G_base q
    
    global J_n_AM J_p_AM J_n_GM J_p_GM J_n_SG J_p_SG vt

    figure(1);
    
    x = (0:1:n_node-1) * dx * L_base * 1e4; % 原始格点
    x2 = (0.5 + (0:1:n_node-2)) * dx * L_base * 1e4; % 中间格点

    % 第一行 电压电场
    subplot(4,4,1)
    plot(x, psi * V_base,'-xk')
    title("Potential")
    drawnow limitrate

    subplot(4,4,2)
    plot(x2, E * E_base, '-b')
    title("Electric Field")
    drawnow limitrate

    % 第二行 载流子浓度
    subplot(4,4,5)
    plot(x, n * N_base, '-r')
    title("Electron Concentration")
    drawnow limitrate

    subplot(4,4,6)
    plot(x, p * N_base, 'b-')
    title("Hole Concentration")
    drawnow limitrate

    subplot(4,4,7)
    hold off;
%     plot(x, p * N_base, '-b')
    semilogy(x, p * N_base, '-b')
    drawnow limitrate
    hold on;
%     plot(x, n * N_base, '-r')
    semilogy(x, n * N_base, '-r')
    title("Electron&Hole Concentration")
    drawnow limitrate

    % 第三行 电流密度
    subplot(4,4,9)
    plot(x2, J_n * J_base,'r-')
    title("Electron Current Density")
    drawnow limitrate

    subplot(4,4,10)
    plot(x2, J_p * J_base,'b-')
    title("Hole Current Density")
    drawnow limitrate

    subplot(4,4,11)
    plot(x2, (J_p + J_n) * J_base, 'm-')
    ylim([min((J_p + J_n) * J_base) - 0.5, max((J_p + J_n) * J_base) + 0.5])
    title("Total Current Density")
    drawnow limitrate
    

%     subplot(4,4,12)
%     hold off;
%     plot(x2, J_n_AM * J_base, 'r-')
%     drawnow limitrate
%     hold on
%     plot(x2, J_p_AM * J_base, 'b-')
%     title("Electron&Hole Current Density (AM)")
%     drawnow limitrate

    subplot(4,4,13)
    hold off;
    plot(x, phi_n * V_base, 'r-')
    drawnow limitrate
    hold on
    plot(x, phi_p * V_base, 'b-')
    title("Electron&Hole QF Level")
    drawnow limitrate

    subplot(4,4,14)
    plot(x, R * G_base, 'g-')
    title("Recombination Rate")
    drawnow limitrate
    
%     subplot(4,4,15)
%     hold off;
%     plot(x2, J_n_GM * J_base, 'r-')
%     drawnow limitrate
%     hold on
%     plot(x2, J_p_GM * J_base, 'b-')
%     title("Electron&Hole Current Density (GM)")
%     drawnow limitrate
    
%     subplot(4,4,16)
%     hold off;
%     plot(x2, J_n_SG * J_base, 'r-')
%     drawnow limitrate
%     hold on
%     plot(x2, J_p_SG * J_base, 'b-')
%     title("Electron&Hole Current Density (SG)")
%     drawnow limitrate




%     subplot(4,4,16)
%     hold off;
% %     plot(x, exp(-phi_n / vt), 'r-')
% %     plot(x, slotboom_n, 'r-')
%     semilogy(x, slotboom_n, 'r-')
%     drawnow limitrate
%     hold on
% %     plot(x, exp(phi_p / vt), 'b-')
% %     plot(x, slotboom_p, 'b-')
%     semilogy(x, slotboom_p, 'b-')
%     title("Slotboom variable (r-n, b-p)")
%     drawnow limitrate


end


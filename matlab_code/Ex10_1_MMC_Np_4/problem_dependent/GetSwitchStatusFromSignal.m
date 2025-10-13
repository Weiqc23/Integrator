% t 当前的t
% x 当前t的y(1*m)
% 得到每個MMC的 [i1;i2] = MMC.PvPxtrue*[Vc,iL] + MMC.rhs_vx
function [] = GetSwitchStatusFromSignal(t, x, signal_list)

matrix_change = false;

global MMC1
% x_index：3(Vc),2(IL)
component_x = x(MMC1.x_index);
% signal_index：1,2
component_signal = signal_list(MMC1.signal_index);
[MMC1, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC1, matrix_change);

global MMC2
component_x = x(MMC2.x_index);
component_signal = signal_list(MMC2.signal_index);
[MMC2, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC2, matrix_change);

global MMC3
component_x = x(MMC3.x_index);
component_signal = signal_list(MMC3.signal_index);
[MMC3, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC3, matrix_change);

global MMC4
component_x = x(MMC4.x_index);
component_signal = signal_list(MMC4.signal_index);
[MMC4, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC4, matrix_change);

global MMC5
component_x = x(MMC5.x_index);
component_signal = signal_list(MMC5.signal_index);
[MMC5, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC5, matrix_change);

global MMC6
component_x = x(MMC6.x_index);
component_signal = signal_list(MMC6.signal_index);
[MMC6, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC6, matrix_change);

global MMC7
component_x = x(MMC7.x_index);
component_signal = signal_list(MMC7.signal_index);
[MMC7, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC7, matrix_change);

global MMC8
component_x = x(MMC8.x_index);
component_signal = signal_list(MMC8.signal_index);
[MMC8, matrix_change] = Half_Bridge_GetStatusFromSignal(component_x, component_signal, MMC8, matrix_change);

%%
global A Bsw
global PuPx_sw Jac_pwl_transpose
global mode
if (matrix_change)
    mode = 1;
    Jac_pwl_transpose = (A + Bsw*PuPx_sw)';
end

end
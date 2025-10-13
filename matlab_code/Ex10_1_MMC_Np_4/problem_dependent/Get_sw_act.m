function change_interval = Get_sw_act(idx_act, idx_no_act)
global MMC1 MMC2 MMC3 MMC4 MMC5 MMC6 MMC7 MMC8
change_interval = [-1,-1];

% 处理过零的元件（更新）
for i = 1:length(idx_act)
    sw_idx = idx_act(i);
    if (sw_idx==MMC1.sw_index)
        MMC1 = Half_Bridge_Update_status(MMC1);
        change_interval = MMC1.change_interval;
        break;
    elseif (sw_idx==MMC2.sw_index)
        MMC2 = Half_Bridge_Update_status(MMC2);
        change_interval = MMC2.change_interval;
        break;
    elseif (sw_idx==MMC3.sw_index)
        MMC3 = Half_Bridge_Update_status(MMC3);
        change_interval = MMC3.change_interval;
        break;
    elseif (sw_idx==MMC4.sw_index)
        MMC4 = Half_Bridge_Update_status(MMC4);
        change_interval = MMC4.change_interval;
        break;
    elseif (sw_idx==MMC5.sw_index)
        MMC5 = Half_Bridge_Update_status(MMC5);
        change_interval = MMC5.change_interval;
        break;
    elseif (sw_idx==MMC6.sw_index)
        MMC6 = Half_Bridge_Update_status(MMC6);
        change_interval = MMC6.change_interval;
        break;
    elseif (sw_idx==MMC7.sw_index)
        MMC7 = Half_Bridge_Update_status(MMC7);
        change_interval = MMC7.change_interval;
        break;
    elseif (sw_idx==MMC8.sw_index)
        MMC8 = Half_Bridge_Update_status(MMC8);
        change_interval = MMC8.change_interval;
        break;
    end
end
% 处理过零的元件，但不是第一个过零的元件（其state没有改变）
for i = 1:length(idx_no_act)
    sw_idx = idx_no_act(i);
    if (sw_idx==MMC1.sw_index)
        MMC1.next_state_prediction = MMC1.state;
        break;
    elseif (sw_idx==MMC2.sw_index)
        MMC2.next_state_prediction = MMC2.state;
        break;
    elseif (sw_idx==MMC3.sw_index)
        MMC3.next_state_prediction = MMC3.state;
        break;
    elseif (sw_idx==MMC4.sw_index)
        MMC4.next_state_prediction = MMC4.state;
        break;
    elseif (sw_idx==MMC5.sw_index)
        MMC5.next_state_prediction = MMC5.state;
        break;
    elseif (sw_idx==MMC6.sw_index)
        MMC6.next_state_prediction = MMC6.state;
        break;
    elseif (sw_idx==MMC7.sw_index)
        MMC7.next_state_prediction = MMC7.state;
        break;
    elseif (sw_idx==MMC8.sw_index)
        MMC8.next_state_prediction = MMC8.state;
        break;
    end
end

%%
global A Bsw
global PuPx_sw Jac_pwl_transpose
Jac_pwl_transpose = (A + Bsw*PuPx_sw)';

global mode
mode = 1;

end
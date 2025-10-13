function [trigger_val, change_interval] = Get_sw_Info(sw_idx)
global MMC1 MMC2 MMC3 MMC4 MMC5 MMC6 MMC7 MMC8
trigger_val = [];
change_interval = [-1,-1];
    if (sw_idx==MMC1.sw_index)
        trigger_val = MMC1.trigger_val;
        change_interval = MMC1.change_interval;
    elseif (sw_idx==MMC2.sw_index)
        trigger_val = MMC2.trigger_val;
        change_interval = MMC2.change_interval;
    elseif (sw_idx==MMC3.sw_index)
        trigger_val = MMC3.trigger_val;
        change_interval = MMC3.change_interval;
    elseif (sw_idx==MMC4.sw_index)
        trigger_val = MMC4.trigger_val;
        change_interval = MMC4.change_interval;
    elseif (sw_idx==MMC5.sw_index)
        trigger_val = MMC5.trigger_val;
        change_interval = MMC5.change_interval;
    elseif (sw_idx==MMC6.sw_index)
        trigger_val = MMC6.trigger_val;
        change_interval = MMC6.change_interval;
    elseif (sw_idx==MMC7.sw_index)
        trigger_val = MMC7.trigger_val;
        change_interval = MMC7.change_interval;
    elseif (sw_idx==MMC8.sw_index)
        trigger_val = MMC8.trigger_val;
        change_interval = MMC8.change_interval;
    end
end
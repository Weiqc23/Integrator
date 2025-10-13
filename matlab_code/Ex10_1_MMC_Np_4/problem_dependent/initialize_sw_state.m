function [] = initialize_sw_state(y0, signal0, n)
%%
global usw_num m 
global PuPx_sw Res_sw
PuPx_sw = zeros(usw_num,m);
Res_sw = zeros(usw_num,1);

%%

% 初始化：信号[0,0]
global MMC1 MMC2 MMC3 MMC4 MMC5 MMC6 MMC7 MMC8
MMC1.signal_index = [1,2];
MMC1.signal = [-100,-100]; %只是为了启动
MMC1.sw_index = 1;
MMC1.sw_inputs_index = MMC1.inputs_index;
MMC1.next_state_prediction = [0,0];
MMC1 = Half_Bridge_Update_status(MMC1);

MMC2.signal_index = [3,4];
MMC2.signal = [-100,-100]; %只是为了启动
MMC2.sw_index = 2;
MMC2.sw_inputs_index = MMC2.inputs_index;
MMC2.next_state_prediction = [0,0];
MMC2 = Half_Bridge_Update_status(MMC2);

MMC3.signal_index = [5,6];
MMC3.signal = [-100,-100]; %只是为了启动
MMC3.sw_index = 3;
MMC3.sw_inputs_index = MMC3.inputs_index;
MMC3.next_state_prediction = [0,0];
MMC3 = Half_Bridge_Update_status(MMC3);

MMC4.signal_index = [7,8];
MMC4.signal = [-100,-100]; %只是为了启动
MMC4.sw_index = 4;
MMC4.sw_inputs_index = MMC4.inputs_index;
MMC4.next_state_prediction = [0,0];
MMC4 = Half_Bridge_Update_status(MMC4);

MMC5.signal_index = [9,10];
MMC5.signal = [-100,-100]; %只是为了启动
MMC5.sw_index = 5;
MMC5.sw_inputs_index = MMC5.inputs_index;
MMC5.next_state_prediction = [0,0];
MMC5 = Half_Bridge_Update_status(MMC5);

MMC6.signal_index = [11,12];
MMC6.signal = [-100,-100]; %只是为了启动
MMC6.sw_index = 6;
MMC6.sw_inputs_index = MMC6.inputs_index;
MMC6.next_state_prediction = [0,0];
MMC6 = Half_Bridge_Update_status(MMC6);

MMC7.signal_index = [13,14];
MMC7.signal = [-100,-100]; %只是为了启动
MMC7.sw_index = 7;
MMC7.sw_inputs_index = MMC7.inputs_index;
MMC7.next_state_prediction = [0,0];
MMC7 = Half_Bridge_Update_status(MMC7);

MMC8.signal_index = [15,16];
MMC8.signal = [-100,-100]; %只是为了启动
MMC8.sw_index = 8;
MMC8.sw_inputs_index = MMC8.inputs_index;
MMC8.next_state_prediction = [0,0];
MMC8 = Half_Bridge_Update_status(MMC8);

global mode
mode  = 1;
global A Bsw
global Jac_pwl_transpose
Jac_pwl_transpose = (A + Bsw*PuPx_sw)';

%%
MMC1 = initialize_interpolation_val(MMC1,n);
MMC2 = initialize_interpolation_val(MMC2,n);
MMC3 = initialize_interpolation_val(MMC3,n);
MMC4 = initialize_interpolation_val(MMC4,n);
MMC5 = initialize_interpolation_val(MMC5,n);
MMC6 = initialize_interpolation_val(MMC6,n);
MMC7 = initialize_interpolation_val(MMC7,n);
MMC8 = initialize_interpolation_val(MMC8,n);

%%
GetSwitchStatusFromSignal(0, y0, signal0);

end

function sw = initialize_interpolation_val(sw,n)
sw.state_judge_wrong = false;
sw.change_interval = zeros(1,2);
sw.trigger_val = zeros(n+1,1);
end
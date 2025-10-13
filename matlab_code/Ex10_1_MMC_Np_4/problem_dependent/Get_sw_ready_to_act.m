function [vec_sw_ready_to_act, some_state_judge_wrong] = Get_sw_ready_to_act(vec_x)

vec_sw_ready_to_act = [];
some_state_judge_wrong = false;

global MMC1
relative_x = vec_x(:, MMC1.x_index);
MMC1 = Half_Bridge_CheckStatus(relative_x, MMC1);
if (~isequal(MMC1.next_state_prediction, MMC1.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC1.sw_index];
end

global MMC2
relative_x = vec_x(:, MMC2.x_index);
MMC2 = Half_Bridge_CheckStatus(relative_x, MMC2);
if (~isequal(MMC2.next_state_prediction, MMC2.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC2.sw_index];
end

global MMC3
relative_x = vec_x(:, MMC3.x_index);
MMC3 = Half_Bridge_CheckStatus(relative_x, MMC3);
if (~isequal(MMC3.next_state_prediction, MMC3.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC3.sw_index];
end

global MMC4
relative_x = vec_x(:, MMC4.x_index);
MMC4 = Half_Bridge_CheckStatus(relative_x, MMC4);
if (~isequal(MMC4.next_state_prediction, MMC4.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC4.sw_index];
end

global MMC5
relative_x = vec_x(:, MMC5.x_index);
MMC5 = Half_Bridge_CheckStatus(relative_x, MMC5);
if (~isequal(MMC5.next_state_prediction, MMC5.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC5.sw_index];
end

global MMC6
relative_x = vec_x(:, MMC6.x_index);
MMC6 = Half_Bridge_CheckStatus(relative_x, MMC6);
if (~isequal(MMC6.next_state_prediction, MMC6.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC6.sw_index];
end

global MMC7
relative_x = vec_x(:, MMC7.x_index);
MMC7 = Half_Bridge_CheckStatus(relative_x, MMC7);
if (~isequal(MMC7.next_state_prediction, MMC7.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC7.sw_index];
end

global MMC8
relative_x = vec_x(:, MMC8.x_index);
MMC8 = Half_Bridge_CheckStatus(relative_x, MMC8);
if (~isequal(MMC8.next_state_prediction, MMC8.state))
    vec_sw_ready_to_act =[vec_sw_ready_to_act, MMC8.sw_index];
end


some_state_judge_wrong = MMC1.state_judge_wrong || MMC2.state_judge_wrong || MMC3.state_judge_wrong || MMC4.state_judge_wrong||...
    MMC5.state_judge_wrong || MMC6.state_judge_wrong || MMC7.state_judge_wrong || MMC8.state_judge_wrong;

end
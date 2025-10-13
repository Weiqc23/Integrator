function [t_end, y_end, some_state_judge_wrong] = detectAndInterpolate_SingleStep(t_end, y_end, t0, y0)

vec_t = [t0, t_end]';
vec_x = [y0; y_end];

%% Problem Dependent
[vec_sw_ready_to_act, some_state_judge_wrong] = Get_sw_ready_to_act(vec_x);

%%
if (some_state_judge_wrong)
    return;
end
if (isempty(vec_sw_ready_to_act))
    return;
end
%%
num_sw_ready_to_act = length(vec_sw_ready_to_act);
vec_sw_t_zero = zeros(1,num_sw_ready_to_act);

for i = 1 : num_sw_ready_to_act
    sw_idx = vec_sw_ready_to_act(i);
    [Tsw_trigger_val, Tsw_change_interval] = Get_sw_Info(sw_idx);

    % vec_sw_t_zero(i) = polyfit_find_t_zero(vec_t, Tsw.trigger_val, Tsw.change_interval);
    vec_sw_t_zero(i) = linear_find_t_zero(vec_t, Tsw_trigger_val, Tsw_change_interval);
end
%%
if length(vec_sw_ready_to_act)>1
    [t_zero, t_zero_index] = find_all_minval(vec_sw_t_zero);
    idx_act = vec_sw_ready_to_act(t_zero_index);
    idx_no_act = setdiff( sort(vec_sw_ready_to_act),  idx_act);
else
    t_zero = vec_sw_t_zero(1);
    idx_act = vec_sw_ready_to_act(1);
    idx_no_act = [];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 这里其实没有把更新的状态传出去，之前不报错是因为:
% % Nonlinear的例子里其实没有被动开关的情况；每次在GetStatusFromSignal()那里重新更新状态了。
% vec_sw_ready_to_act(idx_act).state = vec_sw_ready_to_act(idx_act).next_state_prediction;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if (~isempty(idx_no_act))
%     for i = 1:length(idx_no_act)
%         idx_stay = idx_no_act(i);
%         vec_sw_ready_to_act(idx_stay).next_state_prediction = vec_sw_ready_to_act(idx_stay).state;
%     end
% end

%% Problem Dependent
change_interval = Get_sw_act(idx_act, idx_no_act);

%%

[vec_t, vec_x] = Interpolate(vec_t, vec_x, t_zero, change_interval);
t_end = vec_t(end);
y_end = vec_x(end,:);

end



%%
function [min_val, min_indices] = find_all_minval(vec)
% 找到最小值
min_val = min(vec);

% 找到所有最小值的索引
min_indices = find(vec == min_val);
end

function t_zero = linear_find_t_zero(vec_t, vec_x, ref_interval)
    t_interval = vec_t(ref_interval);
    t1 = t_interval(1);
    t2 = t_interval(2);
    
    x1 = vec_x(ref_interval(1));
    x2 = vec_x(ref_interval(2));
    % 线性插值计算过零时间
    t_zero = t1 - x1 .* (t2 - t1) ./(x2 - x1);
end

function [vec_t, vec_y] = Interpolate(vec_t, vec_y, t_zero, ref_interval)
[vec_t, vec_y] = linear_Interpolate(vec_t, vec_y, t_zero, ref_interval);
end

function [vec_t, vec_y] = linear_Interpolate(vec_t, vec_y, t_zero, ref_interval)
    
    %step1: 线性插值
    t_interval = vec_t(ref_interval);
    t1 = t_interval(1);
    t2 = t_interval(2);
    y_interval = vec_y(ref_interval,:);
    y1 = y_interval(1,:);
    y2 = y_interval(2,:);
    y_interp = y1 + (t_zero - t1) * (y2 - y1) / (t2 - t1);

    % Step2: 截断 vec_t、vec_y、vec_z，只保留 t < t_zero 的部分
    idx_valid = find(vec_t < t_zero);
    vec_t = vec_t(idx_valid);
    vec_y = vec_y(idx_valid,:);

    % Step3: 添加 t_zero 时刻的点
    vec_t = [vec_t; t_zero];
    vec_y = [vec_y; y_interp];


end
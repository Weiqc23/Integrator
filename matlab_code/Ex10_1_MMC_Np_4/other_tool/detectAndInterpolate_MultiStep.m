function [tnode, ynode, some_state_judge_wrong] = detectAndInterpolate_MultiStep(t0, y0, tnode, ynode)

vec_t = [t0; tnode];
% p*m
vec_x = [y0; ynode];

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
    % 第vec_sw_ready_to_act(i)个MMC有过零点
    sw_idx = vec_sw_ready_to_act(i);
    % Tsw_trigger_val ：p个iL
    % Tsw_change_interval ：过零的区间（t节点）
    [Tsw_trigger_val, Tsw_change_interval] = Get_sw_Info(sw_idx);

    % vec_sw_t_zero(i) = polyfit_find_t_zero(vec_t, Tsw_trigger_val, Tsw_change_interval);
    % 计算过0的时间节点
    vec_sw_t_zero(i) = linear_find_t_zero(vec_t, Tsw_trigger_val, Tsw_change_interval);
end

%%

if length(vec_sw_ready_to_act)>1
    
    % 最早过零的时间和节点
    [t_zero, t_zero_index] = find_all_minval(vec_sw_t_zero);
    % 所有器件过零点最早的那个
    idx_act = vec_sw_ready_to_act(t_zero_index);
    % 其它过零的元件，但并不是第一个（第一个过零时，这些元件还没有过零 ）
    idx_no_act = setdiff( sort(vec_sw_ready_to_act),  idx_act);
else
    t_zero = vec_sw_t_zero(1);
    idx_act = vec_sw_ready_to_act(1);
    idx_no_act = [];
end

%% Problem Dependent
% 找到iL第一个过零点（即从负变正或从正变负的区间）所在的区间索引
% [k,k+1]:零点在k,k+1中间
change_interval = Get_sw_act(idx_act, idx_no_act);

%%
[vec_t, vec_x] = Interpolate(vec_t, vec_x, t_zero, change_interval);
% [vec_t, vec_x] = linear_Interpolate(vec_t, vec_x, t_zero, vec_sw_ready_to_act(idx_act).change_interval);

% 更新tnode和ynode只包含过零点之前的节点，最后一个节点为插值计算的过零点。
tnode=vec_t(2:end)';
ynode = vec_x(2:end,:);

end



%%
function [min_val, min_indices] = find_all_minval(vec)
% 找到最小值
min_val = min(vec);

% 找到所有最小值的索引
min_indices = find(vec == min_val);
end


function t_zero = polyfit_find_t_zero(vec_t, vec_x, ref_interval)

    global use_linear_interpolation
    use_linear_interpolation = 0;

    % 多项式拟合，获取最小正的过零点
    % 多项式阶数使用 n-1（最大插值拟合），你可以调整为更低阶数
    n = length(vec_t);
    poly_order = n-1;  % 限制阶数避免过拟合

    % 对每列做拟合并插值
    vec_y_poly_coeff = polyfit(vec_t, vec_x, poly_order);

    % 求根
    roots_detected_var = roots(vec_y_poly_coeff);
    % 提取实数根，且在 vec_t 范围内
    real_roots = roots_detected_var(imag(roots_detected_var)==0);
    valid_roots = real_roots(real_roots > vec_t(ref_interval(1)) & real_roots < vec_t(ref_interval(2)));

    % 找到最小的过零点
    t_zero = min(valid_roots);
    if(isempty(t_zero))
        t_zero = linear_find_t_zero(vec_t, vec_x, ref_interval);
        % disp('Polyfit-finding-roots may be somthing-wrong!');
    end

end

function t_zero = linear_find_t_zero(vec_t, vec_x, ref_interval)
global use_linear_interpolation
use_linear_interpolation = 1;
    t_interval = vec_t(ref_interval);
    t1 = t_interval(1);
    t2 = t_interval(2);

    x1 = vec_x(ref_interval(1));
    x2 = vec_x(ref_interval(2));
    % 线性插值计算过零时间
    t_zero = t1 - x1 .* (t2 - t1) ./(x2 - x1);
end

function [vec_t, vec_y] = Interpolate(vec_t, vec_y, t_zero, ref_interval)

global use_linear_interpolation
if (use_linear_interpolation)
    [vec_t, vec_y] = linear_Interpolate(vec_t, vec_y, t_zero, ref_interval);
else
        m = size(vec_y,2);
    vec_y_interp = zeros(1,m);
    % 多项式阶数使用 n-1（最大插值拟合），你可以调整为更低阶数
    n = length(vec_t);
    poly_order = n-1;  % 限制阶数避免过拟合
    % 对每列做拟合并插值
    for i = 1:m
        vec_y_poly_coeff = polyfit(vec_t, vec_y(:,i), poly_order);
        % 对每列 vec_z 做拟合并插值
        vec_y_interp(i) = polyval(vec_y_poly_coeff, t_zero);
    end 

    % Step 3: 截断 vec_t、vec_y、vec_z，只保留 t < t_zero 的部分
    idx_valid = find(vec_t < t_zero);
    vec_t = vec_t(idx_valid);
    vec_y = vec_y(idx_valid,:);

    % 添加 t_zero 时刻的点
    vec_t = [vec_t; t_zero];
    vec_y = [vec_y; vec_y_interp];
end

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
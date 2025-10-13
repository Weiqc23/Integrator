function [avg_abs_error, avg_rel_error] = Get_error_norm(t, t_ref, y_ref, t_plot, y_plot)
% Get_error_norm: 计算在特定时刻t，y_plot相较于y_ref的误差范数。
%
% 输入:
%   t       - (列向量) 需要计算误差的特定时间点。
%   t_ref   - (列向量) 参考数据y_ref对应的时间向量。
%   y_ref   - (矩阵) 参考数据，每一行对应t_ref中的一个时间点。
%   t_plot  - (列向量) 待比较数据y_plot对应的时间向量。
%   y_plot  - (矩阵) 待比较数据，每一行对应t_plot中的一个时间点。
%
% 输出:
%   rel_error - (列向量) 在每个指定时间点t上的相对误差(二范数)。
%   abs_error - (列向量) 在每个指定时间点t上的绝对误差(二范数)。
%
% 注意:
%   - 两个时间向量 t_ref 和 t_plot 并不需要相同，但必须都包含 t 中的所有时刻。
%   - y_ref 和 y_plot 的行数应分别与 t_ref 和 t_plot 的长度相等。
%   - y_ref 和 y_plot 的列数必须相等。

% 获取指定时间点的数量
num_t = length(t);

% 为输出结果预分配内存，以提高效率
vec_abs_error = zeros(num_t, 1);
vec_rel_error = zeros(num_t, 1);

% 遍历每一个指定的时间点
for i = 1:num_t
    tnow = t(i);

    % % 找到当前时间点在 t_ref 和 t_plot 中的索引
    % % "find(..., 1)"确保即使有重复的时间点，也只取第一个匹配项
    % idx_ref = find(t_ref == tnow, 1);
    % idx_plot = find(t_plot == tnow, 1);
    % 
    % % 检查是否在两个时间向量中都找到了该时间点
    % if isempty(idx_ref) || isempty(idx_plot)
    %     error('时间点 %f 未在某个时间向量中找到，请检查输入。', tnow);
    % end

    % --- 核心修改：使用最近邻查找 ---
    % 找到 t_ref 中与 tnow 最接近的时间点的索引
    [~, idx_ref] = min(abs(t_ref - tnow));

    % 找到 t_plot 中与 tnow 最接近的时间点的索引
    [~, idx_plot] = min(abs(t_plot - tnow));
    % ---------------------------------

    % 根据索引提取对应的 y 值
    % y_ref 和 y_plot 的每一行对应一个时间点
    y_ref_at_t = y_ref(idx_ref, :);
    y_plot_at_t = y_plot(idx_plot, :);

    % 计算绝对误差 (向量差的二范数)
    vec_abs_error(i) = norm(y_plot_at_t - y_ref_at_t);

    % 计算参考值的范数，用于计算相对误差
    norm_y_ref = norm(y_ref_at_t);

    % 计算相对误差，并处理分母为零的特殊情况
    if norm_y_ref == 0
        if vec_abs_error(i) == 0
            % 如果参考值和绝对误差都为零，则相对误差为零
            vec_rel_error(i) = 0;
        else
            % 如果参考值为零但存在误差，则相对误差为无穷大
            vec_rel_error(i) = inf;
        end
    else
        vec_rel_error(i) = vec_abs_error(i) / norm_y_ref;
    end
end

% 计算所有时刻误差的平均值
avg_abs_error = mean(vec_abs_error);
avg_rel_error = mean(vec_rel_error);

end
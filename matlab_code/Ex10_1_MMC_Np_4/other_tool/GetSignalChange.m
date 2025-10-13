function [tseries_g_change, value_g_change] = GetSignalChange(tout_ref, g_ref)
    % tout_ref: 时间序列
    % g_ref   : 每一行对应一个时刻的信号值
    
    % 初始化
    tseries_g_change = tout_ref(1);   % 第一个时间点
    value_g_change   = g_ref(1,:);    % 第一个值
    
    % 遍历检测相邻行是否不同
    for i = 2:size(g_ref,1)
        if any(g_ref(i,:) ~= g_ref(i-1,:))
            % 发现变化 -> 记录这个时间点和对应值
            tseries_g_change(end+1,1) = tout_ref(i);
            value_g_change(end+1,:)   = g_ref(i,:);
        end
    end
    
    % 保证最后一个时间点也在结果里
    if tseries_g_change(end) ~= tout_ref(end)
        tseries_g_change(end+1,1) = tout_ref(end);
    else
        value_g_change = value_g_change(1:end-1,:);
    end

end

% xtue：p*2（Vc,iL）
function sw = Half_Bridge_CheckStatus(xtrue, sw)

sw.state_judge_wrong = false;

%%
% 规定Vc,iL正方向
x = xtrue*sw.mat_PXPXtrue';
iL = x(:,end);
signs = sign(iL);% 计算符号变化（sign(iL) 会返回每个元素的符号：-1, 0, 1）
zero_crossings = find(diff(signs) ~= 0);% 找到符号变化的点（diff(signs) ~= 0 表示符号变化）


if isempty(zero_crossings)

    % iL没有变号
    sw.state = sw.next_state_prediction;
    return;
end

% 找到iL第一个过零点（即从负变正或从正变负的区间）所在的区间索引
change_index = [zero_crossings(1), zero_crossings(1)+1];
sw.next_state_prediction = Get_Half_Bridge_status(xtrue(change_index,:), sw);

if isequal(sw.state, sw.next_state_prediction)
    return;
end


% ？some_state_judge_wrong是什么含义
if ( change_index(1)==1 && abs(iL(1))<=1e-8 )
    sw.state_judge_wrong = true;
    sw.state = sw.next_state_prediction;
    return;
end

sw.change_interval = change_index;

% sw.trigger_val = iL;
% 直接写成sw.trigger_val = iL时会报错: 
% 维度 1 在左侧是固定的，但在右侧是变化的([4 x 1] ~= [:? x 1])。
% sw.trigger_val(1:length(iL)) = iL;
for i = 1:length(iL)
    sw.trigger_val(i) = iL(i);
end


end
function [A_m,tc_m] = tlgr(n)
% function [A,A_m,tc,tc_m] = tlgr(n)
% 获取 Radau IIA 节点的积分矩阵 A 和节点 tc
% 使用 if-elseif-else 结构以兼容 MATLAB Coder (codegen)

% [A,A_m,tc,tc_m] = tlgr4();
[A_m,tc_m] = tlgr4();
% if n == 3
%     [A, tc] = tlgr3();
% elseif n == 4
%     [A, tc] = tlgr4();
% elseif n == 5
%     [A, tc] = tlgr5();
% elseif n == 6
%     [A, tc] = tlgr6();
% elseif n == 7
%     [A, tc] = tlgr7();
% elseif n == 8
%     [A, tc] = tlgr8();
% elseif n == 9
%     [A, tc] = tlgr9();
% elseif n == 10
%     [A, tc] = tlgr10();
% else
%     % 如果 n 的值不受支持，则抛出错误
%     error('Unsupported n value: %d. Supported values for tlgr are 3-10.', n);
% 
% end

end
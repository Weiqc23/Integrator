disp('开始编译 MEX 文件...');

% 1. 定义输入参数的示例值
% (请确保这里的变量在你的工作区中存在且大小正确)
% 例如:


coder.varsize('tseries_g_change', [inf 1], [true false]); % 声明为可变大小
coder.varsize('value_g_change', [inf inf], [true true]); % 声明为可变大小

codegen FE -args {n, m, h_FE, t0, y0, tfinal, tseries_g_change, value_g_change};
disp('FE 编译成功！');
% [tplot_FE, yplot_FE, zplot_FE] = FE_mex(n, m, h_FE, t0, y0, tfinal, tseries_g_change, value_g_change);

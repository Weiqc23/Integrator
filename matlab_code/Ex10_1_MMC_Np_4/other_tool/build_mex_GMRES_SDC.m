% 创建一个 `build_mex.m` 文件，这会让一切都变得清晰明了，极易扩展。
% 
% --- build_mex.m ---
% 脚本文件，用于编译 GMRES_SDC 函数

disp('开始编译 MEX 文件...');

% 1. 定义输入参数的示例值
% (请确保这里的变量在你的工作区中存在且大小正确)
% 例如:


coder.varsize('tseries_g_change', [inf 1], [true false]); % 声明为可变大小
coder.varsize('value_g_change', [inf inf], [true true]); % 声明为可变大小

% 2. 将所有输入参数按顺序放入一个 cell 数组
args_spec = {n, m, h_GMRES, t0, y0, tfinal, tseries_g_change, value_g_change};

% 3. 将所有全局变量及其初始值放入一个 cell 数组
%    注意：这里我们初始化的是 set_kmax 和 set_etol，而不是 kmax 和 etol！
globals_spec = {
    'global use_linear_interpolation', 1
    % ...等等，确保列出所有用 global 关键字声明的变量...
};

% Initialization_for_build_mex(n);
% 4. 运行 codegen 命令，只调用一次 -globals
codegen KDC-ODE-Gauss-matlab\GMRES_SDC.m -args args_spec -report ...

disp('GMRES_SDC 编译成功！');

% [tplot,yplot,zplot] = GMRES_SDC_mex(n, m, h_GMRES, t0, y0, tfinal, tseries_g_change, value_g_change);


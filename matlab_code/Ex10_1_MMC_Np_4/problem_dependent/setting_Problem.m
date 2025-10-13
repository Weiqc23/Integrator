function [y0,electrical_states,inputs,outputs] = setting_Problem()
%% Device Parameters
% Ron导通电阻
% Rs关断电阻
% Vf正向压降
% Vfs反向压降
% branch支路数量
global MMC1
MMC1.Ron = 1e-3; 
MMC1.Rs = 1e6;
MMC1.Vf = 0;
MMC1.Vfd = 0;
MMC1.branch_num = 2;

global MMC2
MMC2.Ron = 1e-3;
MMC2.Rs = 1e6;
MMC2.Vf = 0;
MMC2.Vfd = 0;
MMC2.branch_num = 2;

global MMC3
MMC3.Ron = 1e-3;
MMC3.Rs = 1e6;
MMC3.Vf = 0;
MMC3.Vfd = 0;
MMC3.branch_num = 2;

global MMC4
MMC4.Ron = 1e-3;
MMC4.Rs = 1e6;
MMC4.Vf = 0;
MMC4.Vfd = 0;
MMC4.branch_num = 2;

global MMC5
MMC5.Ron = 1e-3;
MMC5.Rs = 1e6;
MMC5.Vf = 0;
MMC5.Vfd = 0;
MMC5.branch_num = 2;

global MMC6
MMC6.Ron = 1e-3;
MMC6.Rs = 1e6;
MMC6.Vf = 0;
MMC6.Vfd = 0;
MMC6.branch_num = 2;

global MMC7
MMC7.Ron = 1e-3;
MMC7.Rs = 1e6;
MMC7.Vf = 0;
MMC7.Vfd = 0;
MMC7.branch_num = 2;

global MMC8
MMC8.Ron = 1e-3;
MMC8.Rs = 1e6;
MMC8.Vf = 0;
MMC8.Vfd = 0;
MMC8.branch_num = 2;

%% state_variable_equ
global A B C D m mi mo
[A,B,C,D,y0,electrical_states,inputs,outputs] = power_analyze('ModularMultiLevelConverter');
% A:10*10
% b:10*18
m = length(A);
mi = size(D,2);

z_index = [17,18,19]; 
C = C(z_index,:);
D = D(z_index,:);
mo = length(z_index);

y0 = y0';

%%
% 输入信号u1u2
MMC1.inputs_index = [1,2];
MMC1.x_index = [1+2,2]; % [Vc, iL]
MMC1.mat_PXPXtrue = [1,0; 
                                      0,-1];% iL取流出为正

MMC2.inputs_index = [3,4];
MMC2.x_index = [2+2,2]; % [Vc, iL]
MMC2.mat_PXPXtrue = [1,0; 
                                      0,-1];% iL取流出为正

MMC3.inputs_index = [5,6];
MMC3.x_index = [3+2,2]; % [Vc, iL]
MMC3.mat_PXPXtrue = [1,0; 
                                      0,-1];% iL取流出为正

MMC4.inputs_index = [7,8];
MMC4.x_index = [4+2,2]; % [Vc, iL]
MMC4.mat_PXPXtrue = [1,0; 
                                      0,-1];% iL取流出为正

MMC5.inputs_index = [9,10];
MMC5.x_index = [5+2,1]; % [Vc, iL]
MMC5.mat_PXPXtrue = eye(2); % iL取流出为正

MMC6.inputs_index = [11,12];
MMC6.x_index = [6+2,1]; % [Vc, iL]
MMC6.mat_PXPXtrue = eye(2); % iL取流出为正

MMC7.inputs_index = [13,14];
MMC7.x_index = [7+2,1]; % [Vc, iL]
MMC7.mat_PXPXtrue = eye(2); % iL取流出为正

MMC8.inputs_index = [15,16];
MMC8.x_index = [8+2,1]; % [Vc, iL]
MMC8.mat_PXPXtrue = eye(2); % iL取流出为正

global check_interpolate
check_interpolate = true;

global get_signal
get_signal = true;

global index_sw Bsw usw_num
% 开关索引
index_sw = 1:16;
% 开关对应B矩阵
Bsw = B(:, index_sw);
usw_num = length(index_sw);

global index_ls
index_ls = [17,18];


end
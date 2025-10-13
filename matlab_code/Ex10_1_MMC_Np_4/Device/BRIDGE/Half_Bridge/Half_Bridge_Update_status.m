function sw = Half_Bridge_Update_status(sw)
%% 先确定更新状态
sw.state = sw.next_state_prediction;

%%
Vf = sw.Vf;
Vfd = sw.Vfd;
Rs = sw.Rs;
Ron = sw.Ron;

%%
% 导通电导：1/Ron + 1/Rs
% 关断电导：1/Rs
conductor_T1 = Get_Branch_Geq(sw.state(1), Ron, Rs);
conductor_T2 = Get_Branch_Geq(sw.state(2), Ron, Rs);

% sw.state(1)==1且sw.signal(1)==1：IGBT压降
% sw.state(1)==1且sw.signal(1)==0：反向二极管压降压降
Vf_T1 = Get_Branch_Vf(sw.state(1), sw.signal(1), Vf, Vfd);
Vf_T2 = Get_Branch_Vf(sw.state(2), sw.signal(2), Vf, Vfd);

%%

% Av = B[Vc;iL] + b_vf
A = [ conductor_T1,  -conductor_T2;
        1,    1];

% b = [iL+(-Vf_T1+Vf_T2)/Rs; 
%         Vc-Vf_T1-Vf_T2];

% v_vec = A\b;

B = [0,1;
        1,0];

b_vf = [(-Vf_T1+Vf_T2)/Rs; 
        -Vf_T1-Vf_T2];

% v = sw.PvPxtrue*[Vc,iL] + sw.rhs_vx
sw.PvPxtrue = A\(B*sw.mat_PXPXtrue); 
sw.rhs_vx = A\b_vf;

%% i = PiPv* v = PiPv * (sw.PvPxtrue*[Vc,iL] + sw.rhs_vx)
% 把[Vc;iL]映射到v，再把v映射到i
PiPv = Get_Bridge_PiPv(sw.Ron, sw.state);

% i = v = sw.PvPxtrue*[Vc,iL] + sw.rhs_vx
sw.PiPxtrue = PiPv*sw.PvPxtrue;
sw.rhs_ix = PiPv*sw.rhs_vx;

%%
% 用sw.PiPxtrue和sw.rhs_ix更新全局变量
Update_sw_matrix_per_device(sw);

end

function PiPv = Get_Bridge_PiPv(Ron, state)
branch_num = 2;
PiPv = zeros(branch_num);
for i = 1:branch_num
    PiPv(i,i) = Get_Branch_g(Ron, state(i));
end
end
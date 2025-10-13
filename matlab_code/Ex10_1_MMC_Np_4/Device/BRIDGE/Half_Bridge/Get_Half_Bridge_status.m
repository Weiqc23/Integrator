% 判断sw是否两个IGBT导通
% 输出[IGBT1是否导通,IGBT2是否导通]
function state_output = Get_Half_Bridge_status(xtrue, sw)

% 默认Vc大于0，默认iL不会断流
% iL取流出为正

% sw.mat_PXPXtrue [(1,0);(0,-1)]
% iL流出为正
x = xtrue*sw.mat_PXPXtrue';
Vc = x(:,1);
iL = x(:,end);

v_thresh = 2*max(sw.Vf, sw.Vfd);

epsilon = 1e-6;
if ( abs( iL(1) )>epsilon && min(Vc)>v_thresh )
    state_output = Get_state_From_iL(iL(end), sw);
else
    % iL0很小了，有断流的可能。需要用v来判断
    state_output = Get_state_From_v(xtrue(end,:), sw);
end

end

function state_output = Get_state_From_iL(iL, sw)

s1 = sw.signal(1);
s2 = sw.signal(2);

% state的排列顺序: T1 T2
if (iL>0)
    if (s1>0.5)
        % IGBT1导通
        state_output = [1,0];
    else
        % D2导通
        state_output = [0,1];
    end

elseif (iL<0)
    if (s2>0.5)
        % IGBT2导通
         state_output = [0,1];
    else
        % D1导通
        state_output = [1,0];
    end

else
    state_output = sw.state;
end

end

function state_output = Get_state_From_v(xtrue, sw)
s1 = sw.signal(1);
s2 = sw.signal(2);

%% 顺序: T1 T2 T3 T4 D5 D6
v_vec = Get_Bridge_vx(xtrue', sw);
% 求出来v已经仅仅是Ron上面的v了

state_T1 = Judge_IGBTDiode(v_vec(1), s1);
state_T2 = Judge_IGBTDiode(v_vec(2), s2);

if ( (state_T1+state_T2)~=1 )
    state_output = [0,0];
else
    state_output = [state_T1, state_T2];
end

end

function state_output = Judge_IGBTDiode(V, signal)
if (V>=0 && signal>0.5)
    % IGBT导通
    state_output = 1;
elseif ( (-V)>=0)
    % 反向二极管导通
    state_output = 1;
else
    state_output = 0;
end
end
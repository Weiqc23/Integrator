function [t_end,y_end, tnode,ynode,znode] = Collocation_onestep(m, n ,t0, y0, dt, tc, A, B, kmax, etol)

tnode=(dt*tc+t0)';  % set up the Gaussian points in [t0,t1].
t_end=dt+t0;    % t1 is the final time.

% p*m
ynode = newton_raphson_Collocation_onestep(m,y0,tnode,dt,n, A, kmax, etol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  keep the best solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nodeType
if (nodeType==1) % Gauss
    ynode=[ynode; dt*B*fright(m,ynode,tnode)+y0];
    tnode = [tnode; t_end];
end

global check_interpolate
if check_interpolate
    [tnode, ynode, some_state_judge_wrong] = detectAndInterpolate_MultiStep(t0, y0, tnode, ynode);
    znode = GetOutputs_MultiStep(tnode, ynode);
    if (some_state_judge_wrong)
        % ？重新Collocation_onestep
        [t_end,y_end, tnode,ynode,znode] = Collocation_onestep(m, n, t0, y0, dt, tc, A, B, kmax, etol);
        return;
    end
else
    znode = GetOutputs_MultiStep(tnode, ynode);
end

y_end=ynode(end,:);
t_end = tnode(end);

end
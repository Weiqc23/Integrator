function [ynew]=predictor(y0,t0,tnode,m,n)
%% 直接用y0
% ynew=ones(n,1)*y0;
% return;

% BE 或者 TR
ynew=ones(n,1)*y0;
y = y0;
for i = 1:n
    if i == 1
        tl = t0;
        t = tnode(i);
    else
        tl = tnode(i-1);
        t = tnode(i);
    end


    % y = newton_raphson_TR_onestep(y, tl, t, m);
    % y = newton_raphson_FE_onestep(y, tl, t, m);
    % fprintf("i: %d\n",i);
    y = newton_raphson_BE_onestep(y, tl, t, m);

    ynew(i,:) = y;
end

return;

%% FE 先算1轮
ynew=ones(n,1)*y0;
y = y0;
ynew(1,:) = y;

for i = 1:n
    t0 = tnode(i);
    t = tnode(i+1);
    % y_ref = newton_raphson_TR_onestep(y, t0, t, m);
    y = newton_raphson_FE_onestep(y, t0, t, m);
    ynew(i+1,:) = y;
end

end
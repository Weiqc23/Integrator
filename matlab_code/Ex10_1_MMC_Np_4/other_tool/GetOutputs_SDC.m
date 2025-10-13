function [z, znew] = GetOutputs_SDC(t1, y, t, ynew)
global C D nodeType mo
l = size(ynew,1);
znew = zeros(l,mo);
for i = 1:l
    xnow = ynew(i,:);
    tnow = t(i);
    u = Get_u_From_x(xnow', tnow);
    znow = C*xnow'+D*u;
    znew(i,:) = znow';
end
if (nodeType==1) % Gauss
    u = Get_u_From_x(y', t1);
    z = (C*y'+D*u)';
else % RadauII
    z = znew(end,:);
end
end
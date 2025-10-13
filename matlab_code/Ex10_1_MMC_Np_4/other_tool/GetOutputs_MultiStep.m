function znew = GetOutputs_MultiStep(tnew, ynew)
global C D mo
l = size(ynew,1);
znew = zeros(l,mo);
for i = 1:l
    xnow = ynew(i,:);
    tnow = tnew(i);
    u = Get_u_From_x(xnow', tnow);
    znow = C*xnow'+D*u;
    znew(i,:) = znow';
end

end
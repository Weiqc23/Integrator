function znew = GetOutputs_SingleStep(t, ynew)
global C D

% disp(['get size(ynew)=', mat2str(size(ynew))]);
% disp(ynew);
% disp(['get size(C)=', mat2str(size(C))]);
% disp(C);

u = Get_u_From_x(ynew', t);
znow = C*ynew'+D*u;
znew = znow';
end
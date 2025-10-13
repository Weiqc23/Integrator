function znew = GetOutputs_SingleStep(t, ynew)
global C D
u = Get_u_From_x(ynew', t);
znow = C*ynew'+D*u;
znew = znow';
end
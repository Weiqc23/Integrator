function u = Get_u_From_x(xnow, tnow)

global mi
u = zeros(mi,1);

global index_sw
global PuPx_sw Res_sw
u(index_sw) = PuPx_sw*xnow + Res_sw;

global index_ls
u(index_ls) = [VDCm(); VDCp()];

end

function u = VDCp()
u = 2000;
end

function u = VDCm()
u = 2000;
end
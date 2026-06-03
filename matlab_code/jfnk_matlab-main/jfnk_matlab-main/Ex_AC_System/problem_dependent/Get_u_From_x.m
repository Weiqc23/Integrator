function u = Get_u_From_x(xnow, tnow)

u = Vs(tnow);

end

function u = Vs(t)
Vrms = 735e3/sqrt(3);
f = 60;
deg_PhaseShift=0;
u=sqrt(2)*Vrms*sin(2*pi*f*t+deg_PhaseShift/180*pi);
end
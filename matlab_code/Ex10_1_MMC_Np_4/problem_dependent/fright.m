function f=fright(m,x,t)

global A B
l=length(t);
f=zeros(l,m);
for i=1:l
    tnow = t(i);
    xnow = x(i,:);
    
    u = Get_u_From_x(xnow', tnow);
     
    f(i,:) = ( A*xnow'+B*u)';
end

% f = f/1e0;

end

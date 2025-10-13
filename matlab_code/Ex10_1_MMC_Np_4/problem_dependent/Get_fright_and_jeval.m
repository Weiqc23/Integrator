function [f, dfdy] = Get_fright_and_jeval(m,x,t)

global A B 
global Jac_pwl_transpose

l=length(t);
f=zeros(l,m);
dfdy = zeros(m,m,l);

  for j=1:l
      tnow = t(j);
      xnow = x(j,:);

      u = Get_u_From_x(xnow', tnow);     
      f(j,:) = ( A*xnow'+B*u)';

      % dfdy(j,i,k) is df(i)/dy(j) at time t(k)
      dfdy(:,:,j) = Jac_pwl_transpose;
  end

  % f = f/1e0;
  % dfdy = dfdy/1e0;

end
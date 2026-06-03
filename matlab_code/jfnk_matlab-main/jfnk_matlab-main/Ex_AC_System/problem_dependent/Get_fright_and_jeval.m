function [f, dfdy] = Get_fright_and_jeval(m,x,t)

global A B 

l=length(t);
f=zeros(l,m);
dfdy = zeros(m,m,l);

  for j=1:l
      tnow = t(j);
      xnow = x(j,:);

      u = Get_u_From_x(xnow', tnow);     
      f(j,:) = ( A*xnow'+B*u)';

      % dfdy(j,i,k) is df(i)/dy(j) at time t(k)
      dfdy(:,:,j) = A';
  end

end
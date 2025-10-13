function dfdy=jeval(m,x,t)

global Jac_pwl_transpose
nn=length(t);
dfdy = zeros(m,m,nn);
  for j=1:nn
      
      % dfdy(j,i,k) is df(i)/dy(j) at time t(k)
      dfdy(:,:,j) = Jac_pwl_transpose;
  end

% dfdy = dfdy/1e0;

end

function dfdy=jeval(m,x,t)
% global time_jeval
% tic;

global A
nn=length(t);
dfdy = zeros(m,m,nn);
  for j=1:nn
      
      % dfdy(j,i,k) is df(i)/dy(j) at time t(k)
      dfdy(:,:,j) = A';
  end

% time_jeval = time_jeval+toc;
end

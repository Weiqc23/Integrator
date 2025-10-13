function ax = atv(x,t,dt,A,dfdy,n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is the matrix vector product subroutine.
%
% Input:
%   x: the input vector.
%   t: the points where we want the solution.
%   dt: the step size.
%   A: the integration matrix.
%   dfdy: the Jacobian matrix.
%   n: the number of Gaussian points.
%
% Question: (Jingfang)
%   Are the Jacobian matrix the same?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global iformulation;

switch iformulation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: for the Picard formulation, one first multiple the Jacobian
  %       matrix, then the integration matrix.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 0           % The Picard formulation.
    x=reshape(x,n+1,m);

    % --- 新增：预分配 dum0 矩阵 ---
    % 在循环开始前，根据最终的大小创建一个全零的矩阵
    % 这样 codegen 就提前知道了 dum0 的大小和类型
    dum0 = zeros(n+1, m); 

    for i=1:n+1
        dum0(i,:)=x(i,:)*dfdy(:,:,i);
    end

    dum=dt*A*dum0;

    rhs=x-dum;
    delta=updated(rhs,t,dfdy,n,m);  % Apply the preconditioner.

    ax=delta(:); % Output the results.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Note: for the yp formulation, one first multipole the integration matrix,
  %       then the Jacobian matrix.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 1          % the yp formulation.
    x=reshape(x,n+1,m);
    dum=dt*A*x;

    for i=1:n+1
      dum(i,:)=dum(i,:)*dfdy(:,:,i); %Note: the Jacobian matrix is written so we have
      % a right matrix vector product.
    end

    rhs=x-dum;
    delta=updated(rhs,t,dfdy,n,m);  % Apply the preconditioner.

    ax=delta(:); % Output the results.
  otherwise
    error('invalid iformulation'); % <-- 新代码
end

global count_atv
count_atv = count_atv+1;

return
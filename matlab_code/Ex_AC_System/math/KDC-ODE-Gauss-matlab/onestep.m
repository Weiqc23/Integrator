function [t1,y,z, t,ynew,znew, count1]=...
    onestep(m,t0,y0,dt,n,tc,kmax,gtol,etol,k0,A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  one time step for gmres-sdc
%
%  Global parameter: iformulation: yp formulation or Picard formulation
%
%  Input:
%    m: the dimension of the problem.
%    t0: starting time.
%    y0: initial value.
%    dt: time step.
%    n: number of quadrature points.
%    tc: The quadrature points.
%    gtol: Stop tolerance for GMRES.
%    etol: stop tolerance for the residual.
%    k0: GMRES(k0), dimension of the Krylov subspace.
%    A: the integration matrix. 就是paper里的积分矩阵S
%
%  Output:
% 区间 t= [t0, t1]，采用n个内部节点
%    t1: 区间末端点时刻
%    y: 区间末端点时刻的状态变量的取值
%    z: 区间末端点时刻的非状态变量的取值
%    t: 区间
%    ynew: 区间上每个节点的状态变量的取值
%    znew: 区间上每个节点的非状态变量的取值
%    count1: 每个时步的迭代次数
%
%  Note: one may change the following:
%    predictor.m: to use better predictor.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=dt*tc+t0;  % set up the Gaussian points in [t0,t1].
t1=dt+t0;    % t1 is the final time.

y0_block = ones(n,1)*y0;
delta0_block=zeros(m*n,1); % the initial guess
params=[gtol,k0]; %GMRES parameters.

global nodeType use_gmres set_etol
    etol = set_etol;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the Picard formulation.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call predictor. Currently a constant initial guess.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [ynew]=predictor(y0,t0,t,m,n);
    % disp(['one_step']);
    % disp('y0');
    % disp(num2str(y0, '%.12f '));
    % disp("predictor:");
    % disp(num2str(ynew, '%.12f '));
    % ynew = backward_euler(@f_eval,t,y0,)

    % y=ynew(n+1,:);
    if (nodeType==1) % Gauss
          temp = dt*B*fright(m,ynew,t);
          y=temp(1,:)+y0;
    else % Radau-II
          y=ynew(n,:);
    end
    
    % [f, dfdy] = Get_fright_and_jeval(m,ynew,t);
    f = fright(m,ynew,t);

    % disp(['f']);
    % disp(num2str(f, '%.12f '));

    rhs = dt*A*f - ynew + y0_block; % compute the rhs initially.
    % this is the residual defined in the paper.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize monitor vectors.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % errest=max(max(abs(rhs))); % estimate error (residue).
    
    % lerrest=errest; % used to check if residual decays or not.
    lcount=0;       % used to check if the number of nondcay error iterations.
    count1=0;           % the number of GMRES corrections.

    % global time_onestep_while
    % tic;

    while count1<kmax% && errest>etol %kmax: max number of iterations allowed.
      count1=count1+1;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  One more GMRES correction.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dfdy=jeval(m,ynew,t);  % evaluate the jacobian matrix.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Step 1: Compute the right hand side
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [epsilon]=updated(rhs,t0,t,dfdy,n,m);  % apply the preconditioner to the right hand side.
        % disp(['epsilon']);
        % disp(num2str(epsilon, '%.12f '));
        % disp(['epsilon_sdc']);
        % disp(num2str(epsilon_sdc, '%.12f '));
        % fprintf('max epsilon diff = %.3e\n', max(abs(epsilon-epsilon_sdc)));
      % epsilon=updated_mex(rhs,t,dfdy,n,m);

      if use_gmres
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Step 2: Use GMRES and solve the linear equation
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      b=epsilon(:);           % the right hand side.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % main part, solve the preconditioned linear system.
      %  input:
      %    delta: initial guess
      %    b: preconditioned right hand side.
      %    atv: matrix vector product subroutine.
      %    params: parameters for GMRES. 1-gtol. 2. gmres(k0).
      %  ouput:
      %    delta: the solution.
      %    error: error after each gmres correction.
      %    total_iters: total number of iterations.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [delta] = gmres(delta0_block, b, params,t,dt,A,dfdy,n,m);
      % [delta] = gmres_mex(delta, b, params,t,dt,A,dfdy,n,m);
      delta=reshape(delta,n+1,m); %note that gmres is designed for a row vector.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Step 3: update.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ynew=ynew+delta;                               % the new solution.

      else %只采用SDC的话
        v0 = ynew;
        errest = max(vecnorm(epsilon, 2, 1)./ vecnorm(ynew,2,1));
               
        ynew=ynew+epsilon;

        % fprintf("count1: %d\n",count1);
        % disp(['epsilon']);
        % disp(num2str(epsilon, '%.12f '));
        % disp(['epsilon_sdc']);
        % disp(num2str(epsilon_sdc, '%.12f '));
        % fprintf('max epsilon diff = %.3e\n', max(abs(epsilon-epsilon_sdc)'));
        % 
        % disp(['ynew']);
        % disp(num2str(ynew, '%.12f '));
        % fprintf("\n\n");
        % disp(['ynew_sdc']);
        % disp(num2str(ynew_sdc, '%.12f '));
        % fprintf('max y diff = %.3e\n', max(abs(ynew-ynew_sdc)'));
        % fprintf("\n\n");
        % disp("update");
        % disp("v0");
        % disp(num2str(v0, '%.12f '));
        % disp("y");
        % disp(num2str(ynew, '%.12f '));
        % disp("d");
        % disp(num2str(epsilon, '%.12f '));
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %  keep the best solution.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (nodeType==1) % Gauss
          temp = dt*B*fright(m,ynew,t);
          y=temp(1,:)+y0;
      else % RadauII
          y = ynew(n,:);
      end
            % [f, dfdy] = Get_fright_and_jeval(m,ynew,t);

      f = fright(m,ynew,t);
      rhs=dt*A*f-ynew+y0_block; %Compute the new right hand side
      
      if errest <= etol
        break;
      end



      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % monitor the errors.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % errest=max(max(abs(rhs))); % right hand side error.

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if (CheckConvergence_SDC(vec_f, ynew(2:end,:), ones(n,1)*y0, n, dt, errest, etol))
      %     break;
      % end
      if (isnan(errest))
          error('Newton iterations do not converge. Try to reduce the step size!'); % 终止并显示错误信息
      end
      % check how many non-decay iterations have happened.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % if errest>=lerrest*0.5
      %   if lcount==kmax %400
      %     % ierrmsg=1
      %     disp('Warning! Too many non-decay iterations have happened!');
      %     break;
      %   else
      %     lcount=lcount+1;
      %   end
      % else
      %   lerrest=errest; %lcount=0;
      % end

    end

    % time_onestep_while = time_onestep_while+toc;

    if count1>=kmax
         disp('Warning! Exceed the upper limit of Newton iterations!');
    end


    % disp('ynew');
    % disp(num2str(ynew, '%.12f '));
    % fprintf("\n\n");
    % 计算非状态变量z的取值
    [z, znew] = GetOutputs_SDC(t1, y, t, ynew);


return
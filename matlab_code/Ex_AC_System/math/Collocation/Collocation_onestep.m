function [t1,y,z, t,ynew,znew] = Collocation_onestep(m, t0, y0, dt, n, tc, A)

t=dt*tc+t0;  % set up the Gaussian points in [t0,t1].
t1=dt+t0;    % t1 is the final time.

ynew = newton_raphson_Collocation_onestep(m,y0,t,dt,n, A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  keep the best solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nodeType

      if (nodeType==1) % Gauss
          temp = dt*B*fright(m,ynew,t);
          y=temp(1,:)+y0;
      else % RadauII
          y=ynew(n,:);
      end

[z, znew] = GetOutputs_SDC(t1, y, t, ynew);


end
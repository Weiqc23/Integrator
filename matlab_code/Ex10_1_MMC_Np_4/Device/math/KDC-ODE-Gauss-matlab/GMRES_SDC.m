function [tplot,yplot,zplot] = GMRES_SDC(n, m, h0, t0, y0, tfinal, tseries_g_change, value_g_change)
iprob = 1;
% global iformulation
% 
% if (iformulation==0)
%     scolor='b';   % different plot color for different formulation.
% else
%     scolor='r';
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameter settings. 
%  One can change the following parameters:
%    kmax, gtol, etol, n, h0, k0.
%
%  One may also change: 
%     predictor.m: which predictor to use.
%     updated.m: which preconditioner to use.
%     chebnodes.m: which node points to use.
%     lcount==4 (see onestep.m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global set_kmax
kmax=set_kmax;            %  maximum number of outer GMRES iterations
gtol=1e-6;%1e-2;           %  Tolerance for GMRES call(original 1e-11)
                     %  for linear problems, this is
                     %  the total stopping procedure
                     %  For nonlinear problems, it is the
                     %  tolerance for each Newton iteration
global set_etol
etol=set_etol;% 1e-14;          %  Total error tolerance for SDC step(original 1e-13)
                     %  the residual in Picard form must be 
                     %  smaller than etol
    % nk = 1; 
    % figure(3*nk -2); % output figures.
    % subplot(2,2,k);
    % figure(3*nk-1);
    % subplot(2,2,k);
    % figure(3*nk);
    % subplot(2,2,k);
    % figure(10);
    
    % loop for k0, full GMRES: k0 = (n+1)*m
    k0 = (n+1)*m; %20; %(n+1)*m; % k0=max(n+4,81);  % Number of terms in Krylov subspace for the 
                       % restarted GMRES(k0).
                       % Note: instead of GMRES(k0), BiCGStab or
                       %  TFQMR may be applied, which should be more 
                       %  efficient and require less storage.

    %%%%%%%%%%%%%%%%%%%%%%Main Subroutine%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [tplot,yplot,zplot]=sdcgmres(iprob,m,t0,tfinal,y0,h0,n,kmax,gtol,etol,k0, tseries_g_change, value_g_change);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The required parameters of sdcgmres are explained below.
% Output
%    ysolfin: the final solution.
%    res: residue output by GMRES.
%    indres (row vector) index for res : res(indres(i-1)+1:indres(i)), 
%        residues outputs from gmres' during the i_th time step
%    errrhs: the error after each GMRES acceleration.    
%    errtrue: the true error after each gmres acceleration.
%    inderr: (row vector) index for err: err(inderr(i-1)+1:inderr(i)),
%        errors outputs from gmres' during the i_th time step
%    iter (row vector) iteration numbers outputed by each gmres.
%    ierrmsg: error message index.
% 
% Input:
%    iprob: problem index.
%    m: size of the problem
%    t0: initial time
%    tfinal:  fintal time
%    y0: (row vector) initial value for y
%    h0:  step size
%    n: number of grid points used for each time step
%    kmax: maximal number that the gmres is done
%    gtol: tolerance for gmres
%    etol: tolerance for err/res. Note that this is for right hand side.
%    k0: gmresk selecting parameter
%
    
end
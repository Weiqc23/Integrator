clear all;
dt=1;
n=5;
addpath '../../KDC-ODE-Radau-matlab/nodes/'
%tgauss10;
tlgr5
t=dt*tc;
lambda=-100;
y0=1;
yp=inv(eye(n+1)-dt*lambda*A)*(lambda* y0*ones(n+1,1)-lambda*(cos(t))'-(sin(t))');


% the error for the solution yp.
err0=yp+sin(t)'

% the error in the solution y.
err1=(cos(t))'-dt*A*yp-y0*ones(n+1,1)

% the error at the last point.
err2=cos(dt)-y0-dt*B*yp

% << NumericalDifferentialEquationAnalysis`
% GaussianQuadratureWeights[10, 0, 1, 20]
% compared with Mathematica file, everything seemed correct!
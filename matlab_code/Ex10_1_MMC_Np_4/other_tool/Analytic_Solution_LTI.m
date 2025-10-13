function x = Analytic_Solution_LTI(A,B,x0,t0,t)
%ANALYTIC_SOLUTION_LTI Analytic solution of dx/dt = A x + B u(t)
%
%   x = Analytic_Solution_LTI(A,B,x0,t,t0)
%
% Inputs:
%   A   - n x n constant matrix
%   B   - n x m constant matrix
%   x0  - n x 1 initial state vector at time t0
%   t   - scalar time at which to evaluate the solution
%   t0  - initial time (scalar)
%
% Requires:
%   A user-supplied function u(t), defined separately as:
%       function val = u(t)
%           % Example: scalar input
%           val = sin(t);
%           % or vector input, e.g.
%           % val = [sin(t); cos(t)];
%       end
%
% Output:
%   x   - n x 1 state vector at time t

    % Homogeneous solution
    Phi = expm(A * (t - t0));
    x_h = Phi * x0;

    % Particular solution (numerical quadrature)
    integrand = @(s) expm(A * (t - s)) * B * u_res(s);
    if t == t0
        x_p = zeros(size(x0));
    else
        x_p = integral(@(s) integrand(s), t0, t, 'ArrayValued', true);
    end

    % Total solution
    x = x_h + x_p;
end

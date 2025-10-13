function if_converged = CheckConvergence_SDC(vec_f, vec_y, vec_y0, n, h, errest, set_etol)

vec_f0 = ones(n,1)*vec_f(1,:);
vec_ft = vec_f(2:end,:);
vec_residual_TR = h/2*(vec_ft+vec_f0)-(vec_y-vec_y0);

tol = 1e-5;
if_converged = ( max(max(abs(vec_residual_TR))) < tol );

% if_converged = ( errest < set_tol );

end
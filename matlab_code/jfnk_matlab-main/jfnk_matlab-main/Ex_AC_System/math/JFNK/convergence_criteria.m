function [is_converged] = convergence_criteria(d,y,tol);

    value = relative_norm(d,y);

    is_converged = (value <= tol);

end
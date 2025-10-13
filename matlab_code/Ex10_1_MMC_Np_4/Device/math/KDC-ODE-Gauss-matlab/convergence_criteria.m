function [is_converged] = convergence_criteria(d,y,tol)

    ratio = vecnorm(d,2,1)./vecnorm(y,2,1);
    value = max(ratio);

    is_converged = (value < tol);

end
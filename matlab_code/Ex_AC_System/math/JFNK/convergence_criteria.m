function [is_converged] = convergence_criteria(d,y,tol);
    % disp("is_converged");
    % disp("d");
    % disp(num2str(d, '%.12f '));
    % disp("y");
    % disp(num2str(y, '%.12f '));
    
    value = relative_norm(d,y);
    % disp("value");
    % disp(num2str(value, '%.12f '));
    is_converged = (value <= tol);

end
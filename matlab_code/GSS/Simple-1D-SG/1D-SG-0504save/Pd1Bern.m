function [dBdx] = Pd1Bern(x)
    % Bernoulli 函数的导数

    BP0_DBERN = -4.117169706049845e+01;
    BP1_DBERN = -3.742994775023696e+01;
    BP2_DBERN = -1.557983451962405e-02;
    BP3_DBERN = 8.928625524999144e-03;
    BP4_DBERN = 3.742994775023696e+01;
    BP5_DBERN = 7.451332191019408e+02;

    if x <= BP0_DBERN
        dBdx = (-1.0);
    elseif x <= BP1_DBERN
        dBdx = ((1.0 - x) * exp(x) - 1.0);
    elseif x < BP2_DBERN
        y = exp(x);
        z = y - 1.0;
        dBdx = (((1.0 - x) * y - 1.0) / (z * z));
    elseif x <= BP3_DBERN
        dBdx = (-0.5 + x / 6.0 * (1.0 - x * x / 30.0));
    elseif x < BP4_DBERN
        y = exp(-x);
        z = 1 - y;
        dBdx = (((1.0 - x) * y - y * y) / (z * z));
    elseif x < BP5_DBERN
        y = exp(-x);
        dBdx = ((1.0 - x) * y - y * y);
    else
        dBdx = (0.0);
    end
end


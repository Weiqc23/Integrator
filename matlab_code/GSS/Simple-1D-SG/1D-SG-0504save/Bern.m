function [Bx] = Bern(x)
    % Bernoulli function

    BP0_BERN = -3.742994775023696e+01;
    BP1_BERN = -1.983224379137254e-02;
    BP2_BERN = 2.177550998053050e-02;
    BP3_BERN = 3.742994775023696e+01;
    BP4_BERN = 7.451332191019408e+02;

    if x <= BP0_BERN
        Bx = (-x);
    elseif x < BP1_BERN
        Bx = (x / (exp(x) - 1.0));
    elseif x <= BP2_BERN
        Bx = (1.0 - x / 2.0 * (1.0 - x / 6.0 * (1.0 - x * x / 60.0)));
    elseif x < BP3_BERN
        y = exp(-x);
        Bx = ((x * y) / (1.0 - y));
    elseif x < BP4_BERN
        Bx = (x * exp(-x));
    else
        Bx = 0;
    end
end
function [do_stop] = stopping_criteria(t_final, t0)

    T_SMALL = 1e-15;
    EPS = 1e-20;
    t = t0+EPS;

    if ((t_final >= T_SMALL) || (t_final == 0))
        do_stop = (t>=t_final);
    else
        do_stop = ((t_final - t) / t_final) <= 1e-5;
    end

end
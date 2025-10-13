function i_vec = Bridge_ix(xtrue, sw)
i_vec = sw.PiPxtrue*xtrue + sw.rhs_ix;
end
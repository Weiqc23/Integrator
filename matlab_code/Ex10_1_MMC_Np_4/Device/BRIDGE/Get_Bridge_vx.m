function v_vec = Get_Bridge_vx(xtrue, sw)
v_vec = sw.PvPxtrue*xtrue + sw.rhs_vx;
end
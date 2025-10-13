function [] = Update_sw_matrix_per_device(sw)
    global PuPx_sw Res_sw
    PuPx_sw(sw.sw_inputs_index, sw.x_index) = Bridge_PiPx(sw);
    Res_sw(sw.sw_inputs_index) = sw.rhs_ix;
end
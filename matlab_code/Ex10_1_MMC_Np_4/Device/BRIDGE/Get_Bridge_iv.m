function i_vec = Get_Bridge_iv(v_vec, sw)
branch_num = sw.branch_num;
i_vec = zeros(branch_num,1);

for i = 1:branch_num
   i_vec(i) = Get_branch_i(v_vec(i), sw.Ron, sw.state(i));
end
end

function i = Get_branch_i(v, Ron, state_input)
if (state_input==0)
    i = 0;
else
    i = v/Ron;
end

end
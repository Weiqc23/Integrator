function Geq = Get_Branch_Geq(state_single_branch, Ron, Rs)
if (state_single_branch==0)
    Geq = 1/Rs;
else
    Geq = 1/Ron+1/Rs;
end
end
function Vf = Get_Branch_Vf(state_single_branch, signal_single_branch, Vf, Vfd)
if (state_single_branch==1)
    if (signal_single_branch>0.5)
    else
        Vf = -Vfd;
    end
end
end
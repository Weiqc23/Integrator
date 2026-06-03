function [f_total_] = Ftotal()
    % 将三个 rhs 合并成为一个大 rhs
    f_psi_ = Fpsi();
    f_n_ = Fn();
    f_p_ = Fp();
    f_total_ = [f_psi_; f_n_; f_p_];
end


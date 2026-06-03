function [step_size] = BankRoseDamping()
    global n_node
    global J_mat fx delta_x
    global psi phi_n phi_p n p

    % Bank Rose damping
    psi_save = psi;
    n_save = n;
    p_save = p;
    
    fx = Ftotal();
    norm_1 = norm(fx, 2);

    K = 0.0;
    delta = 1e-4;
    norm_2 = 1.0;

    % fprintf('In damping, initial norm = %.3e\n', norm_1);
    while 1
        step_size = 1.0 / (1.0 + K * norm_1);
        if step_size < 1e-10
            fprintf('Bank Rose damping fail\n');
            return
        end

        psi(2:end-1) = psi(2:end-1) + step_size * delta_x(1:n_node-2);
        n(2:end-1) = n(2:end-1) + step_size * delta_x((n_node-2)+1:2*(n_node-2));
        p(2:end-1) = p(2:end-1) + step_size * delta_x(2*(n_node-2)+1:3*(n_node-2));

        fx = Ftotal();
        norm_2 = norm(fx, 2);
        % fprintf('Step size = %.3e, norm after step = %.3e\n', step_size, norm_2);
        if isnan(fx)
            fprintf('NaN occurs while damping\n');
        end
        if ((1.0 - norm_2 / norm_1) / step_size < delta) || isnan(norm_2)
            % Damping 失败
            if (abs(K) < 1e-10) % K 不能严格 == 0
                K = 1;
            else
                K = 10 * K;
                psi = psi_save;
                n = n_save;
                p = p_save;
            end
            continue;
        else
            % Damping 成功
            K = K / 10;
            break;
        end
    end

    % Damping 成功

    % 真正的修正在外面进行吧，这里把 psi phi 都还原回去
    psi = psi_save;
    n = n_save;
    p = p_save;

end


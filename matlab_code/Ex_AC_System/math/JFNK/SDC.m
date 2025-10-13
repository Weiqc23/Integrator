% function[y_sdc,y_sdc_m,Y,D,D_m,is_converged] = SDC(f_eval, t0, t, t_m, y0, y_old, y_old_m, S,S_m, S_p, S_p_m, n_iter_max_sdc, be_tol,sdc_tol)
function[y_sdc,Y,D,is_converged] = SDC(f_eval, t, y0, y_old, S, S_p, n_iter_max_sdc, be_tol,sdc_tol)
    
    Y = {};
    % Y_m = {};
    D = {};
    % D_m = {};
    n = length(t);
    % n_m = length(t_m);
    m = length(y0);
    y_sdc = y_old;
    y0_block = ones(n,1)*y0';
    % disp(['SDC:size(y_old)=', mat2str(size(y_old)), ', numel(y_old)=', num2str(numel(y_old)), ', n=', num2str(n), ', m=', num2str(m)]);

    % y_sdc = reshape(y_old, [n, m]);
    % y_sdc_m = reshape(y_old_m, [n_m, m]);
    
    k = 0;
    is_converged = false;

    F = zeros(n, m);
    % F_m = zeros(n_m,m);
    % disp('y_old');
    % disp(y_old);
    % disp('y_sdc');
    % disp(y_sdc);
    while ((~is_converged) && (k < n_iter_max_sdc))
        
        Y{end+1} = y_sdc;
        % Y_m{end+1} = y_sdc_m;
        F = zeros(n,m);
        % F_m = zeros(n_m,m);
        for (i = 1:n)
            yi = y_sdc(i,:);
            yi = yi(:);
            F(i,:) = f_eval(t(i), yi).';
            % F(i,:) = transpose(f_eval(t(i), transpose(y_sdc(i,:))))
        end
        % for (i = 1:n_m)
        %     yi_m = y_sdc_m(i,:);
        %     yi_m = yi_m(:);
        %     F_m(i,:) = f_eval(t_m(i), yi_m).';
        %     % F(i,:) = transpose(f_eval(t(i), transpose(y_sdc(i,:))))
        % end
        [y_sdc, d] = SDC_sweep(@f_eval,t,y0,y_sdc,F,S,S_p,be_tol);
        % [y_sdc_m, d_m] = SDC_sweep(@f_eval,t0,t_m,y0,y_sdc_m,F_m,S_m,S_p_m,be_tol);
        is_converged = convergence_criteria(d,Y{end},sdc_tol);
        % if is_converged
        %     f = fright(m,y_sdc,t);
        %     f_eval = f_eval(t(1),y_sdc(1,:));
        %     rhs=S*f-y_sdc+y0_block; %Compute the new right hand side
        %     % disp("t");
        %     % disp(num2str(t,'%.12f'));
        %     % disp("y_sdc")
        %     % disp(num2str(y_sdc,'%.17f'));
        %     % disp("f");
        %     % disp(num2str(f,'%.17f'));
        %     % disp("f_eval");
        %     % disp(num2str(f_eval(:)','%.12f'));
        %     % 
        %     % disp("rhs");
        %     % disp(num2str(rhs, '%.12f '));
        %     % fprintf("\n");
        % end
        D{end+1} = d;
        % D_m{end+1} = d_m;
        k = k+1;

    end
    % disp(['sdc']);
    % % disp('y0');
    % % disp(num2str(y0', '%.12f '));
    % % disp('y_old');
    % % disp(num2str(y_old, '%.12f '));
    % disp('y_sdc');
    % disp(num2str(y_sdc, '%.12f '));
    % disp('is_con');
    % disp(is_converged);
    % fprintf("\n\n");
    
end


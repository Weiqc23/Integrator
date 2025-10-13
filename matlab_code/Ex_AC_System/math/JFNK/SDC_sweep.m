function [y,d] = SDC_sweep(f_eval,t,y0,y_old,F,S,S_p,be_tol)
    % disp('y0');
    % disp(num2str(y0', '%.12f '));
    % disp('y_old');
    % disp(num2str(y_old, '%.12f '));
    % disp('F');
    % disp(num2str(F, '%.12f '));

    n = length(t);
    m = length(y0);
    y = zeros(n,m);
    d = zeros(n,m);
    h = 0;
    
    for (i = 1:n)
        % if (i==1)
        %     h = t(i)-t0;
        % else
        %     h = t(i)-t(i-1);
        % end
        h = S_p(i,i);
    
        w = zeros(n,1);
        w(i) = h;
        % do_skip = all(S_p{i} == 0, 'all');
        do_skip = all(S_p(i, :) == 0);

        if ((i==1)&&do_skip) 
            y(1,:) = y0';
        else
            rhs = zeros(m,1);

            if ((i==1) && (~do_skip)) 
                rhs = y0 + ((S(i,:) - w') * F).';
            else
                rhs = y(i-1,:)' + ((S(i,:) - S(i-1,:) - w') * F).';
                % disp(['SDC_sweep:size(y(i-1,:))=', mat2str(size(y(i-1,:)'))]);
                % disp(['SDC_sweep:size(((S(i,:) - S(i-1,:) - w))=', mat2str(size( (S(i,:) - S(i-1,:) - w') ))]);
                % disp(['SDC_sweep:size(F)=', mat2str(size(F))]);
                % disp(['SDC_sweep:size(rhs)=', mat2str(size(rhs)), ', numel(rhs)=', num2str(numel(rhs))]);
            end
            % disp(['sdc_sweep_rhs']);
            % fprintf('列号: ');
            % for col = 1:length(rhs)
            %     fprintf('%19d', col);
            % end
            % fprintf('\n');
            % disp(num2str(rhs', '%.12f '));

            [y_temp, d_temp] = SDC_node(@f_eval,t(i),h,rhs, transpose(y_old(i,:)),be_tol);

            y(i,:) = y_temp;
            d(i,:) = d_temp;

        end
    
    end
    

    % disp("sdc_sweep");
    % disp("y_old");
    % disp(num2str(y_old, '%.12f '));
    % disp("y");
    % disp(num2str(y, '%.12f '));
    % disp("d");
    % disp(num2str(d, '%.12f '));
    % fprintf("\n\n")
end
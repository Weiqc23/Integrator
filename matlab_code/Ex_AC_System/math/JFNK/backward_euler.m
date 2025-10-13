function [y] = backward_euler(f_eval,t,y0,S_p,be_tol)
    n = length(t);
    m = length(y0);
    y = zeros(n,m);
    % disp("start backward_euler")
    for i = 1:n
        h = S_p(i,i);
        % fprintf("i: %d\n",i);
        if (i==1 && h==0)
            y(1,:) = y0';
        else 
            rhs = zeros(m);
            if (i==1)
                rhs = y0;
            else
                rhs = y(i-1,:)';
            end
            % 求解y: y - h*f(t(i),y) - rhs = 0;
            y(i,:) = backward_euler_node(@f_eval, t(i), h, rhs, rhs, be_tol);
            
            % fprintf("i: %d\n",i);
            % disp("rhs");
            % disp(num2str(rhs', '%.12f '));
            % disp("y(i,:)");
            % disp(num2str(backward_euler_node(@f_eval, t(i), h, rhs, rhs, be_tol), '%.12f '));
        end
            
    end
    % disp('y:');
    % disp(num2str(y, '%.12f '));
end

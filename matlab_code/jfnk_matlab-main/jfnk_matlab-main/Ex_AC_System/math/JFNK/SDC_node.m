function [y_new,d] = SDC_node(f_eval,t,h,rhs,y_old,be_tol)
    
    % disp(['SDC_node:size(y_old)=', mat2str(size(y_old)), ', numel(y_old)=', num2str(numel(y_old))]);
    % disp(['SDC_node:size(rhs)=', mat2str(size(rhs)), ', numel(rhs)=', num2str(numel(rhs))]);
    % 
    % disp("SDC_node, rhs");
    % disp(rhs);
    % disp("SDC_node, y_old");
    % disp(y_old);

    % y_new_line = backward_euler_node_sdc(@f_eval,t,h,rhs,y_old,be_tol);
    y_new = backward_euler_node_sdc(@f_eval,t,h,rhs,y_old,be_tol);
    d = y_new - y_old;
    % d_line = y_new_line - y_old;
    % norm_diff = norm(y_new - y_new_line);
    % fprintf('max diff = %.3e\n', max(abs(y_new-y_new_line)));
    % disp("sdc_node");
    % disp("d");
    % disp(num2str(d', '%.12f '));
    % disp("d_line");
    % disp(num2str(d_line', '%.12f '));
    % 
    % 
    % disp("y_new");
    % disp(num2str(y_new', '%.12f '));
    % disp("y_line");
    % disp(num2str(y_new_line', '%.12f '));
    

end
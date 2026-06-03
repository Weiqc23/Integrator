function f = f_eval(t, y)
    global A B C D m mi mo 
    % disp(['size(A)=', mat2str(size(A)), ', size(B)=', mat2str(size(B)), ', size(y)=', mat2str(size(y))]);
    y = y(:);
    u = Get_u_From_x(y,t);
    f = A*y + B*u;
    f = f(:);

    % disp("f");
    % disp(num2str(f','%.12f'));
    % disp("f_eval");
    % disp("y")
    % disp(num2str(y','%.12f'));
    % fprintf("\n");
end
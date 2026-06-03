function f=fright(m,x,t)

global A B
l=length(t);
f=zeros(l,m);
for i=1:l
    tnow = t(i);
    xnow = x(i,:);
    
    u = Get_u_From_x(xnow', tnow);
     
    f(i,:) = ( A*xnow'+B*u)';
    % disp("f");
    % disp(num2str(f','%.12f'));
    % disp("f_eval");
    % disp("y")
    % disp(num2str(xnow','%.12f'));
    % fprintf("\n");
end

end

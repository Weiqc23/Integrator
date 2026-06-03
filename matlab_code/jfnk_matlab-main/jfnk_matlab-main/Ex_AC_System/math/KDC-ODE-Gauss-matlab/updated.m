function [delta]=updated(rhs,t0,t,dfdy,n,m)
global count_updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given right hand side, this subroutine returns
% (1-h*tilde{A}*alpha)^(-1)*rhs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 更快：不初始化所有元素，交给编译器优化
delta = coder.nullcopy(zeros(n,m));
% delta(1,:) = 0;

global ExImplicit
if (ExImplicit) % 隐式
    Mat_unit = eye(m);
    for i = 1:n

        if (i==1) 
            dt = t(i)-t0;
        else 
            dt = t(i)-t(i-1);
        end

        Jac = dfdy(:,:,i);
        
        if i==1
            delta(i,:) = rhs(i,:)/( Mat_unit - dt*Jac ); %隐式
        else
            delta(i,:) = (delta(i-1,:)+rhs(i,:)-rhs(i-1,:))/( Mat_unit - dt*Jac ); %隐式
        end

        % disp("update");
        % disp("d");
        % disp(num2str(delta(i,:), '%.12f '));
        % fprintf('列号: ');
        % for col = 1:length(rhs)
        %     fprintf('%17d', col);
        % end
        % fprintf('\n');
        % disp(['update_rhs']);
        % disp(num2str(rhs, '%.12f '));
    end
    
else
    for i = 1:n
        dt = t(i+1)-t(i);
        Jac = dfdy(:,:,i);
        delta(i+1,:) = dt*delta(i,:)*Jac + delta(i,:) + rhs(i+1,:)-rhs(i,:); %显式
    end

end

% % Mat_unit = eye(m);
% for i = 1:n
%     dt = t(i+1)-t(i);
%     Jac = dfdy(:,:,i);
%     delta(i+1,:) = dt*delta(i,:)*Jac + delta(i,:) + rhs(i+1,:)-rhs(i,:); %显式
%     % delta(i+1,:) = (delta(i,:)+rhs(i+1,:)-rhs(i,:))/( Mat_unit - dt*Jac ); %隐式
% end

count_updated = count_updated+1;

end

%%
% delta(i+1,:) = Get_delta_Explicit(delta(i,:), rhs(i+1,:), rhs(i,:), t(i+1)-t(i), dfdy(:,:,i));
function delta = Get_delta_Explicit(delta0, rhs, rhs0, dt, dfdy)
delta=dt*delta0*dfdy + delta0 + rhs-rhs0;
end

% delta(i+1,:) = Get_delta_Implicit(delta(i,:), rhs(i+1,:), rhs(i,:), t(i+1)-t(i), dfdy(:,:,i+1), m);
function delta = Get_delta_Implicit(delta0, rhs, rhs0, dt, dfdy, m)
Mat_unit=eye(m);
delta=(delta0+rhs-rhs0)/( Mat_unit - dt*dfdy ); 
end

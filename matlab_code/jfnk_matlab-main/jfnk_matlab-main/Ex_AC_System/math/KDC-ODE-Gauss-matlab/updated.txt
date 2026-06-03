function delta=updated(rhs,t,dfdy,n,m)
global count_updated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given right hand side, this subroutine returns
% (1-h*tilde{A}*alpha)^(-1)*rhs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delta=zeros(n+1,m);
% 更快：不初始化所有元素，交给编译器优化
delta = coder.nullcopy(zeros(n+1,m));
delta(1,:) = 0;

global ExImplicit
if (ExImplicit) % 隐式
    Mat_unit = eye(m);
    for i = 1:n
        dt = t(i+1)-t(i);
        Jac = dfdy(:,:,i);
        delta(i+1,:) = (delta(i,:)+rhs(i+1,:)-rhs(i,:))/( Mat_unit - dt*Jac ); %隐式
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

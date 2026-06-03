classdef Params
    properties
        n_nodes
        m
        node_type
        t
        S
        S_p
        S_b
        spectral_radius
    end
    
    methods
        function obj = Params(n, m, node_type)
            % 构造函数
            obj.n_nodes = n;
            obj.m = m;
            obj.node_type = node_type;
            if (node_type==1)% Gauss
                [obj.S, obj.S_b, obj.t] = chebnodes(n);  %set up the integration matrix, and construct tc.
            else %Radau-II
                [obj.S, obj.t] = tlgr(n);  %set up the integration matrix, and construct tc.
                obj.S_b = zeros(n);
            end
            [obj.S_p] = backward_euler_matrix(obj.t);
            obj.spectral_radius = spectral_radius(node_type, obj.S, obj.S_p);
        end
    end
end

function radius = spectral_radius(node_type, S, S_p)
   
    n = size(S, 1);

    C = [];

    if node_type==1
        C = eye(n-1) - inv(S_p(1:n-1, 1:n-1)) * S(1:n-1, 1:n-1);
    else
        C = eye(n) - inv(S_p) * S;
    end

    % 特征值分解
    e_vals = eig(C);

    % 计算谱半径（最大特征值模）
    radius = max(abs(e_vals));
end

function [S_p] = backward_euler_matrix(t)
     do_end_point = (abs(t(end) - 1.0) < 1e-10);
    
     if ~do_end_point
        tau = [t, 1.0];
    else
        tau = t;
     end

     n = length(tau);

     S_p = zeros(n,n);

     for i = 1:n
        if i == 1
            delta_t = tau(1) - 0;
            S_p(:, 1) = delta_t;
        else
            delta_t = tau(i) - tau(i-1);
            S_p(i:end, i) = delta_t;
        end
     end

end
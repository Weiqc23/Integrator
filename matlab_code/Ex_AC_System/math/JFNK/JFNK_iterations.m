function[y,Y_list,D_list,is_converged] = JFNK_iterations(f_eval,t,y0,y_init,S,S_p,n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)

    y = y_init;
    k_sdc = length(t)+1;
    m = length(y0);
    
    Y_list = {};
    D_list = {};
    
    k=0;
    is_converged = false;
    
    while ((~is_converged) && (k<n_iter_max_newton))
        
        [y,Y,D,is_converged] = SDC(@f_eval,t,y0,y,S,S_p, k_sdc, be_tol, sdc_tol);
    
        if (is_converged)
            Y_list = {Y_list{:},Y{1:end}};
            D_list = {D_list{:},D{1:end}};
        else
            Y_list = {Y_list{:},Y{1:end-1}};
            D_list = {D_list{:},D{1:end-1}};
            y = Y{end};
            
            yT = transpose(y);
            DT = cell(size(D));
            for (i = 1:length(D))
                DT{i} = transpose(D{i});
            end
    
            A = cell(1,k_sdc-1);
            for (i = 1:k_sdc-1)
                A{i} = DT{i+1}-DT{i};
            end
        
            M = cell(1,k_sdc-1);
            for (i = 1:k_sdc-1)
                M{i} = DT{i};
            end
    
            for (i = 1:m) 
                B = zeros(length(A), size(A{1}, 2));
                for (j = 1:length(A)) 
                     B(j, :) = A{j}(i,:);
                end
    
                BT = transpose(B);
                rhs = -transpose(DT{end}(i,:));
    
                c = BT\rhs;
    
                V = zeros(length(M),size(M{1},2));
                
                for (j = 1:length(M)) 
                    V(j,:) = M{j}(i,:);
                end
    
                VT = transpose(V);
                dy = VT*c;
    
                yT(i,:) = yT(i,:) + transpose(dy);
            end
            y = transpose(yT);
        end

        k = k+1;
    end

end
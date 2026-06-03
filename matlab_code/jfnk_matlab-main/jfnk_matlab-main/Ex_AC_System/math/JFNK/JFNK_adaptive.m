function [t_all,y_all,z_all] = JFNK_adaptive(f_eval, t_init, t_final, dt, p , y0, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol)

    do_stop = false;
    y0 = y0(:);
    y_all = y0';
    t_all = t_init(:);
    
    z_all = GetOutputs_SingleStep(t_init, y0');
    Y_all = {};
    D_all = {};
    v0 = y0;
    t0 = t_init;
    j = 0;

    % fprintf("do_stop: %d\n", do_stop); 
    % fprintf("t0: %.6f\n", t0);
    % fprintf("t_final: %.6f\n", t_final);
    %%
% 初始化全局计时器变量
    time_JFNK_onestep = 0;
    time_save_results = 0;
    
    % 改为更稳的计时方式
    startTime = datetime('now');
    while ((t0 < t_final) && (~do_stop))
        if (t0 + dt > t_final)
            dt = t_final - t0;
        end
        
        tic;
            [t,y,Y,D]=JFNK_uniform(@f_eval, t0, t0+dt, p, v0, n_iter_max_newton, n_iter_max_sdc, be_tol, sdc_tol);
            [z_end,z] = GetOutputs_SDC(t(end), y(end)', t, y);
            time_JFNK_onestep = time_JFNK_onestep + toc;
        
        tic;

        if (p.t(end)==1)

            y_all = [y_all;y(1:end,:)];
            z_all = [z_all;z(1:end,:)];
            t_all = [t_all;t(1:end)];
            
        else
            y_all = [y_all;y(:,:)];
            z_all = [z_all;z(:,:)];
            t_all = [t_all;t(:)]
        end

        time_save_results = time_save_results + toc;
        
        
        v0 = y(end,:)';
        t0 = t0 + dt;
        Y_all{end+1} = Y;
        D_all{end+1} = D;

        % if (mod(j,100)==0)
        %     fprintf("j: %d\n",j);
        %     fprintf("t: %.6f / %.6f\n", t0, t_final);
        %     disp(['SDC:size(y)=', mat2str(size(y)), ', numel(y)=', num2str(numel(y))]);
        %     disp(y);
        %     fprintf("\n\n");
        %     do_stop = stopping_criteria(t_final,t0);
        % end
        
        j = j+1;
    end
    
    % fprintf('y形状: %d x %d\n', size(y_all,1), size(y_all,2));
    % disp('y(end):');
    % disp(num2str(y_all(end,:), '%.12f '));
    % disp('y(end-1):');
    % disp(num2str(y_all(end-1,:), '%.12f '));
    % 
    % fprintf('  - JFNK_onestep():  %.4f s\n', time_JFNK_onestep);
    % fprintf('  - save_results():  %.4f s\n', time_save_results);

end
% [tplot_ssss, yplot_ssss, zplot_ssss] = JFNK_adaptive(@f_eval,t0,tfinal,h,p,y0,N_ITER_MAX_NEWTON,N_ITER_MAX_SDC,BE_TOL,SDC_TOL);
[tplot,yplot,zplot] = GMRES_SDC(n, m, h_GMRES, t0, y0, tfinal);
[tplot_ssss, yplot_ssss, zplot_ssss] = JFNK_adaptive(@f_eval,t0,tfinal,h,p,y0,N_ITER_MAX_NEWTON,N_ITER_MAX_SDC,BE_TOL,SDC_TOL);
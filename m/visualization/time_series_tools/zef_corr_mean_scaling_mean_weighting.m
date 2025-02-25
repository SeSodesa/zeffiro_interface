function [y_vals, plot_mode] = zef_corr_mean_scaling_mean_weighting(time_series)
%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
%This function processes the N-by-M data array f for N channels and M time
%steps. The other arguments can be controlled via the ZI user interface.
%The desctiption and argument definitions shown in ZI are listed below.
%Description: Correlation, mean scaling, mean weighting

time_series = time_series./max(time_series);
D = diag(sqrt(max(time_series,[],2)));
D = D./max(D(:));
y_vals = corr(time_series');
y_vals(find(isnan(y_vals))) = 1;
y_vals = D*y_vals*D;

plot_mode = 2;

function [y_vals, plot_mode] = zef_corr_max_scaling_max_weighting(time_series)
%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
%This function processes the N-by-M data array f for N channels and M time
%steps. The other arguments can be controlled via the ZI user interface.
%The desctiption and argument definitions shown in ZI are listed below.
%Description: Correlation, max scaling, max weighting

time_series = time_series./mean(time_series);
y_vals = cov(time_series');

plot_mode = 2;

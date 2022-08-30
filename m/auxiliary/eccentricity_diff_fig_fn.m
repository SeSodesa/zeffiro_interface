function [mag_fig, rdm_fig] = eccentricity_diff_fig_fn( ...
    source_points, ...
    mags, ...
    rdms, ...
    legend_labels, ...
    n_intervals, ...
    x_scale, ...
    min_ecc ...
)

    % eccentricity_diff_fig_fn
    %
    % A function that generates two figures comparing the MAG and RDM measures
    % between an analytical and numerical lead fields calculated based on a
    % set of source point at different eccentricities.
    %
    % Input:
    %
    % - source_points
    %
    %   A cell array of 3 × N matrix of source point xyz-coordinates. In other words, the
    %   number of columns in the matrix.
    %
    % - mags
    %
    %   The MAG measures calculated for the groups of source_points.
    %
    % - rdms
    %
    %   The RDM measures calculated for the groups of source_points.
    %
    % - n_intervals
    %
    %   The number of box plot groups one wishes to have displayed in the
    %   final figures.
    %
    % - x_scale
    %
    %   The scale of the x-axis.
    %
    % - min_ecc
    %
    %   The minimum eccentricity at which result display is attempted at.
    %
    % - legend_labels
    %
    %   The labels of the groups of points found in source_points. The size
    %   must match M, the number of source point groups in source_points.
    %

    arguments

        source_points cell { mustBeVector }

        mags cell { mustBeVector }

        rdms cell { mustBeVector }

        legend_labels cell { mustBeText, mustBeVector }

        n_intervals (1,1) double { mustBeInteger, mustBePositive } = 5;

        x_scale (1,1) double = 1.2;

        min_ecc (1,1) double { ...
            mustBeGreaterThanOrEqual(min_ecc, 0), ...
            mustBeLessThan(min_ecc, 1) ...
        } = 0;

    end

    % Set return values.

    rdm_fig = figure(100); set(rdm_fig, 'renderer','painters', 'position', [20,20,800,600]);

    rdm_ax = axes(rdm_fig);

    mag_fig = figure(101); set(mag_fig, 'renderer','painters', 'position', [20,20,800,600]);

    mag_ax = axes(mag_fig);

    % Then start filling the figures in.

    source_points_len = numel(source_points);

    % Iterate over the source matrices in source_points and generate index
    % vectors.

    ind_vec_cell = cell(source_points_len,1);

    for spi = 1 : source_points_len

        source_point_group = source_points{spi};

        group_size = size(source_point_group);

        if mod(group_size(1), 3) ~= 0

            close(rdm_fig);

            close(mag_fig);

            error("The number of rows in source point matrix was not divisible by 3.");

        end

        n_s_p = sqrt(sum(source_point_group.^2,1));
        n_s_p = repmat(n_s_p,3,1);
        n_s_p = n_s_p(:);
        n_s_p = n_s_p/max(n_s_p);

        % Populate index vector.

        index_vec = zeros(size(n_s_p));

        % Lay out tick labels.

        for ind = 1 : n_intervals

            I = find(n_s_p >= min_ecc + (ind-1)*(1-min_ecc)/n_intervals & n_s_p <= min_ecc + ind*(1-min_ecc)/n_intervals);

            index_vec(I) = ind;

            tick_labels{ind} = round(min_ecc + 0.5*(ind+(ind-1))*(1-min_ecc)/n_intervals,3);

        end % for

        ind_vec_cell{spi} = index_vec;

    end % for

    % Set tick label positions.

    tick_labels = cell(n_intervals,1);

    for ind = 1 : n_intervals

        tick_labels{ind} = round(min_ecc + 0.5*(ind+(ind-1))*(1-min_ecc)/n_intervals,3);

    end % for

    % Group indices stored in index vectors.

    flattened_ind_vecs = vertcat(ind_vec_cell{:});

    n_of_inds = numel(flattened_ind_vecs);

    group_vec = zeros(n_of_inds, 1);

    low = 1;

    for ind = 1 : source_points_len

        n_cols = 3 * size(source_points{ind},2);

        high = low + n_cols - 1;

        group_vec(low:high) = ind_vec_cell{ind};

        low = high + 1;

    end

    cat_group_vec = categorical(group_vec);

    % Add graphics to plots.

    rdm_vec = vertcat(rdms{:});

    h_rdm = boxchart( ...
        rdm_ax, ...
        x_scale * flattened_ind_vecs, ...
        log10(rdm_vec), ...
        'GroupByColor', cat_group_vec ...
    );

    for i = 1 : length(h_rdm)
        set(h_rdm(i),'MarkerStyle','.')
        set(h_rdm(i),'MarkerSize',0.5)
        set(h_rdm(i),'JitterOutliers','off')
        set(h_rdm(i),'linewidth',0.25)
    end

    set(rdm_ax,'Xtick',x_scale*[1:n_intervals])
    set(rdm_ax,'XtickLabels',tick_labels)
    set(rdm_ax,'xlim',[x_scale/2 x_scale*(n_intervals+1/2)])
    pbaspect([3 1 1])
    set(rdm_ax,'fontsize',12)
    set(rdm_ax,'linewidth',0.25)
    set(rdm_ax,'xgrid','on')
    set(rdm_ax,'ygrid','on')
    set(rdm_ax,'box','on')
    set(rdm_ax,'ylim',[-3 0])
    set(rdm_ax,'ytick',log10([0.00001 0.00002 0.00004 0.00006 0.00008 0.0001 0.0002 0.0004 0.0006 0.0008 0.001 0.002 0.004 0.006 0.008 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.4 0.6 0.8 1 2 4 6 8 10]))
    set(rdm_ax,'yticklabels',{'1E-5','2E-5','4E-5','6E-5','8E-5','1E-4','2E-4','4E-4','6E-4','8E-4','1E-3','2E-3','4E-3','6E-3','8E-3','1E-2','2E-2','4E-2','6E-2','8E-2','1E-1','2E-1','4E-1','6E-1','8E-1','1','2','4','6','8','10'})
    set(rdm_ax,'ticklength',[0.0100 0.0250]/2)

    legend(legend_labels,'Orientation','Horizontal','Location','Northwest');

    mag_vec = vertcat(mags{:});

    h_mag = boxchart( ...
        mag_ax, ...
        x_scale * flattened_ind_vecs, ...
        log10(abs(mag_vec)), ...
        'GroupByColor',group_vec ...
    );

    for i = 1 : length(h_mag)
        set(h_mag(i),'MarkerStyle','.')
        set(h_mag(i),'MarkerSize',0.5)
        set(h_mag(i),'JitterOutliers','off')
        set(h_mag(i),'linewidth',0.25)
    end

    set(mag_ax,'Xtick',x_scale*[1:n_intervals])
    set(mag_ax,'XtickLabels',tick_labels)
    set(mag_ax,'xlim',[0 x_scale*(n_intervals+1)])

    pbaspect([3 1 1])

    set(mag_ax,'fontsize',12)
    set(mag_ax,'linewidth',0.25)
    set(mag_ax,'xgrid','on')
    set(mag_ax,'ygrid','on')
    set(mag_ax,'box','on')
    set(mag_ax,'ylim',[-3 0.5])
    set(mag_ax,'ytick',log10([0.00001 0.00002 0.00004 0.00006 0.00008 0.0001 0.0002 0.0004 0.0006 0.0008 0.001 0.002 0.004 0.006 0.008 0.01 0.02 0.04 0.06 0.08 0.1 0.2 0.4 0.6 0.8 1 2 4 6 8 10]))
    set(mag_ax,'yticklabels',{'1E-5','2E-5','4E-5','6E-5','8E-5','1E-4','2E-4','4E-4','6E-4','8E-4','1E-3','2E-3','4E-3','6E-3','8E-3','1E-2','2E-2','4E-2','6E-2','8E-2','1E-1','2E-1','4E-1','6E-1','8E-1','1','2','4','6','8','10'})
    set(mag_ax,'ticklength',[0.0100 0.0250]/2)

    legend(legend_labels,'Orientation','Horizontal','Location','Northwest')

end % function

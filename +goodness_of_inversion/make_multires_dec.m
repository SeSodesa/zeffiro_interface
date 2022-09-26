function [multires_dec, multires_ind, multires_count] = make_multires_dec( ...
    zef, ...
    n_decompositions, ...
    n_levels, ...
    multires_sparsity ...
)

    % make_multigrid_dec
    %
    % Copied over from zef_make_multires_dec and modified to accept zef as an
    % argument.

    s_ind = unique(zef.source_interpolation_ind{1});

    for n_rep = 1 : n_decompositions

        source_points = zef.source_positions;

        source_points = source_points(s_ind,:);
        source_points = source_points';
        source_points_aux = source_points;
        size_center_points = size(source_points,2);
        center_points = source_points;

        multires_dec{n_rep}{n_levels} = [1:size_center_points]';
        multires_ind{n_rep}{n_levels} = [1:size_center_points]';
        multires_count{n_rep}{n_levels} = ones(size_center_points,1);

        par_num = zef.parallel_vectors;
        bar_ind = ceil(size_center_points/(50*par_num));

        use_gpu  = zef.use_gpu;
        gpu_num  = zef.gpu_num;

        for k = 1 : n_levels-1

            size_source_points = floor(size_center_points/multires_sparsity^(n_levels - k));
            source_interpolation_aux = zeros(size_source_points,1);

            aux_ind = randperm(size_center_points);
            %aux_ind = sort(aux_ind);
            source_points = source_points_aux(:,aux_ind(1:size_source_points));
            ones_vec = ones(size(source_points,2),1);

            multires_dec{n_rep}{k} = aux_ind(1:size_source_points)';

            MdlKDT = KDTreeSearcher(source_points');
            source_interpolation_aux = knnsearch(MdlKDT,center_points');

            multires_ind{n_rep}{k} = source_interpolation_aux;
            [aux_vec, i_a, i_c] = unique(source_interpolation_aux);
            multires_count{n_rep}{k} = accumarray(i_c,1);

        end % for

        multires_dec{n_rep}{n_levels} = [1:size_center_points]';
        multires_ind{n_rep}{n_levels} = [1:size_center_points]';
        multires_count{n_rep}{n_levels} = ones(size_center_points,1);

    end % for

end % function

%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
function [z,reconstruction_information] = ramus( ...
    zef, ...
    params ...
)

    arguments

        zef (1,1) struct

        params.multiresolution_levels (1,1) double { mustBeInteger, mustBePositive } = 2;

        params.multiresolution_sparsity (1,1) double { mustBeReal, mustBePositive } = 10;

        params.n_of_decompositions (1,1) double { mustBeInteger, mustBePositive } = 10;

        params.hyperprior (1,1) string { mustBeMember( ...
            params.hyperprior, ...
            ["spatially balanced", "spatially constant"]) ...
        } = "spatially balanced";

        params.signal_to_noise_ratio (1,1) double { mustBeReal } = 30;

        params.ias_map_iterations (1,1) double { mustBeInteger, mustBePositive } = 10;

        params.sampling_frequency (1,1) double { mustBePositive } = 1025;

        params.low_cut_frequency (1,1) double { mustBePositive } = 7;

        params.high_cut_frequency (1,1) double { mustBePositive } = 9;

        params.time_start (1,1) double { mustBeNonnegative } = 0;

        params.time_window (1,1) double { mustBeNonnegative } = 0;

        params.n_of_time_steps (1,1) double { mustBeInteger, mustBePositive } = 1;

        params.time_step (1,1) double { mustBeInteger, mustBeNonnegative } = 0;

        params.data_normalization (1,1) string { mustBeMember( ...
            params.data_normalization, ...
            [ ...
                "maximum entry", ...
                "maximum column norm", ...
                "average column norm", ...
                "none" ...
            ] ...
        ) } = "maximum entry";

        params.initial_guess_mode (1,1) string { mustBeMember( ...
            params.initial_guess_mode, ...
            [ ...
                "hyperparameter value", ...
                "hyperprior expectation", ...
            ] ...
        ) } = "hyperparameter value";

        params.n_decompositions (1,1) double { mustBeInteger, mustBePositive } = 10;

        params.n_levels (1,1) double { mustBeInteger, mustBePositive } = 2;

        params.multires_sparsity (1,1) double { mustBeInteger, mustBePositive } = 10;

        params.hyperprior_tail_length_db (1,1) double { mustBeReal, mustBePositive } = 1

        params.hyperprior_weight (1,1) double { mustBeReal, mustBePositive } = 1 / 2

    end

    h = zef_waitbar([0 0 0],['RAMUS iteration.']);

    cleanup_fn = @(go) close(go);

    cleanup_obj = onCleanup(@() cleanup_fn(h));

    [s_ind_1] = unique(zef.source_interpolation_ind{1});
    n_interp = length(s_ind_1);
    n_multires = params.multiresolution_levels;
    ramus_hyperprior = params.hyperprior;
    sparsity_factor = params.multires_sparsity;
    snr_val = params.signal_to_noise_ratio;
    pm_val = zef.inv_prior_over_measurement_db;
    amplitude_db = zef.inv_amplitude_db;
    pm_val = pm_val - amplitude_db;
    sampling_freq = params.sampling_frequency;
    high_pass = params.low_cut_frequency;
    low_pass = params.high_cut_frequency;
    number_of_frames = params.n_of_time_steps;
    time_step = params.time_step;
    source_direction_mode = zef.source_direction_mode;
    source_directions = zef.source_directions;
    n_decompositions = params.n_of_decompositions;
    weight_vec_aux = (sparsity_factor.^[0:n_multires-1]');

    std_lhood = 10^(-snr_val/20);

    reconstruction_information.tag = 'RAMUS';
    reconstruction_information.inv_time_1 = params.time_start;
    reconstruction_information.inv_time_2 = params.time_window;
    reconstruction_information.inv_time_3 = params.time_step;
    reconstruction_information.sampling_freq = params.sampling_frequency;
    reconstruction_information.low_pass = params.high_cut_frequency;
    reconstruction_information.high_pass = params.low_cut_frequency;
    reconstruction_information.number_of_frames = params.n_of_time_steps;
    reconstruction_information.source_direction_mode = zef.source_direction_mode;
    reconstruction_information.source_directions = zef.source_directions;
    reconstruction_information.ias_hyperprior = params.hyperprior;
    reconstruction_information.snr_val = params.signal_to_noise_ratio;
    reconstruction_information.pm_val = zef.inv_prior_over_measurement_db;

    [L,n_interp, procFile] = goodness_of_inversion.process_lead_fields(zef);

    if zef.use_gpu == 1 && zef.gpu_count > 0
        L = gpuArray(L);
    end

    L_aux = L;

    S_mat = std_lhood^2*eye(size(L,1));

    if zef.use_gpu == 1 && zef.gpu_count > 0
        S_mat = gpuArray(S_mat);
    end

    [f_data] = goodness_of_inversion.get_filtered_data(zef);

    tic;

    z_inverse = cell(0);

    tic;

    for f_ind = 1 : number_of_frames

        time_val = toc;

        if f_ind > 1
            date_str = datestr(datevec(now+(number_of_frames/(f_ind-1) - 1)*time_val/86400));
        end

        if source_direction_mode == 1 || source_direction_mode == 2
            z_aux = zeros(size(L_aux,2),1);
        end

        if source_direction_mode == 3
            z_aux = zeros(3*size(L_aux,2),1);
        end

        z_vec = ones(size(L_aux,2),1);

        [f] = goodness_of_inversion.get_time_step( ...
            zef, ...
            f_data, ...
            f_ind, ...
            params.time_start, ...
            params.time_window, ...
            params.time_step, ...
            params.sampling_frequency, ...
            'Optional_averaging_bool', true ...
        );

        if f_ind == 1
            zef_waitbar([0 0 0],h,['IAS MAP iteration. Time step ' int2str(f_ind) ' of ' int2str(number_of_frames) '.']);
        end

        n_ias_map_iter = params.ias_map_iterations;

        if zef.use_gpu == 1 && zef.gpu_count > 0
            f = gpuArray(f);
        end

        multires_dec =  zef.ramus_multires_dec;
        multires_ind =  zef.ramus_multires_ind;
        multires_count = zef.ramus_multires_count;
        n_iter = zef.ramus_multires_n_iter;

        if length(n_iter) < n_multires
            n_iter = n_iter(1)*ones(1,n_multires);
        end

        mr_sparsity = zef.ramus_multires_sparsity;

        z_vec_aux = zeros(size(L_aux,2),1);
        iter_ind = 0;
        source_count_aux = 0;

        for n_rep = 1 : n_decompositions

            for j = 1 : n_multires

                iter_ind = iter_ind + 1;

                n_mr_dec = length(multires_dec{n_rep}{j});

                if source_direction_mode == 1 || source_direction_mode == 2
                    mr_dec = [multires_dec{n_rep}{j}; multires_dec{n_rep}{j}+n_interp ; multires_dec{n_rep}{j} + 2*n_interp];
                    mr_dec = mr_dec(:);
                    mr_ind = [multires_ind{n_rep}{j} ; multires_ind{n_rep}{j} + n_mr_dec ; multires_ind{n_rep}{j} + 2*n_mr_dec];
                    mr_ind = mr_ind(:);
                end

                if source_direction_mode == 3
                    mr_dec = multires_dec{n_rep}{j};
                    mr_dec = mr_dec(:);
                    mr_ind = multires_ind{n_rep}{j};
                    mr_ind = mr_ind(:);
                end

                if n_iter(j) > 0

                    L_aux_2 = L_aux(:,mr_dec);

                    if source_count_aux == 0
                        source_count = size(L_aux_2,2);
                        source_count_aux = 1;
                    end

                    if params.data_normalization == "maximum entry"
                        normalize_data = "maximum";
                    else
                        normalize_data = "L2";
                    end

                    if ramus_hyperprior == "spatially balanced"
                        balance_spatially = true;
                    else
                        balance_spatially = false;
                    end

                    if params.hyperprior == "spatially balanced"

                        [beta, theta0] = goodness_of_inversion.find_ig_hyperprior( ...
                            L_aux_2, ...
                            snr_val-pm_val, ...
                            params.hyperprior_tail_length_db, ...
                            source_count, ...
                            normalize_data, ...
                            balance_spatially, ...
                            params.hyperprior_weight ...
                        );

                    elseif params.hyperprior == "spatially constant"

                        [beta, theta0] = zef_find_g_hyperprior( ...
                            L_aux_2, ...
                            snr_val-pm_val, ...
                            params.hyperprior_tail_length_db, ...
                            source_count, ...
                            params.data_normalization, ...
                            balance_spatially, ...
                            params.hyperprior_weight ...
                        );

                    end

                    if n_rep == 1 || zef.ramus_init_guess_mode == 2

                        if zef.inv_hyperprior == 1

                            if length(theta0) > 1  || length(beta) > 1
                                theta = theta0./(beta-1);
                            else
                                theta = (theta0./(beta-1))*ones(size(L_aux_2,2),1);
                            end

                        elseif zef.inv_hyperprior == 2

                            if length(theta0) > 1  || length(beta) > 1
                                theta = theta0.*beta;
                            else
                                theta = (theta0.*beta)*ones(size(L_aux_2,2),1);
                            end

                        end

                    else

                        theta = theta(mr_dec);

                    end % if

                    for i = 1 : n_iter(j)

                        if f_ind > 1
                            zef_waitbar([i/n_iter(j) j/n_multires n_rep/n_decompositions],h,['Dec. ' int2str(n_rep) ' of ' int2str(n_decompositions) ', Step ' int2str(f_ind) ' of ' int2str(number_of_frames) '. Ready: ' date_str '.' ]);
                        else
                            zef_waitbar([i/n_iter(j) j/n_multires n_rep/n_decompositions],h,['IAS MAP iteration. Dec. ' int2str(n_rep) ' of ' int2str(n_decompositions) ', Time step ' int2str(f_ind) ' of ' int2str(number_of_frames) '.' ]);
                        end

                        d_sqrt = sqrt(theta);

                        if zef.use_gpu == 1 && zef.gpu_count > 0
                            d_sqrt = gpuArray(d_sqrt);
                        end

                        L = L_aux_2.*repmat(d_sqrt',size(L,1),1);

                        z_vec = d_sqrt.*(L'*((L*L' + S_mat)\f));

                        if zef.use_gpu == 1 && zef.gpu_count > 0
                            z_vec = gather(z_vec);
                        end

                        if zef.inv_hyperprior == 1
                            theta = (theta0+0.5*z_vec.^2)./(beta + 1.5);
                        elseif zef.inv_hyperprior == 2
                            theta = theta0.*(beta-1.5 + sqrt((1./(2.*theta0)).*z_vec.^2 + (beta+1.5).^2));
                        end

                    end % for

                    if length(theta0) > 1
                        theta0 = theta0(mr_ind);
                    end

                    if length(beta) > 1
                        beta = beta(mr_ind);
                    end

                    theta = theta(mr_ind);

                    z_vec = z_vec(mr_ind);

                else

                    z_vec = zeros(length(mr_ind),1);

                    weight_vec_aux(j) = 0;

                end % if

                z_vec_aux = z_vec_aux + z_vec;

            end % for

        end % for

        z_vec = z_vec_aux/(n_multires*n_decompositions*sum(weight_vec_aux));

        z_inverse{f_ind} = z_vec;

    end % for

    [z] = goodness_of_inversion.postprocess_reconstruction(z_inverse, procFile);
    [z] = goodness_of_inversion.normalize_reconstruction(z);

end % function

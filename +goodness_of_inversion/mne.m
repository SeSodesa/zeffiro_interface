function [rec, info] = mne(zef, mne_type, params)

    arguments

        zef struct

        mne_type (1,1) string { mustBeMember(mne_type, ["MNE", "sLORETA", "dSPM"]) }

        params.low_cut_frequency (1,1) double { mustBePositive } = 7;

        params.high_cut_frequency (1,1) double { mustBePositive } = 9;

        params.normalize_data (1,1) string { mustBeMember( ...
            params.normalize_data, ...
            [ "maximum entry", "maximum column norm", "average column norm", "none" ] ...
        ) } = "maximum entry";

        params.number_of_frames (1,1) double { mustBeInteger, mustBePositive } = 1;

        params.prior (1,1) string { mustBeMember( ...
            params.prior, ...
            [ "balanced", "constant" ] ...
        ) } = "balanced";

        params.sampling_frequency (1,1) double { mustBeReal, mustBePositive } = 1025;

        params.time_start (1,1) double { mustBeReal, mustBeNonnegative } = 0;

        params.time_window (1,1) double { mustBeReal, mustBeNonnegative } = 0;

        params.time_step (1,1) double { mustBeReal, mustBePositive } = 1;

        params.signal_to_noise_ratio (1,1) double { mustBeReal, mustBePositive } = 30;

        params.inv_amplitude_db (1,1) double = 20;

        params.inv_prior_over_measurement_db (1,1) double = 20;

    end

    inverse_gamma_ind = [1:4];

    gamma_ind = [5:10];

    h = zef_waitbar(0,['MNE Reconstruction.']);

    cleanup_fn = @(graphics) close(graphics);

    cleanup_obj = onCleanup(@() cleanup_fn(h));

    [procFile.s_ind_1] = zef.source_interpolation_ind{1};

    n_interp = length(procFile.s_ind_1);

    pm_val = params.inv_prior_over_measurement_db;
    amplitude_db = params.inv_amplitude_db;
    pm_val = pm_val - amplitude_db;

    snr_val = params.signal_to_noise_ratio;
    mne_type = mne_type;
    mne_prior = params.prior;
    std_lhood = 10^(-snr_val/20);

    zef.inv_sampling_frequency = params.sampling_frequency;
    zef.inv_high_pass = params.low_cut_frequency;
    zef.inv_low_pass = params.high_cut_frequency;
    zef.number_of_frames = params.number_of_frames;
    zef.inv_time_1 = params.time_start;
    zef.inv_time_2 = params.time_window;
    zef.inv_time_3 = params.time_step;

    source_direction_mode = zef.source_direction_mode;

    source_directions = zef.source_directions;

    info=[];
    info.tag='mne_type';
    info.type='mne_type';
    info.std_lhood=std_lhood;

    info.snr_val = params.signal_to_noise_ratio;
    info.mne_type = mne_type;
    info.mne_prior = params.prior;
    info.sampling_freq = params.sampling_frequency;
    info.high_pass = params.low_cut_frequency;
    info.low_pass = params.high_cut_frequency;
    info.number_of_frames = zef.mne_number_of_frames;
    info.time_step = params.time_step;
    info.source_direction_mode = zef.source_direction_mode;
    info.source_directions = zef.source_directions;

    [L,n_interp, procFile] = zef_processLeadfields(zef);

    source_count = n_interp;

    if params.prior == "balanced"

        balance_spatially = 1;

    else

        balance_spatially = 0;

    end

    [theta0] = goodness_of_inversion.find_gaussian_prior( ...
        snr_val-pm_val, ...
        L, ...
        "source_space_size", size(L,2), ...
        "normalization_type", params.normalize_data, ...
        "balance_snr", balance_spatially ...
    );

    if zef.use_gpu == 1 & zef.gpu_count > 0

        L = gpuArray(L);

    end

    S_mat = std_lhood^2*eye(size(L,1));

    if zef.use_gpu == 1 & zef.gpu_count > 0

        S_mat = gpuArray(S_mat);

    end

    if params.number_of_frames > 1

        rec = cell(params.number_of_frames,1);

    else

        zef.number_of_frames = 1;

    end


    [f_data] = goodness_of_inversion.get_filtered_data( ...
        zef, ...
        false, ...
        "low_cut_frequency", params.low_cut_frequency, ...
        "high_cut_frequency", params.high_cut_frequency ...
    );

    tic;

    for f_ind = 1 : zef.number_of_frames

        time_val = toc;

        if f_ind > 1

            date_str = datestr(datevec(now+(zef.number_of_frames/(f_ind-1) - 1)*time_val/86400));

        end

        if ismember(source_direction_mode, [1,2])

            z_aux = zeros(size(L,2),1);

        end

        if source_direction_mode == 3

            z_aux = zeros(3*size(L,2),1);

        end

        z_vec = ones(size(L,2),1);

        [f] = zef_getTimeStep(f_data, f_ind, zef);

        if f_ind == 1

            zef_waitbar(0,h,['MNE reconstruction. Time step ' int2str(f_ind) ' of ' int2str(zef.number_of_frames) '.']);

        end

        if zef.use_gpu == 1 & zef.gpu_count > 0

            f = gpuArray(f);

        end

        if f_ind > 1

            zef_waitbar(f_ind/zef.number_of_frames,h,['Step ' int2str(f_ind) ' of ' int2str(zef.number_of_frames) '. Ready: ' date_str '.' ]);

        else

            zef_waitbar(f_ind/zef.number_of_frames,h,['MNE reconstruction. Time step ' int2str(f_ind) ' of ' int2str(zef.number_of_frames) '.' ]);

        end

        m_max = sqrt(size(L,2));
        u = zeros(length(z_vec),1);
        z_vec = zeros(length(z_vec),1);

        if length(theta0)==1

            d_sqrt = sqrt(theta0)*ones(size(z_vec));

        else

            d_sqrt = sqrt(theta0);

        end

        if zef.use_gpu == 1 & zef.gpu_count > 0

            d_sqrt = gpuArray(d_sqrt);

        end

        L_inv = L.*repmat(d_sqrt',size(L,1),1);
        L_inv = d_sqrt.*(L_inv'*(inv(L_inv*L_inv' + S_mat)));

        if mne_type == "dSPM"

            aux_vec = sum(L_inv.^2, 2);
            aux_vec = sqrt(aux_vec);
            L_inv = L_inv./aux_vec(:,ones(size(L_inv,2),1));

        elseif mne_type == "sLORETA"

            aux_vec = sqrt(sum(L_inv.*L', 2));
            L_inv = L_inv./aux_vec(:,ones(size(L_inv,2),1));

        end

        z_vec = L_inv * f;

        if zef.use_gpu == 1 & zef.gpu_count > 0

            z_vec = gather(z_vec);

        end

        z_inverse{f_ind} = z_vec;

    end % for

    [rec] = goodness_of_inversion.postprocess_reconstruction(z_inverse, procFile);

end % function

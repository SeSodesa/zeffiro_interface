function zef = init_inverse_method(zef, inverse_method_name, invp, mne, beamformer, ramus, sesame)

    % init_inverse_method
    %
    % Set required values in zef based on the given inverse method name and
    % related given arguments.
    %
    % Inputs:
    %
    % - zef
    %
    %   The zef struct that will be passed to the inverse methods and from
    %   which the necessary fields will be read.
    %
    % - inverse_method_name
    %
    %   This determines which zef values will be set by this function.
    %
    % - mne
    %
    %   The MNE-related key–value pairs given to the main function, stored in
    %   a struct.
    %
    % - beamformer
    %
    %   The Beamformer-related key–value pairs given to the main function, stored in
    %   a struct.
    %
    % - ramus
    %
    %   The RAMUS-related key–value pairs given to the main function, stored in
    %   a struct.
    %
    % - TODO sesame
    %
    %   The SESAME-related key–value pairs given to the main function, stored
    %   in a struct.
    %
    % Outputs:
    %
    % - zef
    %
    %   The input zef struct with the needed inverse method fields set.
    %

    arguments

        zef struct

        inverse_method_name (1,1) string { mustBeMember(inverse_method_name, ["MNE", "Beamformer", "RAMUS", "SESAME"]) }

        invp struct

        mne struct

        beamformer struct

        ramus struct

        sesame struct

    end

    if strcmp(inverse_method_name, "MNE")

        zef.mne_low_cut_frequency = mne.mne_low_cut_frequency;

        zef.mne_high_cut_frequency = mne.mne_high_cut_frequency;

        zef.mne_normalize_data = mne.mne_normalize_data;

        zef.mne_number_of_frames = mne.mne_number_of_frames;

        zef.mne_prior = mne.mne_prior;

        zef.mne_sampling_frequency = mne.mne_sampling_frequency;

        zef.mne_time_start = mne.mne_time_start;

        zef.mne_time_window = mne.mne_time_window;

        zef.mne_time_step = mne.mne_time_step;

    end

    if strcmp(inverse_method_name, "Beamformer")

        zef.beamformer_estimate_type = beamformer.beamformer_estimate_type;

        zef.beamformer_estimated_attribute = beamformer.beamformer_estimated_attribute;

        zef.beamformer_covariance_regularization_parameter = beamformer.beamformer_covariance_regularization_parameter;

        zef.beamformer_covariance_mode = beamformer.beamformer_covariance_mode;

        zef.beamformer_lead_field_regularization_parameter = beamformer.beamformer_lead_field_regularization_parameter;

        zef.beamformer_lead_field_regularization_procedure = beamformer.beamformer_lead_field_regularization_procedure;

        zef.beamformer_signal_to_noise_ratio = beamformer.beamformer_signal_to_noise_ratio;

        zef.beamformer_sampling_frequency = beamformer.beamformer_sampling_frequency;

        zef.beamformer_low_cut_frequency = beamformer.beamformer_low_cut_frequency;

        zef.beamformer_high_cut_frequency = beamformer.beamformer_high_cut_frequency;

        zef.beamformer_time_start = beamformer.beamformer_time_start;

        zef.beamformer_time_window = beamformer.beamformer_time_window;

        zef.beamformer_n_of_time_steps = beamformer.beamformer_n_of_time_steps;

        zef.beamformer_time_step = beamformer.beamformer_time_step;

        zef.beamformer_data_segment = beamformer.beamformer_data_segment;

        zef.beamformer_data_normalization = beamformer.beamformer_data_normalization;

        zef.beamformer_lead_field_normalization = beamformer.beamformer_lead_field_normalization;

        zef.beamformer_inv_cov_lambda = beamformer.beamformer_inv_cov_lambda;

        zef.beamformer_inv_lead_field_lambda = beamformer.beamformer_inv_lead_field_lambda;

    end

    if strcmp(inverse_method_name, "RAMUS")

        zef.ramus_multiresolution_levels = ramus.ramus_multiresolution_levels;

        zef.ramus_multiresolution_sparsity = ramus.ramus_multiresolution_sparsity;

        zef.ramus_n_of_decompositions = ramus.ramus_n_of_decompositions;

        zef.ramus_hyperprior = ramus.ramus_hyperprior;

        zef.ramus_signal_to_noise_ratio = ramus.ramus_signal_to_noise_ratio;

        zef.ramus_ias_map_iterations = ramus.ramus_ias_map_iterations;

        zef.ramus_sampling_frequency = ramus.ramus_sampling_frequency;

        zef.ramus_low_cut_frequency = ramus.ramus_low_cut_frequency;

        zef.ramus_high_cut_frequency = ramus.ramus_high_cut_frequency;

        zef.ramus_time_start = ramus.ramus_time_start;

        zef.ramus_time_window = ramus.ramus_time_window;

        zef.ramus_n_of_time_steps = ramus.ramus_n_of_time_steps;

        zef.ramus_time_step = ramus.ramus_time_step;

        zef.ramus_data_normalization = ramus.ramus_data_normalization;

        zef.ramus_initial_guess_mode = ramus.ramus_initial_guess_mode;

        [decomposition, indices, count] = goodness_of_inversion.make_multires_dec( ...
            zef, ...
            ramus.ramus_n_decompositions, ...
            ramus.ramus_n_levels, ...
            ramus.ramus_multires_sparsity ...
        );

        zef.ramus_multires_dec = decomposition;

        zef.ramus_multires_ind = indices;

        zef.ramus_multires_count = count;

    end

    if strcmp(inverse_method_name, "SESAME")

        error("SESAME is still to be implemented. Aborting.")

    end
end

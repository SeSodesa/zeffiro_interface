function zef = init_inverse_method(zef, inverse_method_name, mne, beamformer, ramus, sesame)

    % init_inverse_method
    %
    % Set required values in zef based on the given inverse method name and
    % related given arguments.
    %

    arguments

        zef struct

        inverse_method_name (1,1) string { mustBeMember(inverse_method_name, ["MNE", "Beamformer", "RAMUS", "SESAME"]) }

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

    end

    if strcmp(inverse_method_name, "RAMUS")
        % zef_init_ramus_inversion_tool;
        % zef_update_ramus_inversion_tool;
    end

    if strcmp(inverse_method_name, "SESAME")
        % SESAME_App_run;
        % zef_SESAME_init;
        % zef_update_SESAME;
    end
end

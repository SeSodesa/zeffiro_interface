function [zef, rec_vec_position, rec_vec_angle, rec_vec_magnitude] = main( ...
    zef,           ...
    n_subset,      ...
    roi_radius,    ...
    inverse_method_name, ...
    mne, ...
    beamformer, ...
    ramus, ...
    sesame ...
)

% goodness_of_inversion.main
%
% Calculates a reconstruction of a dipole source from a lead field L inside a
% Region of Interest (RoI) around the source dipole. Returns the
% reconstruction dipole positions, angles and magnitudes.
%
%
% NOTE: The return values are empty, if an error occurred.
%
% Input:
%
% - zef
%
%   An instance of zef for reading values like model source positions.
%
% - n_subset
%
%   The number of subdivisions of the region of interest.
%
% - roi_radius
%
%   The radius of the region of interest.
%
% - inverse_method_name
%
%   A string that indicates the inverse method used to calculate the
%   reconstruction from the lead field. Currently supported inversion methods
%   are
%
%   - "mne"
%   - "beamformer"
%   - "ramus"
%   - "sesame"
%
%   In order to work, the following lists of fields need to exist within zef
%   for the inverse methods to function at all:
%
%   - MNE
%
%     - compartment_tags
%     - gpu_count
%     - inv_prior_over_measurement_db
%     - inv_amplitude_db
%     - inv_snr
%     - L
%     - measurements
%     - mne_data_segment
%     - mne_high_cut_frequency
%     - mne_low_cut_frequency
%     - mne_normalize_data
%     - mne_number_of_frames
%     - mne_prior
%     - mne_sampling_frequency
%     - mne_time_1
%     - mne_time_2
%     - mne_time_3
%     - mne_type
%     - reuna_t
%     - reuna_p
%     - source_direction_mode
%     - source_directions
%     - source_interpolation_ind
%     - use_gpu
%
%   - Beamformer
%
%     - beamformer.normalize_leadfield.Value
%     - bf_type
%     - cov_type
%     - inv_cov_lambda
%     - inv_high_cut_frequency
%     - inv_leadfield_lambda
%     - inv_low_cut_frequency
%     - inv_sampling_frequency
%     - inv_snr
%     - inv_time_1
%     - inv_time_2
%     - inv_time_3
%     - L_reg_type
%     - number_of_frames
%     - source_direction_mode
%     - source_directions
%     - source_interpolation_ind
%
% - mne
%
%   A set of key–value argument pairs for setting MNE parameters. These are
%   all prefixed with "mne_". These are also used by sLORETA, which is just
%   another version of MNE. The options and their default values are as
%   follows:
%
%   Key                         Value
%   -------------------------   ----------------------------------------
%
%   mne_low_cut_frequency       A positive real number, default = 7.
%
%   mne_high_cut_frequency      A positive real number, default = 9.
%
%   mne_normalize_data          "maximum entry" (default), "maximum column norm",
%                               "average column norm" or "none".
%
%   mne_number_of_frames        A non-negative whole number, default = 1.
%
%   mne_prior                   "balanced" (default) or "constant"
%
%   mne_sampling_frequency      A positive real number, default = 1025.
%
%   mne_time_start              A positive real number, default = 0.
%
%   mne_time_window             A positive real number, default = 0
%
%   mne_time_step               A positive real number, default = 1
%
% - beamformer
%
%   A set of key–value argument pairs for setting Beamformer parameters. These
%   are all prefixed with "beamformer_".
%
% - ramus
%
%   A set of key–value argument pairs for setting RAMUS parameters. These are
%   all prefixed with "ramus_".
%
% - sesame
%
%   A set of key–value argument pairs for setting SESAME parameters. These are
%   all prefixed with "sesame_".
%
% Output:
%
% - rec_vec_position
%
%   The dipole positions in the reconstruction.
%
% - rec_vec_angle
%
%   The angles of the dipoles in the reconstruction.
%
% - rec_vec_mag
%
%   The magnitudes of the dipole moments in the reconstruction.
%

    arguments

        zef struct

        n_subset (1,1) double { mustBeInteger, mustBePositive }

        roi_radius (1,1) double { mustBeReal, mustBePositive }

        inverse_method_name (1,1) string { mustBeMember(inverse_method_name, ["MNE", "sLORETA", "Beamformer", "RAMUS", "SESAME"]) }

        % MNE parameters.

        mne.mne_low_cut_frequency (1,1) double { mustBePositive } = 7;

        mne.mne_high_cut_frequency (1,1) double { mustBePositive } = 9;

        mne.mne_normalize_data (1,1) string { mustBeMember( ...
                mne.mne_normalize_data, ...
                [ "maximum entry", "maximum column norm", "average column norm", "none" ] ...
            ) } = "maximum entry";

        mne.mne_number_of_frames (1,1) double { mustBeInteger, mustBePositive } = 1;

        mne.mne_prior (1,1) string { mustBeMember( ...
                mne.mne_prior, ...
                [ "balanced", "constant" ] ...
            ) } = "balanced";

        mne.mne_sampling_frequency (1,1) double { mustBeReal, mustBePositive } = 1025.

        mne.mne_time_start (1,1) double { mustBeReal, mustBeNonnegative } = 0.

        mne.mne_time_window (1,1) double { mustBeReal, mustBeNonnegative } = 0

        mne.mne_time_step (1,1) double { mustBeReal, mustBePositive } = 1

        % TODO Beamformer parameters.

        beamformer.beamformer_estimate_type (1,1) string { ...
            mustBeMember( ...
                beamformer.beamformer_estimate_type, ...
                [ ...
                    "linearly constrained minimum variance", ...
                    "unit noise gain", ...
                    "unit-gain-constrained", ...
                    "unit noise gain scalar" ...
                ] ...
            ) ...
        } = "linearly constrained minimum variance";

        beamformer.beamformer_estimated_attribute (1,1) string { ...
            mustBeMember( ...
                beamformer.beamformer_estimated_attribute, ...
                ["dipole moments", "locations", "both"] ...
            ) ...
        } = "dipole moments";

        beamformer.beamformer_covariance_regularization_parameter (1,1) double = 0.05;

        beamformer.beamformer_covariance_mode (1,1) string { ...
            mustBeMember( ...
                beamformer.beamformer_covariance_mode, ...
                [ ...
                    "full data, measurement based", ...
                    "full data, basic", ...
                    "pointwise, measurement based", ...
                    "pointwise, basic" ...
                ] ...
            ) ...
        }= "full data, measurement based";

        beamformer.beamformer_lead_field_regularization_parameter (1,1) double = 0.001;

        beamformer.beamformer_lead_field_regularization_procedure (1,1) string { ...
            mustBeMember( ...
                beamformer.beamformer_lead_field_regularization_procedure, ...
                ["pseudoinverse","basic"] ...
            ) ...
        } = "basic";

        beamformer.beamformer_signal_to_noise_ratio (1,1) double = 30;

        beamformer.beamformer_sampling_frequency (1,1) double { mustBePositive } = 1025;

        beamformer.beamformer_low_cut_frequency (1,1) double { mustBePositive } = 7;

        beamformer.beamformer_high_cut_frequency (1,1) double { mustBePositive } = 9;

        beamformer.beamformer_time_start (1,1) double { mustBeNonnegative } = 0;

        beamformer.beamformer_time_window (1,1) double { mustBeNonnegative } = 0;

        beamformer.beamformer_n_of_time_steps (1,1) double { mustBeInteger, mustBePositive } = 1;

        beamformer.beamformer_time_step (1,1) double { mustBeInteger, mustBeNonnegative } = 0;

        beamformer.beamformer_data_segment (1,1) double { mustBeInteger, mustBePositive } = 1;

        beamformer.beamformer_data_normalization (1,1) string { mustBeMember( ...
            beamformer.beamformer_data_normalization, ...
            [ ...
                "maximum entry", ...
                "maximum column norm", ...
                "average column norm", ...
                "none" ...
            ] ...
        ) } = "maximum entry";

        beamformer.beamformer_lead_field_normalization (1,1) string { mustBeMember( ...
            beamformer.beamformer_lead_field_normalization, ...
            [ ...
                "matrix norm", ...
                "column norm", ...
                "row norm", ...
                "none" ...
            ] ...
        ) } = "matrix norm";

        % TODO RAMUS parameters.

        ramus.ramus_multiresolution_levels (1,1) double { mustBeInteger, mustBePositive } = 2;

        ramus.ramus_multiresolution_sparsity (1,1) double { mustBeReal, mustBePositive } = 10;

        ramus.ramus_n_of_decompositions (1,1) double { mustBeInteger, mustBePositive } = 10;

        ramus.ramus_hyperprior (1,1) string { mustBeMember( ...
            ramus.ramus_hyperprior, ...
            ["spatially balanced", "spatially constant"]) ...
        } = "spatially balanced";

        ramus.ramus_signal_to_noise_ratio (1,1) double { mustBeReal } = 30;

        ramus.ramus_ias_map_iterations (1,1) double { mustBeInteger, mustBePositive } = 10;

        ramus.ramus_sampling_frequency (1,1) double { mustBePositive } = 1025;

        ramus.ramus_low_cut_frequency (1,1) double { mustBePositive } = 7;

        ramus.ramus_high_cut_frequency (1,1) double { mustBePositive } = 9;

        ramus.ramus_time_start (1,1) double { mustBeNonnegative } = 0;

        ramus.ramus_time_window (1,1) double { mustBeNonnegative } = 0;

        ramus.ramus_n_of_time_steps (1,1) double { mustBeInteger, mustBePositive } = 1;

        ramus.ramus_time_step (1,1) double { mustBeInteger, mustBeNonnegative } = 0;

        ramus.ramus_data_normalization (1,1) string { mustBeMember( ...
            ramus.ramus_data_normalization, ...
            [ ...
                "maximum entry", ...
                "maximum column norm", ...
                "average column norm", ...
                "none" ...
            ] ...
        ) } = "maximum entry";

        ramus.ramus_initial_guess_mode (1,1) string { mustBeMember( ...
            ramus.ramus_initial_guess_mode, ...
            [ ...
                "hyperparameter value", ...
                "hyperprior expectation", ...
            ] ...
        ) } = "hyperparameter value";

        % TODO SESAME parameters.

        sesame.sesame = []

    end % arguments

    % Initialize return values as empty.

    rec_vec_position = [];
    rec_vec_angle = [];
    rec_vec_magnitude = [];

    % Perform a multigrid decomposition around given source dipole positions.

    [multigrid_dec, multigrid_ind, multigrid_perm] = make_multigrid_dec(zef.source_positions, n_subset, 1, 1);

    multigrid_dec = multigrid_dec{1}{1}{1};
    multigrid_perm = multigrid_perm{3};

    % Reconstruction vectors.

    rec_vec_position = zeros(3*length(multigrid_dec),1);
    rec_vec_angle = zeros(3*length(multigrid_dec),1);
    rec_vec_magnitude = zeros(3*length(multigrid_dec),1);

    % Use a standard Cartesian basis.

    source_dirs = eye(3);

    % Initialize waitbar

    wbtitle = "Dipole localization";

    wb = waitbar(0, wbtitle + ".");

    n_of_iters = length(multigrid_dec);

    % Create cleanup object for handling exceptions.

    cleanupObj = onCleanup(@() cleanup_dipole_localization(wb));

    % Set graphics root visibility to off, so waitbars created by the
    % inverse_methods are not shown.

    set(groot, 'DefaultFigureVisible', 'off');

    % Initialize inverse method fields in global zef instance.

    waitbar(0, wb, wbtitle + ": inverse method initialization.");

    zef = goodness_of_inversion.init_inverse_method(zef, inverse_method_name, mne, beamformer, ramus, sesame);

    waitbar(1, wb);

    % Go over the subsets of the domain and find out the differences between the
    % reconstruction and the (possibly synthetic) comparison dipole data.

    init_time_val = now;

    waitbar(0, wb, wbtitle + ": dipole localization.");

    for i = 1 : n_of_iters

        for j = 1 : 3

            source_dir = source_dirs(j,:);

            source_pos = zef.source_positions(multigrid_dec(i),:);

            zef.inv_synth_source(1,1:3) = source_pos;
            zef.inv_synth_source(1,4:6) = source_dir;

            % Perform a minimum norm estimate and save them.

            meas_data = find_source_legacy_fn(zef);

            zef.measurements = meas_data;

            [zef, rec] = goodness_of_inversion.call_inverse_method(zef, inverse_method_name);

            if isempty(rec)
                rec_vec_position = [];
                rec_vec_angle = [];
                rec_vec_magnitude = [];
                return
            end

            if iscell(rec)
                rec = rec{1};
            end

            reconstruction = zeros(size(zef.source_positions, 1) * 3, 1);

            reconstruction(1:length(rec(:))) = rec;

            % Calculate the differences between sources and reconstruction.

            [position_diff, angle_diff, magnitude_diff] = goodness_of_inversion.roi_dipole_accuracy( ...
                zef.inv_synth_source(1,1:3)                                        ...
            ,                                                                      ...
                zef.inv_synth_source(1,4:6)                                        ...
            ,                                                                      ...
                reconstruction                                                     ...
            ,                                                                      ...
                zef.source_positions                                               ...
            ,                                                                      ...
                roi_radius                                                         ...
            );

            rec_vec_position(3*(i-1)+j) = position_diff;
            rec_vec_angle(3*(i-1)+j) = angle_diff;
            rec_vec_magnitude(3*(i-1)+j) = magnitude_diff;

        end

        if mod(i, ceil(n_of_iters/5000))==0

            time_val = now - init_time_val;

            waitbar( ...
                i/n_of_iters, ...
                wb, ...
                wbtitle ...
                + ": " ...
                + num2str(i) ...
                + " / " ...
                + num2str(n_of_iters) ...
                + ". Ready approx: " ...
                + datestr(datevec(now+(n_of_iters/i - 1)*time_val)) ...
                + "." ...
            );

        end

    end

    % Reshape the data to make it fit into Zeffiro's fields.

    rec_vec_position = reshape(rec_vec_position/sqrt(3),3,length(multigrid_dec));
    rec_vec_angle = reshape(rec_vec_angle/sqrt(3),3,length(multigrid_dec));
    rec_vec_magnitude = reshape(rec_vec_magnitude/sqrt(3),3,length(multigrid_dec));

    % These get returned after this reordering

    rec_vec_position = rec_vec_position(:, multigrid_perm);
    rec_vec_angle = rec_vec_angle(:, multigrid_perm);
    rec_vec_magnitude = rec_vec_magnitude(:, multigrid_perm);

end

%% Helper functions

function cleanup_dipole_localization(wb)

% This is called when the cleabup object in the symbol table of
% main, when the function finishes running.

    close(wb);
    set(groot, 'DefaultFigureVisible', 'on');

end


function [meas_data] = find_source_legacy_fn(zef)

% find_source_legacy_fn
%
% Copied and turned into a proper function from find_source. Generates
% synthetic measurement data from the lead field contained in a iven zef
% instance.
%

    source_positions = zef.source_positions;

    noise_level = zef.inv_synth_source(1,8);

    s_p = zef.inv_synth_source(:,1:3);

    s_o = zef.inv_synth_source(:,4:6);

    s_o = s_o./repmat(sqrt(sum(s_o.^2,2)),1,3);

    s_a = zef.inv_synth_source(:,7);

    s_f = 1e-3*repmat(s_a,1,3).*s_o;

    L = zef.L;

    meas_data = zeros(size(L(:,1),1),1);

    for i = 1 : size(s_p,1)

        [s_min,s_ind] = min(sqrt(sum((source_positions - repmat(s_p(i,:),size(source_positions,1),1)).^2,2)));

        meas_data = meas_data + s_f(i,1)*L(:,3*(s_ind-1)+1) + s_f(i,2)*L(:,3*(s_ind-1)+2) + s_f(i,3)*L(:,3*(s_ind-1)+3);
    end

    n_val = max(abs(meas_data));

    meas_data = meas_data + noise_level*max(abs(meas_data)).*randn(size(meas_data));

end

function [multigrid_dec, multigrid_ind, multigrid_perm] = make_multigrid_dec(center_points,n_subset,n_decs,n_levels)

% Creates a multigrid decomposition, indices and permutations from given center
% points, number of subdivisions of the domain, number of deccompositions and
% the amouns of times the subsets are broken down into further decompositions.

    multigrid_dec = cell(0);
    multigrid_ind = cell(0);
    c_ind = [1 : size(center_points,1)]';
    multigrid_perm_1 = zeros(size(center_points,1),1);
    multigrid_perm_2 = zeros(n_decs*n_levels*n_subset,1);
    ind_counter = 0;

    for j = 1 : n_levels

        s_ind = cell(0);

        for i = 1 : n_decs

            randperm_vec = randperm(length(c_ind));

            if j == 1
                randperm_ind = randperm_vec(1:n_subset);
            else
                randperm_ind = randperm_vec(1:(n_subset-1)*n_subset^(j-1));
            end

            s_ind{i} = c_ind(randperm_ind);
            multigrid_ind{i}{j} = s_ind{i};
            c_ind = setdiff(c_ind, s_ind{i});
            multigrid_perm_1(s_ind{i}) = [ind_counter+1:ind_counter + length(s_ind{i})]';
            multigrid_perm_2(ind_counter+1:ind_counter + length(s_ind{i})) = s_ind{i};
            ind_counter = ind_counter + length(s_ind{i});
        end

        aux_cell = cell(0);

        for i = 1 : n_decs

            if j == 1

                aux_cell{1} = s_ind{i};
                multigrid_dec{i}{j} = aux_cell;

            else

                MdlKDT = KDTreeSearcher(center_points(s_ind_old{i},:));
                aux_ind = knnsearch(MdlKDT,center_points(s_ind{i},:));
                unique_aux_ind = [1 : length(s_ind_old{i})]';

                for k = 1 : length(unique_aux_ind)
                    aux_cell_ind = find(aux_ind == unique_aux_ind(k));
                    aux_cell{k} = s_ind{i}(aux_cell_ind);
                end

                multigrid_dec{i}{j} = aux_cell;

            end
        end

        s_ind_old = s_ind;

    end

    MdlKDT = KDTreeSearcher(center_points(multigrid_perm_2,:));
    multigrid_perm_3 = knnsearch(MdlKDT,center_points);

    multigrid_perm = cell(0);
    multigrid_perm{1} = multigrid_perm_1;
    multigrid_perm{2} = multigrid_perm_2;
    multigrid_perm{3} = multigrid_perm_3;

end

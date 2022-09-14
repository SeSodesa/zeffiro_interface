function [zef, rec_vec_position, rec_vec_angle, rec_vec_magnitude] = main( ...
    zef,           ...
    n_subset,      ...
    roi_radius,    ...
    inverse_method ...
)

% zef_dipole_localization_map
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
% - inverse_method
%
%   A function handle to the inverse method used to calculate the
%   reconstruction from the lead field. Currently supported inversion methods
%   are
%
%   - @zef_find_mne_reconstruction
%   - @zef_beamformer
%   - @zef_ramus_iteration
%
%   Methods that might be supported later include
%
%   - @SESAME_inversion (problems with randomization (?) causing too large vector indices)
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

        inverse_method (1,1) function_handle

    end

    % Initialize return values as empty.

    rec_vec_position = [];
    rec_vec_angle = [];
    rec_vec_magnitude = [];

    % Check for initial errors

    ACCEPTED_METHODS = {@zef_find_mne_reconstruction, @zef_beamformer, @zef_ramus_iteration};

    unknown_method = true;

    for mind = 1 : length(ACCEPTED_METHODS)

        if isequal(inverse_method, ACCEPTED_METHODS{mind})

            unknown_method = false;

            break

        end

    end

    if unknown_method

        warning("zef_dipole_localization_map was given an unknown inversion method. Aborting.");

        return

    end

    % Perform a multigrid decomposition around given source dipole positions.

    [multigrid_dec, multigrid_ind, multigrid_perm] = zef_make_multigrid_dec(zef.source_positions, n_subset, 1, 1);

    multigrid_dec = multigrid_dec{1}{1}{1};
    multigrid_perm = multigrid_perm{3};

    % Reconstruction vectors.

    rec_vec_position = zeros(3*length(multigrid_dec),1);
    rec_vec_angle = zeros(3*length(multigrid_dec),1);
    rec_vec_magnitude = zeros(3*length(multigrid_dec),1);

    % Use a standard Cartesian basis.

    source_dirs = eye(3);

    % Initialize waitbar

    wbtitle = 'Dipole localization';

    wb = waitbar(0, [wbtitle '.']);

    n_of_iters = length(multigrid_dec);

    % Create cleanup object for handling exceptions.

    cleanupObj = onCleanup(@() zef_cleanup_dipole_localization(wb));

    % Set graphics root visibility to off, so waitbars created by the
    % inverse_methods are not shown.

    set(groot, 'DefaultFigureVisible', 'off');

    % Initialize inverse method fields in global zef instance.

    waitbar(0, wb, strcat(wbtitle, ': inverse method initialization.'));

    zef = zef_init_inverse_method(zef, inverse_method);

    waitbar(1, wb);

    % Go over the subsets of the domain and find out the differences between the
    % reconstruction and the (possibly synthetic) comparison dipole data.

    init_time_val = now;

    waitbar(0, wb, strcat(wbtitle, ': dipole localization.'));

    for i = 1 : n_of_iters

        for j = 1 : 3

            source_dir = source_dirs(j,:);

            source_pos = zef.source_positions(multigrid_dec(i),:);

            zef.inv_synth_source(1,1:3) = source_pos;
            zef.inv_synth_source(1,4:6) = source_dir;

            % Perform a minimum norm estimate and save them.

            meas_data = find_source_legacy_fn(zef);

            zef.measurements = meas_data;

            [zef, rec] = zef_call_inverse_method(zef, inverse_method);

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

            [position_diff, angle_diff, magnitude_diff] = zef_roi_dipole_accuracy( ...
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
                strcat( ...
                    wbtitle, ...
                    ': ', ...
                    num2str(i), ...
                    ' / ', ...
                    num2str(n_of_iters), ...
                    '. Ready approx: ', ...
                    datestr(datevec(now+(n_of_iters/i - 1)*time_val)), ...
                    '.' ...
                ) ...
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

function zef_cleanup_dipole_localization(wb)

% This is called when the cleabup object in the symbol table of
% zef_dipole_localization_map, when the function finishes running.

    close(wb);
    set(groot, 'DefaultFigureVisible', 'on');

end


function [meas_data] = find_source_legacy_fn(zef)

% find_source_legacy_fn
%
% Copied and turned into a proper function from zef_find_source. Generates
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

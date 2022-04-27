function [rec_vec_position, rec_vec_angle, rec_vec_magnitude] = zef_dipole_localization_map( ...
    zef,           ...
    n_subset,      ...
    roi_radius,    ...
    inverse_method ...
)

% Calculates a reconstruction of a dipole source from a lead field L inside a
% Region of Interest (RoI) around the source dipole. Returns the reconstruction
% dipole positions, angles and magnitudes.
%
% Input: zef for reading values like model soruce positions, number of
% subdivisions of the RoI, the radius of the RoI and the inverse method used to
% calculate the reconstruction from the lead field. Currently supported
% inversion methods are
%
% - @zef_find_mne_reconstruction
% - @zef_beamformer
% - @zef_ramus_iteration
%
% Methods that might be supported later include
%
% - @SESAME_inversion (problems with randomization (?) causing too large vector indices)
%
% Output: three reconstruction arrays: dipole positions, angles and magnitudes.
% These will be empty in case of errors, so please check for that when calling
% the routine.

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

    if nargout ~= 3
        warning('zef_dipole_localization_map must be called with 3 output arguments');
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

    zef_init_inverse_method(inverse_method);

    waitbar(1, wb);

    % Go over the subsets of the domain and find out the differences between the
    % reconstruction and the (possibly synthetic) comparison dipole data.

    init_time_val = now;

    waitbar(0, wb, strcat(wbtitle, ': dipole localization.'));

    for i = 1 : n_of_iters

        for j = 1 : 3

            source_dir = source_dirs(j,:);
            source_pos = zef.source_positions(multigrid_dec(i),:);

            % This kind of global state manipulation is from the Satan. Do not
            % do this, ever.

            assignin('base', 'sd', source_dir);
            assignin('base', 'sp', source_pos);
            evalin('base', 'zef.inv_synth_source(1,4:6) = sd;');
            evalin('base', 'zef.inv_synth_source(1,1:3) = sp;');

            zef.inv_synth_source(1,4:6) = source_dir;
            zef.inv_synth_source(1,1:3) = source_pos;

            % Perform a minimum norm estimate.

            meas_data = zef_find_source_legacy;

            % Again, global state manipulation. Only the devil could have
            % though up something so evil.

            assignin('base', 'meas_data', meas_data);
            evalin('base', 'zef.measurements = meas_data;');

            [rec] = zef_call_inverse_method(zef, inverse_method);

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

            % Whoo boy. We really are in hell, aren't we. Global state
            % manipulation, again.

            assignin('base', 'reconstruction', reconstruction);
            evalin('base', 'zef.reconstruction = reconstruction;');

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

function zef_cleanup_dipole_localization(wb)

% This is called when the cleabup object in the symbol table of
% zef_dipole_localization_map, when the function finishes running.

    close(wb);
    set(groot, 'DefaultFigureVisible', 'on');

end

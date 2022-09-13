function [position_diff, angle_diff, magnitude_diff] = zef_roi_dipole_accuracy( ...
    dipole_pos       ...
,                    ...
    dipole_vec       ...
,                    ...
    reconstruction   ...
,                    ...
    source_positions ...
,                    ...
    roi_radius       ...
)

% Calculates the difference between positions, angles and magnitudes of a given
% dipole potential field reconstruction (formed via some inverse method) and an
% actual dipole potential field, in a region of interest (RoI).

    arguments

        dipole_pos (:,3) double

        dipole_vec (:,3) double

        reconstruction (:,1) double

        source_positions (:,3) double

        roi_radius (1,1) double { mustBeReal, mustBePositive }

    end

    % Determine the ball that makes up the RoI around source_positions, based
    % on given dipole positions and RoI radius.

    roi_positions = rangesearch(source_positions, dipole_pos, roi_radius);
    roi_positions = roi_positions{1};

    % For whatever reason, the given reconstruction might be in a cell array
    % (a.k.a. tuple), so decompose appropriately.

    if iscell(reconstruction)
        reconstruction = reshape(reconstruction{1}, 3, length(reconstruction(:))/3);
    else
        reconstruction = reshape(reconstruction, 3, length(reconstruction(:))/3);
    end

    % Find out which parts of reconstruction are within the region of interest.

    reconstruction = reconstruction(:,roi_positions);

    % Find out reconstruction positions, directions and magnitudes by averaging
    % within the region of interest.

    rec_pos = sum(source_positions(roi_positions,:)'.*abs(reconstruction),2)./sum(abs(reconstruction),2);
    rec_dir = mean(reconstruction,2);
    rec_mag = norm(rec_dir,2);
    rec_dir = rec_dir / rec_mag;

    % From the missing comparison quantities for the "actual" dipoles.

    dipole_mag = norm(dipole_vec(:),2);
    dipole_dir = dipole_vec/dipole_mag;

    % Form the differences being asked for

    position_diff = norm(rec_pos - dipole_pos(:));
    angle_diff = acos(sum(rec_dir.*dipole_dir(:)))*180/pi;
    magnitude_diff = 1 - rec_mag / dipole_mag;

end

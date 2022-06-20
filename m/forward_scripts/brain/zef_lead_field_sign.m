function Lsign = zef_lead_field_sign( ...
    source_positions, ...
    electrodes, ...
    L ...
)

    % Calculates the sign of the given lead field L based on source positions,
    % their directions, the electrode positions and the lead field itself.

    arguments
        source_positions (:, 3) double
        electrodes
        L
    end

    % Source closest to origin (SCO).

    source_distances = zef_L2_norm(source_positions, 2);

    [~, scoind] = min(source_distances);

    scopos = source_positions(scoind,:);

    % Electrode with greatest z-coordinate (EGZ)

    [~, egzind] = max(electrodes(:,3));

    % Indices of L that match the SCO

    sco_xyz_inds = 3 * (scoind-1) + 1 : 3 * scoind;

    % Part of L that matches both SCO and EGZ

    scoegzL = L(egzind, sco_xyz_inds);

    % Scalar potential of SCO at EGZ

    scoegzu = scoegzL * scopos';

    % Calculate the sign

    try
        Lsign = sign(scoegzL(3) / scoegzu);
    catch
        warning('Sign of L could not be determined because of zero division. Setting it as (+)')
        Lsign = 1;
    end

end

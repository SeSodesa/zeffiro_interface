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

    % Electrode whose projection onto the xy-plane is closest to that of SCO,
    % as in is "most below/above it". Call it projectively closest electrode
    % (PCE).

    sco_ele_diffs = scopos - electrodes;

    sco_ele_diffs_xy = sco_ele_diffs(:, 1:2);

    pce_dists = zef_L2_norm(sco_ele_diffs_xy, 2);

    [~, pceind] = min(pce_dists);

    % Indices of L that match the SCO.

    sco_xyz_inds = 3 * (scoind-1) + 1 : 3 * scoind;

    % Part of L that matches both SCO and PCE.

    scopceL = L(pceind, sco_xyz_inds);

    % Scalar potential of SCO at EGZ.

    scopceu = scopceL * scopos';

    % Calculate the sign

    try
        Lsign = sign(scopceL(3) / scopceu);
    catch
        warning('Sign of L could not be determined because of zero division. Setting it as (+)');
        Lsign = 1;
    end

end

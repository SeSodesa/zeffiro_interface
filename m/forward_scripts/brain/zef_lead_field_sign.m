function Lsign = zef_lead_field_sign( ...
    source_positions, ...
    source_directions, ...
    electrodes, ...
    L ...
)

    % Calculates the sign of the given lead field L based on source positions,
    % their directions, the electrode positions and the lead field itself.

    arguments
        source_positions (:, 3) double
        source_directions (:, 3) double
        electrodes
        L
    end

    % Source closest to origin (SCO).

    source_distances = zef_L2_norm(source_positions, 2);

    [~, scoind] = min(source_distances);

    scopos = source_positions(scoind,:);
    scodir = source_directions(scoind, :);

    % Electrode with greatest z-coordinate (EGZ)

    [~, egzind] = max(electrodes(:,3));

    % Indices of L that match the SCO

    scoLinds = 3 * (scoind-1)+1 : 3 * scoind;

    % Part of L that matches both SCO and EGZ

    scoegzL = L(egzind, scoLinds);

    % Scalar potential of SCO at EGZ

    scoegzu = scoegzL * scopos';

    % Calculate the sign

    try
        Lsign = scodir(3) / scoegzu;
    catch
        warning('Sign of L could not be determined because of zero division. Setting it as (+)')
        Lsign = 1;
    end

end

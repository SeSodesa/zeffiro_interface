function [G, interpolation_positions] = zef_hdiv_interpolation( ...
    p_nodes, ...
    p_tetrahedra, ...
    p_brain_inds, ...
    p_intended_source_inds ...
)

    % Produces a lead field interpolation matrix G with position-based
    % optimization (PBO), based on the H(div) (face-intersecting + edgewise)
    % source model. Also returns the related interpolation positions.

    arguments
        p_nodes (:,3) double {mustBeNonNan}
        p_tetrahedra (:,4) double {mustBeInteger, mustBeNonnegative}
        p_brain_inds (:,1) double {mustBeInteger, mustBeNonnegative}
        p_intended_source_inds (:,1) double {mustBeInteger, mustBeNonnegative}
    end

    G = [];

    % Open up a waitbar.

    wbtitle = 'Lead field H(div) interpolation';
    wb = waitbar(0, wbtitle);

    % Define cleanup operations, in case of an interruption.

    cleanupfn = @(wb) close(wb);
    cleanupobj = onCleanup(@() cleanupfn(wb));

    % Dipoles and their adjacency and weight matrices T and G.

    [T_fi, G_fi, ~, fi_source_directions, fi_source_positions, ~] = zef_fi_dipoles( ...
        p_nodes, ...
        p_tetrahedra, ...
        p_brain_inds ...
    );

    [T_ew, G_ew, ~, ew_source_directions, ew_source_positions, ~] = zef_ew_dipoles( ...
        p_nodes, ...
        p_tetrahedra, ...
        p_brain_inds ...
    );

    % Form local environment indices based on adjacency matrix T_fi. TODO:
    % change this to use some other condition than .

    valid_source_inds = full(find(sum(T_fi) >= 4))';
    valid_source_inds = intersect(valid_source_inds, p_intended_source_inds);

    % Restrict stensils to intended source positions.

    T_fi = T_fi(:, valid_source_inds);
    T_ew = T_ew(:, valid_source_inds);

    % Form interpolation positions (barycenters of tetrahedra).

    source_tetra = p_tetrahedra(valid_source_inds,:);
    interpolation_positions = zef_tetra_barycentra(p_nodes, source_tetra);

    % Form initial values based on given nodes, tetrahedra and lead field.

    n_of_iterations = size(valid_source_inds, 1);

    print_interval = ceil(n_of_iterations / 100);

    % Initialize weight matrix.

    G_rows = size(p_nodes, 1);
    G_cols = 3 * size(interpolation_positions, 1);

    G = sparse(G_rows, G_cols, 0);

    % Start iteration over the source positions of interest.

    tic;

    % Generate coefficient matrix indices.

    S_fi_cell = cell(0);
    S_ew_cell = cell(0);

    for i = 1 : n_of_iterations

        % Find local neighbour indices.

        fi_neighbour_inds = full(find(T_fi(:,i)));
        ew_neighbour_inds = full(find(T_ew(:,i)));

        % N of non-zero object function coefficients.

        n_coeff_fi = numel(fi_neighbour_inds);
        n_coeff_ew = numel(ew_neighbour_inds);
        n_coeff = n_coeff_fi + n_coeff_ew;

        % Dipole locations and directions.

        dir_mat = [ ...
            fi_source_directions(fi_neighbour_inds,:) ; ...
            ew_source_directions(ew_neighbour_inds,:) ...
        ];

        loc_mat = [ ...
            fi_source_positions(fi_neighbour_inds,:) ; ...
            ew_source_positions(ew_neighbour_inds,:) ...
        ];

        % PBO weigth coefficients from differences between barycentra
        % (interpolation positions) and dipole positions.

        interp_pos = interpolation_positions(i, :);
        interp_pos = repmat(interp_pos, n_coeff, 1);

        pos_diffs = loc_mat - interp_pos;

        weight_coefs = zef_L2_norm(pos_diffs, 2);

        % Position-based optimization matrix.

        PBO_mat = [                     ...
            diag(weight_coefs) dir_mat; ...
            dir_mat' zeros(3,3)         ...
        ];

        % Solve for Lagrangian multipliers.

        Coeff_mat = PBO_mat \ [zeros(n_coeff,3); eye(3)];

        % Accumulate interpolation matrix.

        col_inds = 3 * (i-1) + 1 : 3 * i;

        % G(:, col_inds) = ...
        %     G_fi(:,fi_neighbour_inds) ...
        %     * ...
        %     Coeff_mat(1:n_coeff_fi,:) ...
        %     + ...
        %     G_ew(:,ew_neighbour_inds) ...
        %     * ...
        %     Coeff_mat(n_coeff_fi+1:n_coeff,:) ...
        % ;

        % Row indices

        S_fi_cell{i,1} = repmat(fi_neighbour_inds, 1, size(Coeff_mat,2));
        S_ew_cell{i,1} = repmat(ew_neighbour_inds, 1, size(Coeff_mat,2));

        % Column indices

        S_fi_cell{i,2} = repmat(col_inds', length(fi_neighbour_inds), 1);
        S_ew_cell{i,2} = repmat(col_inds', length(ew_neighbour_inds), 1);

        % Values

        S_fi_cell{i,3} = Coeff_mat(1 : n_coeff_fi, :);
        S_ew_cell{i,3} = Coeff_mat(n_coeff_fi+1 : n_coeff, :);

        % Update waitbar.

        if mod(i, print_interval) == 0

            time_val = toc;

            waitbar( ...
                i/n_of_iterations, ...
                wb, ...
                [ ...
                    wbtitle, ...
                    ' (', ...
                    num2str(i), ...
                    ' / ', ...
                    num2str(n_of_iterations), ...
                    '). Ready: ' datestr(datevec(now+(n_of_iterations/i - 1)*time_val/86400)) '.' ...
                ] ...
            );

        end
    end

    % Calculate how many indices we got per column.

    entry_counter_fi = 0;
    entry_counter_ew = 0;

    for i = 1 : n_of_iterations

        n_of_fi_vals = numel(S_fi_cell{i, 3});
        n_of_ew_vals = numel(S_ew_cell{i, 3});

        entry_counter_fi = entry_counter_fi + n_of_fi_vals;
        entry_counter_ew = entry_counter_ew + n_of_ew_vals;

    end

    n_of_entries_fi = entry_counter_fi;
    n_of_entries_ew = entry_counter_ew;

    entry_counter_fi = 0;
    entry_counter_ew = 0;

    I_fi = zeros(n_of_entries_fi,1);
    J_fi = zeros(n_of_entries_fi,1);
    K_fi = zeros(n_of_entries_fi,1);

    I_ew = zeros(n_of_entries_ew,1);
    J_ew = zeros(n_of_entries_ew,1);
    K_ew = zeros(n_of_entries_ew,1);

    for i = 1 : n_of_iterations

        I_fi(entry_counter_fi + 1 : entry_counter_fi + numel(S_fi_cell{i,1})) = S_fi_cell{i,1}(:);
        J_fi(entry_counter_fi + 1 : entry_counter_fi + numel(S_fi_cell{i,2})) = S_fi_cell{i,2}(:);
        K_fi(entry_counter_fi + 1 : entry_counter_fi + numel(S_fi_cell{i,3})) = S_fi_cell{i,3}(:);

        entry_counter_fi = entry_counter_fi + numel(S_fi_cell{i,1});

        I_ew(entry_counter_ew + 1 : entry_counter_ew + numel(S_ew_cell{i,1})) = S_ew_cell{i,1}(:);
        J_ew(entry_counter_ew + 1 : entry_counter_ew + numel(S_ew_cell{i,2})) = S_ew_cell{i,2}(:);
        K_ew(entry_counter_ew + 1 : entry_counter_ew + numel(S_ew_cell{i,3})) = S_ew_cell{i,3}(:);


        entry_counter_ew = entry_counter_ew + numel(S_ew_cell{i,1});

    end

    S_fi = sparse(I_fi, J_fi, K_fi, size(G_fi,2), G_cols);
    S_ew = sparse(I_ew, J_ew, K_ew, size(G_ew,2), G_cols);

    G = G_fi * S_fi + G_ew * S_ew;

    waitbar(1, wb);

end

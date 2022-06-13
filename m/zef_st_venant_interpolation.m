function [G, interpolation_positions] = zef_st_venant_interpolation( ...
    p_nodes, ...
    p_tetrahedra, ...
    p_brain_ind, ...
    p_intended_source_inds, ...
    p_regparam ...
)

    % Produces an interpolation matrix G in a tetrahedral finite element mesh
    % with the St. Venant method.

    G = [];
    interpolation_positions = [];

    % Open up a waitbar

    wbtitle = 'Lead field interpolation (St. Venant)';
    wb = waitbar(0, wbtitle);

    % Define cleanup operations, in case of an interruption.

    cleanupfn = @(handle) close(handle);
    cleanupobj = onCleanup(@() cleanupfn(wb));

    % Define adjacency matrix for tetrahedra

    adjacency_mat = zef_adjacency_matrix(p_nodes, p_tetrahedra(p_brain_ind,:));

    %% First find nodes closest to the given positions.

    % Nearest nodes for each interpolation position with KDTree search

    MdlKDT = KDTreeSearcher(p_nodes);

    interpolation_positions = zef_tetra_barycentra(p_nodes, p_tetrahedra(p_intended_source_inds, :));

    center_node_inds = knnsearch(MdlKDT, interpolation_positions);

    % Storage for the interpolation results.

    n_of_iters = numel(center_node_inds);

    %% Interpolation for each position

    waitbar(0, wb, [wbtitle, ': interpolation']);

    wbtitleloop = [wbtitle, ': interpolation'];

    % Initialize interpolation weight matrix G

    G = sparse(size(p_nodes,1), 3 * size(interpolation_positions, 1), 0);

    % Cartesian directions for interpolation

    basis = eye(3);

    for ind = 1 : n_of_iters

        waitbar(ind / n_of_iters, wb, [wbtitleloop, num2str(ind), ' / ', num2str(n_of_iters)]);

        % Fetch reference node coordinates

        refnode_ind = center_node_inds(ind);
        refnode = p_nodes(refnode_ind, :);

        % Calculate the distances from refnode with Matlab's broadcasting
        % mechanism and save them to the preallocated distance matrix.

        neighbour_inds = find(adjacency_mat(:, refnode_ind));

        % Some nodes apparently have no neighbours. Go figure...

        if isempty(neighbour_inds)
            continue
        end

        neighbour_inds = setdiff(neighbour_inds, refnode_ind);
        n_of_neighbours = numel(neighbour_inds); % + 1;
        neighbours = p_nodes(neighbour_inds, :);
        neighbour_diffs = neighbours - refnode;

        % Calculate the distances and longest edge.

        dists = zef_L2_norm(neighbour_diffs, 2);
        longest_edge_len = max(dists, [], 'all');

        %% Construct restriction matrices P, b and regularization matrix D.

        P = zeros(9, n_of_neighbours);

        % Conservation of charge (rows of ones)

        P(1:3:end,:) = repmat(ones(1, size(P, 2)), 3, 1);

        % Dipole moment approximations

        try
            P(2:3:end,:) = (1/longest_edge_len) * neighbour_diffs';
        catch
            keyboard;
        end
        % Suppression of higher order moments

        suppressor_b = dists'.^2 / longest_edge_len^2;
        suppressor = repmat(suppressor_b, 3, 1);

        P(3:3:end,:) = suppressor;

        b = zeros(9,3);
        b(2:3:end, :) = basis / longest_edge_len;

        D = diag(sum(dists.^2, 2));

        % Form the monopolar loads m via least squares approximation
        % m = inv(P'P + aD)P'b.

        m = inv(P' * P + p_regparam * D) * P' * b;

        for iind = 1 : 3

            % G([refnode_ind; neighbour_inds], 3 * (ind -1 ) + iind) = m(:, iind);
            G(neighbour_inds, 3 * (ind -1 ) + iind) = m(:, iind);

        end

    end

end

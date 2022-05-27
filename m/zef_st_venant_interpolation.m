function G = zef_st_venant_interpolation( ...
    p_nodes, ...
    p_tetrahedra, ...
    p_interpolation_positions, ...
    p_regparam, ...
    p_n_of_electrodes, ...
    p_source_nonzero_inds ...
    p_brain_ind ...
)

    % Interpolates a weight matrix G in a tetrahedral finite element mesh with
    % the St. Venant method. Returns a matrix in case of error.

    G = [];

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

    MdlKDT = KDTreeSearcher(nodes);
    center_node_inds = knnsearch(MdlKDT, interpolation_positions);

    % Storage for the interpolation results.

    n_of_iters = numel(cener_node_inds);

    weights = cell(n_of_iters);

    %% Interpolation for each position

    waitbar(0, wb, [wbtitle, ': interpolation'])

    wbtitleloop = [wbtitle, ': interpolation']

    % Initialize interpolation weight matrix G

    G = sparse(size(nodes,1), 3 * size(p_interpolation_positions, 1), 0);

    % Cartesian directions for interpolation

    basis = eye(3);

    for ind = 1 : n_of_iters

        waitbar(ind / n_of_iters, wb, [wbtitleloop, num2str(ind), ' / ', num2str(n_of_iters)]);

        % Fetch reference node coordinates

        refnode_ind = center_node_inds(ind);
        refnode = nodes(refnode_ind, :);

        % Calculate the distances from refnode with Matlab's broadcasting
        % mechanism and save them to the preallocated distance matrix.

        neighbour_inds = find(adjacency_mat(:, refnode_ind));
        neighbour_inds = setdiff(neighbour_inds, refnode_ind);
        n_of_neighbours = numel(neighbour_inds) + 1;
        neighbours = nodes(neighbour_inds, :);
        neighbour_diffs = neighbours - refnode;

        % Calculate the distances and longest edge.

        dists = sqrt(sum(neighbour_diffs.^2, 2));
        longest_edge_len = max(dists, [], 'all');

        % Restriction matrices P, b and regularization matrix D.

        P = zeros(7, n_of_neighbours);
        P(1,:) = ones(1, size(P, 2)); % sum of m = 0
        P(2:4,2:end) = (1/longest_edge_len) * neighbour_diffs; % Moment conditions (3.44)
        P(5:7,2:end) = dists.^2 / longest_edge_len; % Moment conditions (3.45)

        b = zeros(9,3);
        b(2:3:end, :) = basis / longest_edge_len;

        D = diag(sum(diffs.^2, 2));

        % Form the monopolar loads m via least squares approximation
        % m = inv(P'P + aD)P'b.

        m = inv(P' * P + p_regparam * D) * P' * b;

        for iind = 1 : 3

            G([refnode_inds; neighbour_inds], 3 * (ind -1 ) + iind) = m(:, iind);

        end


    end

end

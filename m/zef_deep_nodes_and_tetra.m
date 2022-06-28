function [ ...
    out_deep_nodes, ...
    out_deep_node_inds, ...
    out_deep_tetra, ...
    out_deep_tetra_inds ...
] = zef_deep_nodes_and_tetra( ...
    in_nodes, ...
    in_tetra, ...
    in_volume_inds, ...
    in_acceptable_depth_mm ...
)

    % Documentation
    %
    % Produces the nodes and tetra which are deep enough within a given
    % volume, and also their indices in the given global tetra and node data
    % structures.
    %
    % Input:
    %
    % - in_nodes: finite elements nodes.
    %
    % - in_tetra: finite element tetrahedra (quadruples of node indices)
    %   constructed from above nodes.
    %
    % - in_volume_inds: the indices of the tetrahedra that form the volume
    %   under observation.
    %
    % - in_acceptable_depth_mm: the depth in millimetres, within which the
    %   deep tetra are located.
    %
    % Output:
    %
    % - out_tetra_ind: the indices of the tetrahedra that are deep enough
    %   inside the given brain segment.

    arguments
        in_nodes (:,3) double
        in_tetra (:,4) double {mustBeInteger, mustBePositive}
        in_volume_inds (:,1) double {mustBeInteger, mustBePositive}
        in_acceptable_depth_mm (1,1) double {mustBeNonnegative}
    end

    % Set empty return values.

    out_deep_node_inds = [];
    out_deep_nodes = [];
    out_deep_tetra_inds = [];
    out_deep_tetra = [];

    % Find out the boundary elements (triangles) of the volume.

    volume_tetra = in_tetra(in_volume_inds, :);

    volume_node_inds = volume_tetra(:);

    unique_volume_node_inds = unique(volume_node_inds);

    volume_nodes = in_nodes(unique_volume_node_inds, :);

    [~, surface_nodes] = zef_surface_mesh(volume_tetra, in_nodes);

    % Find out non-surface nodes deep enough within the volume with
    % rangesearch.

    non_surface_nodes = setdiff(volume_nodes, surface_nodes, 'rows');

    node_inds_too_near_to_surface = zef_nearest_points( ...
        surface_nodes, ...
        non_surface_nodes, ...
        in_acceptable_depth_mm, ...
        'range' ...
    );

    out_deep_node_inds = setdiff( ...
        volume_node_inds, ...
        node_inds_too_near_to_surface, ...
        'rows' ...
    );

    out_deep_nodes = in_nodes(out_deep_node_inds, :);

    % Find tetra in the volume which only contain acceptable (deep) nodes.

    tetra_deep_node_info = ismember(in_tetra, out_deep_node_inds);

    tetra_info_row_sums = sum(tetra_deep_node_info, 2);

    out_deep_tetra_inds = find(tetra_info_row_sums == 4);

    out_deep_tetra = in_tetra(out_deep_tetra_inds, :);

end % zef_deep_tetra

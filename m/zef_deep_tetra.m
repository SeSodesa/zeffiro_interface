function out_tetra_ind = zef_deep_tetra( ...
    in_nodes, ...
    in_tetra, ...
    in_volume_inds, ...
    in_acceptable_depth_mm ...
)

    % Produces subset of the FEM tetrahedra, whose barycenters are deep enough
    % within a given volume.
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

    % Set empty return value.

    out_tetra_ind = [];

    % Find out the boundary elements (triangles) of the volume.

    volume_tetra = in_tetra(in_volume_inds, :);

    volume_node_indices = unique(volume_tetra(:));

    volume_nodes = in_nodes(volume_node_indices, :);

    [surface_node_inds, surface_nodes] = zef_surface_mesh(volume_tetra, in_nodes);

    % Find out non-surface nodes deep enough within the volume with
    % rangesearch.

    non_surface_nodes = setdiff(volume_nodes, surface_nodes, 'rows');

    node_inds_too_near_to_surface = zef_nearest_points( ...
        surface_nodes, ...
        non_surface_nodes, ...
        in_acceptable_depth_mm, ...
        'range' ...
    );

    nodes_too_near_to_surface = non_surface_nodes(node_inds_too_near_to_surface, :);

    acceptable_non_surface_nodes = setdiff( ...
        non_surface_nodes, ...
        nodes_too_near_to_surface, ...
        'rows' ...
    );

    % TODO: Find tetra in the volume which only contain acceptable nodes.

end % zef_deep_tetra

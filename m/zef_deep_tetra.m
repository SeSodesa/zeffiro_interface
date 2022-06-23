function out_tetra = zef_deep_tetra( ...
    in_nodes, ...
    in_volume_tetra, ...
    in_volume_inds, ...
    in_acceptable_depth_mm ...
)

    % Produces subset of the FEM tetrahedra, whose barycenters are deep enough
    % within a given volume.
    %
    % Input:
    %
    % - in_nodes: finite elements nodes.
    % - in_tetra: finite element tetrahedra (quadruples of node indices)
    %   constructed from above nodes.
    % - in_volume_inds: the indices of the tetrahedra that form the volume
    %   under observation.
    % - in_acceptable_depth_mm: the depth in millimetres, within which the deep
    %   tetra are located.
    %
    % Output:
    %
    % - o_tetra: the tetrahedra that are deep enough inside the given brain
    %   segment.

    arguments
        in_nodes (:,3) double
        in_volume_tetra (:,4) double {mustBeInteger, mustBePositive}
        in_volume_inds (:,1) double {mustBeInteger, mustBePositive}
        in_acceptable_depth_mm (1,1) double {mustBeNonnegative}
    end

    % Set empty return value.

    out_tetra = [];

    % First find out the surface triangles (triples of node indices) of the
    % volume determined by the volume indices.

    volume_tetra = in_volume_tetra(in_volume_inds, :);

    [volume_surface_triangles, ~, surf_tetra_inds] = zef_surface_mesh(volume_tetra);

    volume_surface_tetra = volume_tetra(surf_tetra_inds ,:);

    % An abbreviation variable for clearer indexing.

    vst = volume_surface_triangles;

    % Calculate the normed surface normals of the triangles.

    triangle_edges_1 =  in_nodes(vst(:, 2), :) - in_nodes(vst(:, 1), :);
    triangle_edges_2 =  in_nodes(vst(:, 3), :) - in_nodes(vst(:, 1), :);

    surface_normals = cross(triangle_edges_1, triangle_edges_2);
    surface_normals = surface_normals ./ zef_L2_norm(surface_normals, 2);

end % zef_deep_tetra

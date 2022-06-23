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
    % - p_nodes: finite elements nodes.
    % - p_tetra: finite element tetrahedra (quadruples of node indices)
    %   constructed from above nodes.
    % - p_volume_inds: the indices of the tetrahedra that form the volume
    %   under observation.
    % - p_acceptable_depth_mm: the depth in millimetres, within which the deep
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

    % An abbreviation variable for clearer indexing.

    vst = volume_surface_triangles;

    % Get surface tetrahedra

    surface_tetra = surface_tetra_fn(vst, volume_tetra);

    % Calculate the normed surface normals of the triangles.

    triangle_edges_1 =  in_nodes(vst(:, 2), :) - in_nodes(vst(:, 1), :);
    triangle_edges_2 =  in_nodes(vst(:, 3), :) - in_nodes(vst(:, 1), :);

    surface_normals = cross(triangle_edges_1, triangle_edges_2);
    surface_normals = surface_normals ./ zef_L2_norm(surface_normals, 2);

end % zef_deep_tetra

function out_surface_tetra = surface_tetra_fn(in_surface_triangles, in_volume_tetra)

    % Finds out which tetrahedra the given surface triangles participate in.
    %
    % Input:
    %
    % - in_surface_triangles: the surface triangles.
    % - in_volume_tetra: the tetra from which the surface triangles are being sought.
    %   These must not contain any tetra not within the volume defined by the
    %   surface triangles.
    %
    % Output:
    %
    % - out_surface_tetra: the tetrahedra that participate in the given
    %   surface triangles.

    arguments
        in_surface_triangles (:,3) double { mustBeInteger, mustBePositive }
        in_volume_tetra (:,4) double { mustBeInteger, mustBePositive }
    end

    out_surface_tetra = [];

    % Sizes

    n_of_triangles = size(in_surface_triangles, 1);
    n_of_tetra = size(in_volume_tetra, 1);

    % Construct surface tetrahedra and sort containers for binary search.

    out_surface_tetra = zeros(n_of_triangles, 4);

    % Iterate over triangles and find the tetrahedra that share all 3 nodes
    % except one.

    for i = 1 : n_of_triangles

        disp(['Iteration ' num2str(i)]);

        triangle_i = in_surface_triangles(i,:);

        % A broadcasting (row-wise) check for set inclusion.

        tri_nodes_in_tetra = ismember(in_volume_tetra, triangle_i);

        % Row sums to determine membership.

        row_sums = sum(tri_nodes_in_tetra, 2);

        % Tetra indices that have contained the 3 nodes from the current
        % triangle.

        surface_tetra_inds = find(row_sums == 3);
        sti = surface_tetra_inds;

        % Size check to make sure the tetra fit into the dense output array.

        max_ind = max(sti, [], 'all');

        n_of_st = size(out_surface_tetra, 1);

        if max_ind > n_of_st

            warning('Allocating more space for surface tetra.')

            add_needed_space = max_ind - n_of_st;

            out_surface_tetra = [ out_surface_tetra ; zeros(add_needed_space, 4) ];

        end

        % After size check, insert tetra into output array

        out_surface_tetra(sti,:) = in_volume_tetra(sti,:);

    end % for

end % surface_tetra_fn

function o_tetra = zef_deep_tetra( ...
    p_nodes, ...
    p_tetra, ...
    p_volume_inds ...
)

    % Produces subset of the FEM tetrahedra, whose barycenters are deep enough
    % within a given volume.
    %
    % Input:
    %
    % - p_nodes: finite elements nodes.
    % - p_tetra: finite element tetrahedra constructed from above nodes.
    % - p_volume_inds: the indices of the tetrahedra that form the volume
    %   under observation
    %
    % Output:
    %
    % - o_tetra: the tetrahedra that are deep enough inside the given brain
    %   segment.

    arguments
        p_nodes (:,3) double
        p_tetra (:,4) double {mustBeInteger, mustBePositive}
        p_volume_inds (:,1) double {mustBeInteger, mustBePositive}
    end

end

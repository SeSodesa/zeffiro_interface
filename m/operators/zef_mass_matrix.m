function A = zef_mass_matrix(   ...
    nodes                       ...
,                               ...
    tetrahedra                  ...
,                               ...
    electrode_areas             ...
,                               ...
    electrode_indices           ...
,                               ...
    electrode_impedances        ...
,                               ...
    electrodes                  ...
)

% zef_mass_matrix: a two-dimensional version of zef_stiffness_matrix which,
% among other things, models the effects of any electrodes on the surface of a
% head model Ω. Here the mass matrix A is given by
%
%   ∑[ℓ=1,𝐿] 1 / (𝑍ℓ |𝑒ℓ|) ⋅ ∫[𝑒ℓ] ψᵢψⱼ d𝑆,
%
% where 𝑒ℓ is the set of 𝐿 electrodes, each with a surface area d𝑆, 𝑍ℓ is the
% impedance of an electrode, and ψᵢ and ψⱼ are scalar-valued basis functions
% of the discretized domain Ω.

    % Preallocate mass matrix A

    N = size(nodes,1);
    A = spalloc(N,N,0);

    % Calculate |𝑒ℓ| once

    n_of_el = numel(electrodes);

    n_of_tetra_nodes = 4;

    % Effective Contact Impedance

    ECI = electrode_impedances .* n_of_el;

    % Sum up individual integrands to form the mass matrix A

    for ii = 1 : n_of_tetra_nodes

        for jj = ii : n_of_tetra_nodes

            % Product of nodal basis functions ψᵢψⱼ = 1 ⟺ i = j, as the value
            % of a basis function ψₖ = 1 at only one node in a tetrahedron and
            % 0 at the other nodes.

            if isequal(ii,jj)
                integrand = 1;
            else
                integrand = 0;
            end

            % If integrand was one, select area element

            dS = integrand .* electrode_areas;

        end
    end
end

function A = zef_mass_matrix()

% zef_mass_matrix: a two-dimensional version of zef_stiffness_matrix which,
% among other things, models the effects of any electrodes on the surface of a
% head model Ω. Here the mass matrix A is given by
%
%   ∑[𝑒ℓ=1,𝐿] 1 / (𝑍ℓ |𝑒ℓ|) ⋅ ∫[𝑒ℓ] ψᵢψⱼ d𝑆,
%
% where each 𝑒ℓ is an electrode with a surface area d𝑆, 𝑍ℓ is the impedance of
% an electrode, and ψᵢ and ψⱼ are scalar-valued basis functions of the
% discretized domain Ω.

end

function A = zef_stiffness_matrix(basis, tensor)

% The stiffness matrix 𝐴 of a discretized function 𝑢ₕ = ∑ᵢ𝑧ᵢψᵢ, with each 𝑧ᵢ
% being a coordinate and ψᵢ a linear basis function, is defined by
%
%   𝐴[i,j] = ∫[Ω] ∇ψⱼ ⋅ (𝑇∇ψᵢ) d𝑉 ,
%
% where Ω is the entire domain and d𝑉 is a finite volume element from which
% the domain consists of, adn 𝑇 is a tensor, for which the equation
%
%   (𝑇∇u)⋅𝑛⃗ = 0
%
% holds for the non-discretized function 𝑢, on the boundary ∂Ω of the domain.

    error("Not yet implemented…")

end

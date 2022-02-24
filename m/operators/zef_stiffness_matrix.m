function A = zef_stiffness_matrix(basis, tensor)

% The stiffness matrix 𝐴 of a discretized scalar function 𝑢ₕ = ∑ᵢ𝑧ᵢψᵢ, with
% each 𝑧ᵢ being a coordinate and ψᵢ a linear basis function, is defined by
%
%   𝐴[i,j] = ∫[Ω] ∇ψⱼ ⋅ (𝑇∇ψᵢ) d𝑉 ,
%
% where Ω is the entire domain and d𝑉 is a finite volume element from which
% the domain consists of. Here 𝑇 is a tensor, for which the equation
%
%   (𝑇∇u)⋅𝑛⃗ = 0
%
% holds with 𝑢 being the non-discretized scalar function 𝑢 on the boundary ∂Ω
% of the domain and 𝑛⃗ is the surface normal on ∂Ω.

    error("Not yet implemented…")

end

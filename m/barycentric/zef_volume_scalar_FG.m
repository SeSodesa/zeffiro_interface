function M = zef_volume_scalar_FG(nodes,tetra,g_i_ind,scalar_field)

weighting = zef_barycentric_weighting('FG');
M = zef_volume_scalar_D(nodes, tetra, g_i_ind, scalar_field, weighting);

end
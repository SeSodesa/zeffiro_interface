warning('off');
zef.lead_field_type = 2;
zef.imaging_method = 2;
zef.source_ind = [];
zef = zef_process_meshes(zef);
zef = zef_lead_field_matrix(zef);
[zef.L,zef.source_positions,zef.source_directions] = zef_lead_field_filter(zef.L,zef.source_positions,zef.source_directions,zef.lead_field_filter_quantile);
if zef.source_interpolation_on
    zef_source_interpolation;
end
warning('on');

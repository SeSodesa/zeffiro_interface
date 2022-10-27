warning('off');
zef.lead_field_type = 3;
zef.imaging_method = 3;
zef_delete_original_field;
zef_process_meshes;
zef_lead_field_matrix;
[zef.L,zef.source_positions,zef.source_directions] = zef_lead_field_filter(zef.L,zef.source_positions,zef.source_directions,zef.lead_field_filter_quantile);
if zef.source_interpolation_on
    zef_source_interpolation;
end
warning('on');

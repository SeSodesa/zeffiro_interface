zef = zef_add_bounding_box(zef);
zef.exclude_box = 1;
zef.max_surface_face_count = 4;
zef.mesh_smoothing_on = 1;
zef.refinement_on = 1;
zef.refinement_surface_on = 1;
zef.refinement_surface_compartments = [-1 5 6];
zef_mesh_tool;

zef.source_direction_mode = 1; 
zef.inv_sampling_frequency = 2400;
zef.inv_low_cut_frequency = 0;
zef.inv_high_cut_frequency = 0;
evalin('base',['zef.' zef.compartment_tags{1} '_sources = 3;']);
evalin('base',['zef.' zef.compartment_tags{2} '_sources = 2;']);
zef_build_compartment_table;
zeffiro_interface_mesh_tool;
zef_process_meshes;
zef_downsample_surfaces;
zef_process_meshes;
zef_source_interpolation;
zef_ES_update_parameter_values;

zef.ES_search_method                    = get(zef.h_ES_search_method,'Value');
zef.ES_search_type                      = get(zef.h_ES_search_type,'Value');
zef.ES_solver_package                   = zef.h_ES_search_type.Items{ismember(zef.h_ES_search_type.ItemsData,zef.ES_search_type)};
zef.ES_inv_colormap                     = get(zef.h_ES_inv_colormap,'Value');
zef.ES_plot_type                        = get(zef.h_ES_plot_type,'Value');
zef.ES_obj_fun                          = get(zef.h_ES_obj_fun,'Value');
zef.ES_obj_fun_2                        = get(zef.h_ES_obj_fun_2,'Value');
zef.ES_fixed_active_electrodes          = get(zef.h_ES_fixed_active_electrodes,'Value');
zef.ES_threshold_condition              = zef.h_ES_threshold_condition.Value;
zef.ES_algorithm                        = lower(zef.ES_algorithm_list{zef.h_ES_algorithm.Value});

if not(ismember(zef.ES_search_method,zef.h_ES_search_method.ItemsData))
    zef.ES_search_method = zef.h_ES_search_method.ItemsData(1);
end

zef_ES_init_parameter_table;

zef.ES_search_method                    = get(zef.h_ES_search_method,'Value');
zef.ES_search_type                      = get(zef.h_ES_search_type,'Value');
zef.ES_solver_package                   = zef.h_ES_search_type.Items{ismember(zef.h_ES_search_type.ItemsData,zef.ES_search_type)};
zef.ES_inv_colormap                     = get(zef.h_ES_inv_colormap,'Value');
zef.ES_plot_type                        = get(zef.h_ES_plot_type,'Value');
zef.ES_obj_fun                          = get(zef.h_ES_obj_fun,'Value');
zef.ES_obj_fun_2                        = get(zef.h_ES_obj_fun_2,'Value');
zef.ES_fixed_active_electrodes          = get(zef.h_ES_fixed_active_electrodes,'Value');
zef.ES_threshold_condition              = zef.h_ES_threshold_condition.Value;
zef.ES_algorithm                        = lower(zef.ES_algorithm_list{zef.h_ES_algorithm.Value});

if not(ismember(zef.ES_search_method,zef.h_ES_search_method.ItemsData))
    zef.ES_search_method = zef.h_ES_search_method.ItemsData(1);
end
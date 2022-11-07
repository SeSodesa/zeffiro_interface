%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface
function zef = zef_init_ias_roi(zef)

zef = zef_ias_map_estimation_roi_window(zef);

set(zef.h_iasroi_map_estimation,'Position',[ 0.5764    0.2944    0.15    0.50]);

set(zef.h_iasroi_map_estimation,'Name','ZEFFIRO Interface: IAS ROI MAP estimation');
set(findobj(zef.h_iasroi_map_estimation.Children,'-property','FontUnits'),'FontUnits','pixels')
set(findobj(zef.h_iasroi_map_estimation.Children,'-property','FontSize'),'FontSize',zef.font_size);

if not(isfield(zef,'iasroi_rec_source'));
    zef.iasroi_rec_source = [0 0 0 0 0 0 0 3 1];
end;
if not(isfield(zef,'iasroi_roi_mode'));
    zef.iasroi_roi_mode = 3;
end;
if not(isfield(zef,'iasroi_roi_sphere'));
    zef.iasroi_roi_sphere = [0 0 0 15];
end;
if not(isfield(zef,'iasroi_roi_threshold'));
    zef.iasroi_roi_threshold = 0.5;
end;
if not(isfield(zef,'iasroi_hyperprior'));
    zef.iasroi_hyperprior = 2;
end;
if not(isfield(zef,'iasroi_n_map_iterations'));
    zef.iasroi_n_map_iterations = 25;
end;
if not(isfield(zef,'iasroi_pcg_tol'));
    zef.iasroi_pcg_tol = 1e-8;
end;

    zef.iasroi_sampling_frequency = zef.inv_sampling_frequency;
    zef.iasroi_low_cut_frequency = zef.inv_low_cut_frequency;
    zef.iasroi_high_cut_frequency = zef.inv_high_cut_frequency;

if not(isfield(zef,'iasroi_data_segment'));
    zef.iasroi_data_segment = 1;
end;
if not(isfield(zef,'iasroi_normalize_data'));
    zef.iasroi_normalize_data = 1;
end;

    zef.iasroi_time_1 = zef.inv_time_1;
    zef.iasroi_time_2 = zef.inv_time_2;
    zef.iasroi_time_3 = zef.inv_time_3;
    zef.iasroi_number_of_frames = zef.number_of_frames;

zef.iasroi_snr = zef.inv_snr;

set(zef.h_iasroi_rec_source_8 ,'string',num2str(zef.iasroi_rec_source(1,8)));
set(zef.h_iasroi_rec_source_9 ,'value',zef.iasroi_rec_source(1,9));
set(zef.h_iasroi_roi_mode ,'value',zef.iasroi_roi_mode);
set(zef.h_iasroi_roi_sphere_1 ,'string',num2str(zef.iasroi_roi_sphere(:,1)'));
set(zef.h_iasroi_roi_sphere_2 ,'string',num2str(zef.iasroi_roi_sphere(:,2)'));
set(zef.h_iasroi_roi_sphere_3 ,'string',num2str(zef.iasroi_roi_sphere(:,3)'));
set(zef.h_iasroi_roi_sphere_4 ,'string',num2str(zef.iasroi_roi_sphere(:,4)'));
set(zef.h_iasroi_hyperprior ,'value',zef.iasroi_hyperprior);
set(zef.h_iasroi_snr ,'string',num2str(zef.iasroi_snr));
set(zef.h_iasroi_n_map_iterations ,'string',num2str(zef.iasroi_n_map_iterations ));
set(zef.h_iasroi_sampling_frequency ,'string',num2str(zef.iasroi_sampling_frequency));
set(zef.h_iasroi_low_cut_frequency ,'string',num2str(zef.iasroi_low_cut_frequency));
set(zef.h_iasroi_high_cut_frequency ,'string',num2str(zef.iasroi_high_cut_frequency));
set(zef.h_iasroi_normalize_data ,'value',zef.iasroi_normalize_data);
set(zef.h_iasroi_time_1 ,'string',num2str(zef.iasroi_time_1));
set(zef.h_iasroi_time_2 ,'string',num2str(zef.iasroi_time_2));
set(zef.h_iasroi_time_3 ,'string',num2str(zef.iasroi_time_3));
set(zef.h_iasroi_number_of_frames ,'string',num2str(zef.iasroi_number_of_frames));

if get(zef.h_iasroi_roi_mode,'value')==1;
    set(zef.h_iasroi_plot_roi,'enable','on');
    set(zef.h_iasroi_plot_source,'enable','on');
    set(zef.h_iasroi_roi_sphere_1,'enable','on');
    set(zef.h_iasroi_roi_sphere_2,'enable','on');
    set(zef.h_iasroi_roi_sphere_3,'enable','on');
    set(zef.h_iasroi_roi_sphere_4,'enable','on');
    set(zef.h_iasroi_rec_source_8,'enable','on');
    set(zef.h_iasroi_rec_source_9,'enable','on');
    set(zef.h_iasroi_roi_threshold,'enable','off');
else;
    set(zef.h_iasroi_roi_sphere_1,'enable','off');
    set(zef.h_iasroi_roi_sphere_2,'enable','off');
    set(zef.h_iasroi_roi_sphere_3,'enable','off');
    set(zef.h_iasroi_roi_sphere_4,'enable','off');
    set(zef.h_iasroi_rec_source_8,'enable','off');
    set(zef.h_iasroi_rec_source_9,'enable','off');
    set(zef.h_iasroi_plot_roi,'enable','off');
    set(zef.h_iasroi_plot_source,'enable','off');
end;
if get(zef.h_iasroi_roi_mode,'value')==2;
    set(zef.h_iasroi_roi_threshold,'enable','on');
else;
    set(zef.h_iasroi_roi_threshold,'enable','off');
end;

uistack(flipud([zef.h_iasroi_roi_mode; zef.h_iasroi_roi_sphere_1;  zef.h_iasroi_roi_sphere_2;  zef.h_iasroi_roi_sphere_3;
zef.h_iasroi_roi_sphere_4;
zef.h_iasroi_roi_threshold; zef.h_iasroi_hyperprior;
    zef.h_iasroi_snr ; zef.h_iasroi_n_map_iterations ;
    zef.h_iasroi_sampling_frequency ; zef.h_iasroi_low_cut_frequency ;
    zef.h_iasroi_high_cut_frequency ; zef.h_iasroi_time_1 ; zef.h_iasroi_time_2;
    zef.h_iasroi_number_of_frames; zef.h_iasroi_time_3 ; zef.h_iasroi_cancel ;
    zef.h_iasroi_apply; zef.h_iasroi_start ; zef.h_iasroi_plot_source; zef.h_iasroi_rec_source_8;
    zef.h_iasroi_rec_source_9  ]),'top');

end
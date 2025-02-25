%Copyright © 2018- Sampsa Pursiainen & ZI Development Team
%See: https://github.com/sampsapursiainen/zeffiro_interface


if zef.h_segmentation_tool_toggle == 1
    
zef.size_temp = [zef.segmentation_tool_default_position(1) + 0.505*zef.segmentation_tool_default_position(3), ...
                          zef.segmentation_tool_default_position(2),...
                          0.75*0.505*zef.segmentation_tool_default_position(3),...
                          0.75*zef.segmentation_tool_default_position(4)];
                          
else

zef.size_temp = [zef.segmentation_tool_default_position(1) + zef.segmentation_tool_default_position(3), ...
                          zef.segmentation_tool_default_position(2),...
                          0.75*zef.segmentation_tool_default_position(3),...
                          0.6*zef.segmentation_tool_default_position(4)]; 

end

zef.h_zeffiro = figure(...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Units','Pixels',...
'OuterPosition',zef.size_temp,...
'Renderer',get(0,'defaultfigureRenderer'),...
'Visible',zef.use_display,...
'Color',get(0,'defaultfigureColor'),...
'CloseRequestFcn','closereq;',...
'CurrentAxesMode','manual',...
'IntegerHandle','off',...
'NextPlot',get(0,'defaultfigureNextPlot'),...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'DoubleBuffer','off',...
'MenuBar','figure',...
'ToolBar','figure',...
'Name','ZEFFIRO Interface: Figure tool',...
'NumberTitle','off',...
'HandleVisibility','callback',...
'DeleteFcn','zef_reopen_figure;',...
'Tag','figure_tool',...
'UserData',[],...
'WindowStyle',get(0,'defaultfigureWindowStyle'),...
'Resize',get(0,'defaultfigureResize'),...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.99999864 29.69999902],...
'PaperType',get(0,'defaultfigurePaperType'),...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'ScreenPixelsPerInchMode','manual' );


zef.h_zeffiro.Units = 'normalized';

zef.h_zeffiro.ContextMenu = uicontextmenu(zef.h_zeffiro);
uimenu(zef.h_zeffiro.ContextMenu,'Text','Axes pop-up','MenuSelectedFcn','zef_axes_popup;');

addToolbarExplorationButtons(zef.h_zeffiro);

zef.stop_movie = 0;
zef.h_axes1 = uiaxes('Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.05 0.34 0.60 0.60],'FontSize',0.587962962962963,'Tag','axes1');
uicontrol('Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Loop on/count:','HorizontalAlignment','left','Position',[0.68 0.10 0.12 0.03]);

uicontrol('Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','','HorizontalAlignment','left','Position',[0.03 0.95 0.3 0.03],'Tag','time_text');

uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Time:','HorizontalAlignment','left','Position',[0.68 0.90 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Color min:','HorizontalAlignment','left','Position',[0.68 0.85 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Color max:','HorizontalAlignment','left','Position',[0.68 0.80 0.12 0.03]);

%Start controls
zef.h_toggle_controls = uicontrol('Tag','togglecontrolsbutton','UserData',1,'Style','pushbutton','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.68 0.95 0.14 0.05],'String','Toggle controls','Callback','zef_toggle_figure_controls;');
zef.h_toggle_edges = uicontrol('Tag','toggleedgesbutton','UserData',1,'Style','pushbutton','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.83 0.95 0.14 0.05],'String','Toggle edges','Callback','zef_toggle_edges;');

zef.h_colorscale_min_slider = uicontrol('Tag','colorscale_min_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.85 0.17 0.03],'Min',-1,'Max',1,'Value',0,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.colorscale_min_slider = zef_update_colorscale_min; else; zef_update_colorscale_min(gcf); end;');
zef.h_colorscale_min_slider.UserData = zef.colorscale_min_slider;
zef.h_colorscale_max_slider = uicontrol('Tag','colorscale_max_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.80 0.17 0.03],'Min',-1,'Max',1,'Value',0,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.colorscale_max_slider = zef_update_colorscale_max; else; zef_update_colorscale_max(gcf); end;');
zef.h_colorscale_max_slider.UserData = zef.colorscale_max_slider;
zef.h_update_colormap = uicontrol('Tag','colormapselection','Style','popupmenu','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.15 0.08 0.03],'string',zef.colormap_items,'Value',zef.update_colormap,'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_colormap = zef.h_update_colormap.Value; end; zef.h_aux = findobj(get(gcf,''Children''),''Tag'',''axes1'').Colormap; zef.h_aux = zef_colormap(zef.h_update_colormap.Value); zef = rmfield(zef,''h_aux''); zef_update_contrast_and_brightness(gcf);');

zef.h_update_colorscale = uicontrol('Tag','colorscaleselection','Style','popupmenu','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.89 0.15 0.08 0.03],'string',{'Linear','Logarithmic'},'Value',zef.update_colorscale,'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_colorscale = zef_update_colorscale; else; zef_update_colorscale(gcf); end;');
zef.h_update_zoom = uicontrol('Tag','update_zoom_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.75 0.17 0.03],'Min',0.1,'Max',100,'Value',zef.update_zoom,'Sliderstep',[0.001 0.001],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_zoom = zef_update_zoom; else; zef_update_zoom(gcf); end;');
zef.h_update_transparency_reconstruction = uicontrol('Tag','transparency_reconstruction_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.70 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_transparency_reconstruction,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_transparency_reconstruction = zef_update_transparency_reconstruction; else; zef_update_transparency_reconstruction(gcf); end;');
zef.h_update_transparency_surface = uicontrol('Tag','transparency_surface_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.65 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_transparency_surface,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_transparency_surface = zef_update_transparency_surface; else; zef_update_transparency_surface(gcf); end;');
zef.h_update_transparency_sensor = uicontrol('Tag','transparency_sensor_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.60 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_transparency_sensor,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_transparency_sensor = zef_update_transparency_sensor; else; zef_update_transparency_sensor(gcf); end;');
zef.h_update_transparency_cones = uicontrol('Tag','transparency_cones_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.55 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_transparency_cones,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_transparency_cones = zef_update_transparency_cones; else; zef_update_transparency_cones(gcf); end;');
zef.h_update_transparency_additional = uicontrol('Tag','transparency_additional_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.50 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_transparency_additional,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_transparency_additional = zef_update_transparency_additional; else; zef_update_transparency_additional(gcf); end;');
zef.h_update_brightness = uicontrol('Tag','update_brightness_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.45 0.17 0.03],'Min',0,'Max',5,'Value',zef.update_brightness,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); [zef.update_contrast, zef.update_brightness] = zef_update_contrast_and_brightness; else; zef_update_contrast_and_brightness(gcf); end;');
zef.h_update_contrast = uicontrol('Tag','update_contrast_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.40 0.17 0.03],'Min',-1,'Max',1,'Value',zef.update_contrast,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); [zef.update_contrast, zef.update_brightness] = zef_update_contrast_and_brightness; else; zef_update_contrast_and_brightness(gcf); end;');
zef.h_update_ambience = uicontrol('Tag','update_ambience_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.35 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_ambience,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_ambience = zef_update_ambience; else; zef_update_ambience(gcf);end');
zef.h_update_diffusion = uicontrol('Tag','update_diffusion_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.30 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_diffusion,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_diffusion = zef_update_diffusion; else; zef_update_diffusion(gcf);end ');
zef.h_update_specular = uicontrol('Tag','update_specular_slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.25 0.17 0.03],'Min',0,'Max',1,'Value',zef.update_specular,'Sliderstep',[0.01 0.01],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_specular = zef_update_specular; else; zef_update_specular(gcf);end');
zef.h_update_lights = uicontrol('Tag','lightsselection','Style','popupmenu','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.20 0.17 0.03],'string',{'Default (vertical)','Lights off','Add x-lights','Add y-lights','Add z-lights','Add headlight'},'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.update_lights = zef_update_lights; else; zef.update_lights = zef_update_lights(gcf);end');

zef.h_reset_figure_tool_sliders = uicontrol('Style','togglebutton','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.68 0.03 0.065 0.05],'String','Reset','Callback','zef = zef_set_figure_tool_sliders(zef,0);zef.h_reset_figure_tool_sliders.Value=0;');
zef.h_play_movie = uicontrol('Style','pushbutton','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.755 0.03 0.065 0.05],'String','Play','Callback','zef_play_cdata(max(1,double(get(findobj(get(gcf,''Children''),''Tag'',''loop_movie''),''UserData''))*get(findobj(get(gcf,''Children''),''Tag'',''loop_count''),''UserData'')));');
zef.h_stop_movie = uicontrol('Style','togglebutton','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.83 0.03 0.065 0.05],'String','Stop','Callback',@zef_callbackstop);
zef.h_pause_movie = uicontrol('Style','togglebutton','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.905 0.03 0.065 0.05],'String','Pause','Callback',@zef_callbackpause);

zef.h_slider=uicontrol('Tag','slider','Style','slider','Parent',zef.h_zeffiro,'Units','normalized','Position',[0.80 0.90 0.17 0.03],'Min',1e-6,'Max',1,'Value',1e-5,'Sliderstep',[0.01 0.01],'Callback','zef_slidding_callback;');
zef.h_loop_movie = uicontrol('Style','Checkbox','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.89 0.10 0.03 0.03],'Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.loop_movie = get(gcbo,''value''); end;  set(gcbo,''UserData'',get(gcbo,''value''));','HorizontalAlignment','left','Tag','loop_movie');
zef.h_loop_movie_count = uicontrol('Style','Edit','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.94 0.10 0.03 0.03],'String','Loop visualization','Callback','if isequal(get(gca,''Parent''), zef.h_zeffiro); zef.loop_movie_count = str2num(get(gcbo,''string'')); end; set(gcbo,''UserData'',str2num(get(gcbo,''string'')));','HorizontalAlignment','right','Tag','loop_count');
%End controls

uicontrol('Tag','colormapselectiontext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Colormap:','HorizontalAlignment','left','Position',[0.68 0.15 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Distance:','HorizontalAlignment','left','Position',[0.68 0.75 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Transp. rec.:','HorizontalAlignment','left','Position',[0.68 0.70 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Transp. surf.:','HorizontalAlignment','left','Position',[0.68 0.65 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Transp. sens.:','HorizontalAlignment','left','Position',[0.68 0.60 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Transp. cones:','HorizontalAlignment','left','Position',[0.68 0.55 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Transp. add.:','HorizontalAlignment','left','Position',[0.68 0.50 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Brightness:','HorizontalAlignment','left','Position',[0.68 0.45 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Contrast:','HorizontalAlignment','left','Position',[0.68 0.40 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Ambience:','HorizontalAlignment','left','Position',[0.68 0.35 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Diffusion:','HorizontalAlignment','left','Position',[0.68 0.30 0.12 0.03]);
uicontrol('Tag','slidertext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Specular exp.:','HorizontalAlignment','left','Position',[0.68 0.25 0.12 0.03]);
uicontrol('Tag','lightstext','Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Lights:','HorizontalAlignment','left','Position',[0.68 0.20 0.12 0.03]);

set(zef.h_loop_movie_count,'string',num2str(zef.loop_movie_count));
imagesc(zef.h_axes1,flipud(imread('zeffiro_interface_compass.png')));
set(zef.h_axes1,'YDir','normal');
axis(zef.h_axes1,'auto');
zef.h_axes1.YLim = [0 1961];
zef.h_axes1.XLim = [0 2592];


%***********************

zef.h_compartment_visible_color = uicontrol('Style','ListBox','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.03 0.03 0.20 0.20],'Tag','compartment_visible_color');
zef.h_sensor_visible_color = uicontrol('Style','ListBox','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.24 0.03 0.20 0.20],'Tag','sensor_visible_color');
zef.h_system_information = uicontrol('Style','ListBox','Parent',zef.h_zeffiro,'visible','on','Units','normalized','Position',[0.45 0.03 0.20 0.20],'Tag','system_information');
set(zef.h_compartment_visible_color,'ButtonDownFcn','zef_set_compartment_color; zef_update;');
set(zef.h_sensor_visible_color,'ButtonDownFcn','zef_set_sensor_color; zef_update;');

uicontrol(...
'Parent',zef.h_zeffiro,...
'Units','normalized',...
'FontUnits','normalized',...
'HorizontalAlignment','left',...
'String','Compartments:',...
'Style','text',...
'Position',[0.03 0.23 0.20 0.03],...
'Children',[],'FontSize',0.461256944444443);

uicontrol(...
'Parent',zef.h_zeffiro,...
'Units','normalized',...
'FontUnits','normalized',...
'HorizontalAlignment','left',...
'String','Sensors:',...
'Style','text',...
'Position',[0.24 0.23 0.20 0.03],...
'Children',[],'FontSize',0.461256944444443);

uicontrol(...
'Parent',zef.h_zeffiro,...
'Units','normalized',...
'FontUnits','normalized',...
'HorizontalAlignment','left',...
'String','Details:',...
'Style','text',...
'Position',[0.45 0.23 0.20 0.03],...
'Children',[],'FontSize',0.461256944444443);

uicontrol('Style','text','Parent',zef.h_zeffiro,'Units','normalized','String','Copyright © 2018- Sampsa Pursiainen & ZI Development Team. See: https://github.com/sampsapursiainen/zeffiro_interface','HorizontalAlignment','left','Position',[0.03 0.005 0.6 0.02],'Tag','copyright_text');

zef = zef_update_fig_details(zef);

zef.h = get(zef.h_zeffiro,'Children');
zef.h = zef.h(isprop(zef.h,'Units'));
set(zef.h,'Units','pixels');
zef = rmfield(zef,'h');

set(zef.h_zeffiro,'handlevisibility','on');
set(zef.h_zeffiro,'WindowButtonDownFcn','zef.h_zeffiro = zef.h_zeffiro; zef.h_axes1 = findobj(get(zef.h_zeffiro,''Children''),''Tag'',''axes1'');')

zef.h_zeffiro.GraphicsSmoothing = 'off';

set(findobj(get(gcf,'Children'),'Tag','loop_count'),'UserData',str2num(get(findobj(get(gcf,'Children'),'Tag','loop_count'),'String')));
set(findobj(get(gcf,'Children'),'Tag','loop_movie'),'UserData',get(findobj(get(gcf,'Children'),'Tag','loop_count'),'Value'));

set(zef.h_zeffiro,'paperposition',[0 0 zef.snapshot_horizontal_resolution/200 zef.snapshot_vertical_resolution/200]);
set(zef.h_zeffiro,'papersize',[zef.snapshot_vertical_resolution/200 zef.snapshot_horizontal_resolution/200]);

if zef.clear_axes1
zef.h_colorbar = findobj(zef.h_zeffiro,'tag','Colorbar');
if not(isempty(zef.h_colorbar))
colorbar(zef.h_colorbar,'delete');
end
else
zef.clear_axes1 = 1;
end

if isfield(zef,'zeffiro_current_size')
if not(iscell(zef.zeffiro_current_size))
zef = rmfield(zef,'zeffiro_current_size');
end
end

set(zef.h_zeffiro,'Name',[get(zef.h_zeffiro,'Name') ' ' num2str(zef_fig_num)]);
set(zef.h_zeffiro,'AutoResizeChildren','off');

if not(ismember('ZefFig',properties(zef.h_zeffiro)))
addprop(zef.h_zeffiro,'ZefFig');
end
set(zef.h_zeffiro,'ZefFig',zef_fig_num);


set(findobj(zef.h_zeffiro.Children,'-property','FontUnits'),'FontUnits','pixels')
set(findobj(zef.h_zeffiro.Children,'-property','FontSize'),'FontSize',zef.font_size);
set(findobj(zef.h_zeffiro.Children,'Tag','copyright_text'),'FontSize',0.5*zef.font_size);

zef.h_zeffiro.SizeChangedFcn = '';
zef.h_zeffiro.Units = 'normalized';
zef_set_size_change_function(zef.h_zeffiro,1,[],['{''rightColorbar'',''leftColorbar'',''legend''}']);
zef = rmfield(zef,'size_temp');


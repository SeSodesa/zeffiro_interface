function zef = zeffiro_interface(args)
%This fuction starts Zeffiro Interface. It can be run with a variable
%number of arguments, which can be called as a list of name-value pairs as
%follows:
%
%zeffiro_interface('property name 1','property value 1','property name 2','property value 2');
%
%This will enable running Zeffiro with or without a display and performing
%different operations. The list of properties (and their values) is the
%following:
%
%Property: 'restart'                     Value: none
%Property: 'start_mode'                  Value: 'display' or 'nodisplay'
%Propertu: 'open_project'                Value: <project file name>,
%Property: 'import_to_new_project'       Value: <file name>,
%Property: 'import_to_existing_project'  Value: <file name>,
%Property: 'save_project'                Value: <file name>,
%Property: 'export_fem_mesh'             Value: <file name>,
%Property: 'open_figure'                 Value: <file name>,
%Property: 'open_figure_folder'          Value: <file name>,
%Property: 'run_script'                  Value: <file name>,
%Property: 'exit_zeffiro'                Value: none
%Property: 'quit_matlab'                 Vaule: none
%Property: 'use_github'                  Value: 1 (yes) or 0 (no)
%Property: 'use_gpu'                     Value: 1 (yes) or 0 (no)
%Property: 'use_gpu_graphic'             Value: 1 (yes) or 0 (no)
%Property: 'gpu_num'                     Value: <gpu device number>
%Property: 'use_display'                 Value: 1 (yes) or 0 (no)
%Property: 'parallel_processes'          Value: <parallel pool size>
%Property: 'verbose_mode'                Value: 1 (yes) or 0 (no)
%Property: 'use_waitbar'                 Value: 1 (yes) or 0 (no)
%Property: 'use_log'                     Value: 1 (yes) or 0 (no)
%Property: 'log_file_name'               Value: <log file name>

    arguments

        args.restart (1,1) logical = false;

        args.start_mode (1,1) string { mustBeMember(args.start_mode, ["display", "nodisplay", "default"]) } = "default";

        args.open_project (1,1) string = "";

        args.import_to_new_project (1,1) string = "";

        args.import_to_existing_project (1,1) string = "";

        args.save_project (1,1) string = "";

        args.export_fem_mesh (1,1) string = "";

        args.open_figure (1,1) string = "";

        args.open_figure_folder (1,1) string = "";

        args.run_script (1,1) string = "";

        args.exit_zeffiro (1,1) logical = false;

        args.quit_matlab (1,1) logical = false;

        args.use_github (1,1) logical = false;

        args.use_gpu (1,1) logical = false;

        args.use_gpu_graphic (1,1) logical = false;

        args.gpu_num (1,1) double { mustBeInteger, mustBeNonnegative } = 0;

        args.use_display (1,1) logical = false;

        args.parallel_processes (1,1) double {mustBePositive, mustBeInteger} = 1;

        args.verbose_mode (1,1) logical = false;

        args.use_waibar (1,1) logical = false;

        args.use_log (1,1) logical = false;

        args.log_file_name (1,1) string = "";

    end

    %% Set zef fields based on name–value arguments.

    % TODO: add code below to handle different values of the arguments.

    zef = struct;

    zef.zeffiro_restart = args.restart;

    zef.start_mode = args.start_mode;

    zef.open_project = args.open_project;

    zef.import_to_new_project = args.import_to_new_project;

    zef.import_to_existing_project = args.import_to_existing_project;

    zef.save_project = args.save_project;

    zef.export_fem_mesh = args.export_fem_mesh;

    zef.open_figure = args.open_figure;

    zef.open_figure_folder = args.open_figure_folder;

    zef.run_script = args.run_script;

    zef.exit_zeffiro = args.exit_zeffiro;

    zef.quit_matlab = args.quit_matlab;

    zef.use_github = args.use_github;

    zef.use_gpu = args.use_gpu;

    zef.use_gpu_graphic = args.use_gpu_graphic;

    zef.gpu_num = args.gpu_num;

    zef.use_display = args.use_display;

    zef.parallel_processes = args.parallel_processes;

    zef.verbose_mode = args.verbose_mode;

    zef.use_waibar = args.use_waibar;

    zef.use_log = args.use_log;

    zef.log_file_name = args.log_file_name;

    %% Then do initial preparations

    program_path_aux = mfilename("fullpath");

    [program_path, ~] = fileparts(program_path_aux);

    run(program_path + filesep + "m/zef_close_all.m"]);

    zef.zeffiro_task_id = 0;

    zef.zeffiro_restart_time = now;

    zef.zeffiro_restart = zeffiro_restart;

    zef.program_path = program_path;

    zef.code_path = zef.program_path + filesep + "m";

    zef.cluster_path =  zef.program_path + filesep + "cluster";

    addpath(zef.program_path);
    addpath(zef.code_path);
    addpath(zef.program_path);
    addpath(genpath(zef.code_path));
    addpath(genpath(zef.cluster_path));
    addpath(genpath(zef.program_path + filesep + "mlapp"));
    addpath(genpath(zef.program_path + filesep + "fig"));
    addpath(genpath(zef.program_path + filesep + "plugins"));
    addpath(genpath(zef.program_path + filesep + "profile"));
    addpath(genpath(zef.program_path + filesep + "scripts"));

    addpath(zef.program_path + filesep + "external");

    zef.start_mode = "default";

    if exist("zef_start_config.m","file")
        eval("zef_start_config");
    end

    %% Do things based on input arguments.

    % Prevent starting of Zeffiro, if there is an existing value of zef.

    if not(zef.zeffiro_restart) && evalin('base','exist(''zef'',''var'');')

        error('It looks like that another instance of Zeffiro interface is already open. To enable this script, close Zeffiro Interface by command ''zef_close_all'' or clear zef by command ''clear zef''.')

    end

    % Open new project if given.

    if not(open_project == "")

        open_project_file = open_project;

        [file_path, file_1, file_2] = fileparts(open_project_file);

        file_path = file_path + filesep;

        if isempty(file_path)
            file_path = "./data/";
        end

        if isempty(file_2)
            file_2 = ".mat";
        end

        zef.file_path = [file_path];

        zef.file = [file_1 file_2];

        zef = zef_load(zef,zef.file,zef.file_path);

    end % if

    % Importing given file contents to a new project.

    if not(import_to_new_project == "")

        import_segmentation_file = import_to_new_project;

        [file_path, file_1, file_2] = fileparts(import_segmentation_file);

        file_path = [file_path filesep];

        if isempty(file_path)
            file_path = "./data/";
        end

        if isempty(file_2)
            file_2 = ".mat";
        end

        zef.new_empty_project = 1;
        zef_start_new_project;
        zef.file_path = [file_path];
        zef.file = [file_1 file_2];
        zef = zef_import_segmentation(zef);
        zef = zef_build_compartment_table(zef);

    end % if

    % Import given file contents into an existing project.

    if not(import_to_existing_project == "")

        import_segmentation_file = import_to_existing_project;

        [file_path, file_1, file_2] = fileparts(import_segmentation_file);

        file_path = file_path + filesep;

        if isempty(file_path)
            file_path = "./data/";
        end

        if isempty(file_2)
            file_2 = ".mat";
        end

        zef.file_path = file_path;

        zef.file = file_1 + file_2;

        zef.new_empty_project = 0;

        zef = zef_import_segmentation(zef);

        zef = zef_build_compartment_table(zef);

    end % if

    % Choose GPU device, if available.

    zef.gpu_count = gpuDeviceCount;

    if zef.gpu_count > 0 && zef.use_gpu

        try

            gpuDevice(zef.gpu_num);

        catch

            warning("Tried using GPU with index " + zef.gpu_num + " but no such device was found. Starting without GPU...");

        end

    end % if

    if not(export_fem_mesh == "")

        export_fem_mesh_file = export_fem_mesh;

        [file_path, file_1, file_2] = fileparts(export_fem_mesh_file );

        file_path = file_path + filesep;

        if isempty(file_path)
            file_path = './data/';
        end

        if isempty(file_2)
            file_2 = '.mat';
        end

        zef.file_path = file_path;

        zef.file = file_1 + file_2;

        zef.save_switch = 1;

        zef_export_fem_mesh_as(zef);

        option_counter = option_counter + 2;

    end % if

    if not(open_figure == "")

        open_figure_file = open_figure;

        if not(iscell(open_figure_file))

            open_figure_file_aux = open_figure_file;

            open_figure_file = cell(0);

            open_figure_file{1} = open_figure_file_aux;

        end

        for i = 1 : length(open_figure_file)

            [file_path, file_1, file_2] = fileparts(open_figure_file{i});

            file_path = file_path + filesep;

            if isempty(file_path)
                file_path = './fig/';
            end

            if isempty(file_2)
                file_2 = '.fig';
            end

            zef.file_path = file_path;

            zef.file = file_1 + file_2;

            zef.save_switch = 1;

            zef = zef_import_figure(zef);

        end % for

    end % if

    if not(open_figure_folder == "")

        file_path = open_figure_folder;

        dir_aux = dir(fullfile(zef_data.program_path,file_path));

        for i = 3 : length(dir_aux)

            [~, file_1, file_2] = fileparts(dir_aux(i).name);

            if isequal(file_2, ".fig")

                zef.file_path = file_path;

                zef.file = file_1 + file_2;

                zef.save_switch = 1;

                zef = zef_import_figure(zef);

            end % if

        end % for

    end % if

    if not(run_script == "")

        run_script_name = run_script;

        if not(iscell(run_script_name))

           run_script_name_aux = run_script_name;

           run_script_name = cell(0);

           run_script_name{1} = run_script_name_aux;

        end

        for i = 1 : length(run_script_name)
            eval(run_script_name{i});
        end

    end % if

    if zef.exit_zeffiro
        zef_close_all;
    end

    if zef.quit_matlab
        quit force;
    end

%% Old code below. TODO: remove this along with the varargin.

warning off;
option_counter = 1;
zeffiro_restart = 0;

if not(isempty(varargin))

    if isequal(varargin{1},'restart')

        zeffiro_restart = 1;
        option_counter = option_counter + 1;

    end

end

if nargout == 0

    if isequal(zeffiro_restart,0)

        if evalin("base',"exist('zef','var');")
            error("It looks like that another instance of Zeffiro interface is already open. To enable this script, close Zeffiro Interface by command ''zef_close_all'' or clear zef by command ''clear zef''.')
        end

    end

end

    program_path_aux = mfilename("fullpath");
    [program_path, ~] = fileparts(program_path_aux);
    run(program_path + filesep "m/zef_close_all.m"]);
    zef = struct;

    zef.zeffiro_task_id = 0;
    zef.zeffiro_restart_time = now;
    zef.zeffiro_restart = zeffiro_restart;
    zef.program_path = program_path;
    zef.code_path = zef.program_path + filesep + "m"];
    zef.cluster_path =  zef.program_path + filesep + "cluster"];

    addpath(zef.program_path);
    addpath(zef.code_path);
    addpath(zef.program_path);
    addpath(genpath([zef.code_path]));
    addpath(genpath([zef.cluster_path]));
    addpath(genpath([zef.program_path filesep "mlapp"]));
    addpath(genpath([zef.program_path filesep "fig"]));
    addpath(genpath([zef.program_path filesep "plugins"]));
    addpath(genpath([zef.program_path filesep "profile"]));
    addpath(genpath([zef.program_path filesep "scripts"]));
    addpath(zef.program_path + filesep + "external");

    zef.start_mode = "default";

    if exist("zef_start_config.m","file")
      eval("zef_start_config");
    end

    if not(isempty(varargin))

        start_mode = 'display';

     while option_counter <= length(varargin)

            if ischar(varargin{option_counter})
            if ismember(lower(varargin{option_counter}),{lower('display'),lower('nodisplay')})
                start_mode = lower(varargin{option_counter});
                option_counter = option_counter + 1;
            elseif  ismember(lower(varargin{option_counter}),{'start_mode'})
                start_mode = lower(varargin{option_counter+1});
                option_counter = option_counter + 2;
            elseif ismember(varargin{option_counter},lower('profile_name'))
                zef.ini_cell_mod = {'Profile name',varargin{option_counter+1},'profile_name','string'};
                option_counter = option_counter + 2;

            elseif isequal(varargin{option_counter},lower('use_github'))

                use_github = varargin{option_counter+1};

                option_counter = option_counter + 2;

                elseif isequal(varargin{option_counter},lower('use_gpu'))

                use_gpu = varargin{option_counter+1};

                option_counter = option_counter + 2;

                elseif isequal(varargin{option_counter},lower('use_gpu_graphic'))

                use_gpu_graphic = varargin{option_counter+1};

                option_counter = option_counter + 2;

                elseif isequal(varargin{option_counter},lower('gpu_num'))

                gpu_num = varargin{option_counter+1};

                option_counter = option_counter + 2;

                elseif isequal(varargin{option_counter},lower('use_display'))

                use_display = varargin{option_counter+1};

                option_counter = option_counter + 2;

               elseif isequal(varargin{option_counter},lower('parallel_processes'))

               parallel_processes = varargin{option_counter+1};

                option_counter = option_counter + 2;

               elseif isequal(varargin{option_counter},lower('verbose_mode'))

               verbose_mode = varargin{option_counter+1};

                option_counter = option_counter + 2;

                elseif isequal(varargin{option_counter},lower('use_log'))

               use_log = varargin{option_counter+1};

                option_counter = option_counter + 2;

                  elseif isequal(varargin{option_counter},lower('log_file_name'))

               log_file_name = varargin{option_counter+1};

                option_counter = option_counter + 2;

            elseif isequal(varargin{option_counter},lower('use_waitbar'))

               use_waitbar = varargin{option_counter+1};

                option_counter = option_counter + 2;

            else
                option_counter = option_counter + 1;
            end
            else
               option_counter = option_counter + 1;
            end
     end

        zef.start_mode = 'nodisplay';
        zef = zef_start(zef);


        if isequal(zef.zeffiro_restart,0) && isequal(exist([zef.program_path filesep 'data' filesep 'default_project.mat']),2)
        zef = zef_load(zef,'default_project.mat',[zef.program_path filesep 'data' filesep]);
        end

        if exist('use_github','var')
         zef.use_github = use_github;
        end
         if exist('use_gpu','var')
         zef.use_gpu = use_gpu;
          end
         if exist('use_gpu_graphic','var')
           zef.use_gpu_graphic = use_gpu_graphic;
          end
             if exist('gpu_num','var')
           zef.gpu_num = gpu_num;
       end
         if exist('parallel_processes','var')
         zef.parallel_processes = parallel_processes;
         end
        if exist('verbose_mode','var')
         zef.zeffiro_verbose_mode = verbose_mode;
        end
        if exist('use_log','var')
         zef.use_log = use_log;
        end
        if exist('log_file_name','var')
         zef.zeffiro_log_file_name = log_file_name;
        end
         if exist('use_waitbar','var')
         zef.use_waitbar = use_waitbar;
         end
           if exist('use_display','var')
               zef.use_display =  use_display;
           end

         zef = zef_start_log(zef);

       if and(zef.gpu_count > 0, zef.use_gpu)
      gpuDevice(zef.gpu_num);
       end

        option_counter = 1;

        while option_counter <= length(varargin)

            if ischar(varargin{option_counter})

            if isequal(varargin{option_counter},lower('open_project'))

                open_project_file = varargin{option_counter+1};
                [file_path, file_1, file_2] = fileparts(open_project_file);
                file_path = [file_path filesep];

                if isempty(file_path)
                    file_path = './data/';
                end

                if isempty(file_2)
                    file_2 = '.mat';
                end

                zef.file_path = [file_path];
                zef.file = [file_1 file_2];
                zef = zef_load(zef,zef.file,zef.file_path);
                option_counter = option_counter + 2;

               if exist('use_github','var')
                zef.use_github = use_github;
               end
               if exist('use_gpu','var')
                zef.use_gpu = use_gpu;
               end
               if exist('use_gpu_graphic','var')
                   zef.use_gpu_graphic = use_gpu_graphic;
               end
                if exist('gpu_num','var')
                zef.gpu_num = gpu_num;
                end
                 if exist('parallel_processes','var')
                 zef.parallel_processes = parallel_processes;
                 end
                 if exist('verbose_mode','var')
                 zef.zeffiro_verbose_mode = verbose_mode;
                 end
                 if exist('use_log','var')
                 zef.use_log = use_log;
                end
                if exist('log_file_name','var')
                 zef.zeffiro_log_file_name = log_file_name;
                end
                if exist('use_waitbar','var')
                 zef.use_waitbar = use_waitbar;
                 end
                   if exist('use_display','var')
                       zef.use_display =  use_display;
                   end
                   if and(zef.gpu_count > 0, zef.use_gpu)
                  gpuDevice(zef.gpu_num);
                   end

            elseif isequal(varargin{option_counter},lower('import_to_new_project'))

                import_segmentation_file = varargin{option_counter+1};
                [file_path, file_1, file_2] = fileparts(import_segmentation_file);
                file_path = [file_path filesep];

                if isempty(file_path)
                    file_path = './data/';
                end

                if isempty(file_2)
                    file_2 = '.mat';
                end

                zef.new_empty_project = 1;
                zef_start_new_project;
                zef.file_path = [file_path];
                zef.file = [file_1 file_2];
                zef = zef_import_segmentation(zef);
                zef = zef_build_compartment_table(zef);
                option_counter = option_counter + 2;

            elseif isequal(varargin{option_counter},lower('import_to_existing_project'))

                import_segmentation_file = varargin{option_counter+1};
                [file_path, file_1, file_2] = fileparts(import_segmentation_file);
                file_path = [file_path filesep];

                if isempty(file_path)
                    file_path = './data/';
                end

                if isempty(file_2)
                    file_2 = '.mat';
                end

                zef.file_path = [file_path];
                zef.file = [file_1 file_2];
                zef.new_empty_project = 0;
                zef = zef_import_segmentation(zef);
                zef = zef_build_compartment_table(zef);
                option_counter = option_counter + 2;


            elseif isequal(varargin{option_counter},lower('save_project'))

                save_project_file = varargin{option_counter+1};
                [file_path, file_1, file_2] = fileparts(save_project_file);
                file_path = [file_path filesep];

                if isempty(file_path)
                    file_path = './data/';
                end

                if isempty(file_2)
                    file_2 = '.mat';
                end

                zef.file_path = [file_path];
                zef.file = [file_1 file_2];
                zef.save_switch = 1;
                zef = zef_save(zef);
                option_counter = option_counter + 2;

            elseif isequal(varargin{option_counter},lower('export_fem_mesh'))

                export_fem_mesh_file = varargin{option_counter+1};
                [file_path, file_1, file_2] = fileparts(export_fem_mesh_file );
                file_path = [file_path filesep];

                if isempty(file_path)
                    file_path = './data/';
                end

                if isempty(file_2)
                    file_2 = '.mat';
                end

                zef.file_path = [file_path];
                zef.file = [file_1 file_2];
                zef.save_switch = 1;
                zef_export_fem_mesh_as(zef);
                option_counter = option_counter + 2;

            elseif ismember(varargin{option_counter},lower('open_figure'))

                open_figure_file = varargin{option_counter+1};

                if not(iscell(open_figure_file))
                    open_figure_file_aux = open_figure_file;
                    open_figure_file = cell(0);
                    open_figure_file{1} = open_figure_file_aux;
                end

                for i = 1 : length(open_figure_file)

                    [file_path, file_1, file_2] = fileparts(open_figure_file{i});
                    file_path = [file_path filesep];

                    if isempty(file_path)
                        file_path = './fig/';
                    end

                    if isempty(file_2)
                        file_2 = '.fig';
                    end

                    zef.file_path = [file_path];
                    zef.file = [file_1 file_2];
                    zef.save_switch = 1;
                    zef = zef_import_figure(zef);
                    option_counter = option_counter + 2;
                end

            elseif ismember(varargin{option_counter},lower('open_figure_folder'))

                file_path = varargin{option_counter+1};
                dir_aux = dir(fullfile(zef_data.program_path,file_path));

                for i = 3 : length(dir_aux)

                    [~,file_1,file_2] = fileparts(dir_aux(i).name);

                    if isequal(file_2,'.fig')
                        zef.file_path = [file_path];
                        zef.file = [file_1 file_2];
                        zef.save_switch = 1;
                        zef = zef_import_figure(zef);
                        option_counter = option_counter + 2;
                    end
                end

                option_counter = option_counter + 2;

            elseif ismember(varargin{option_counter},lower('run_script'))

                run_script_name = varargin{option_counter+1};

                if not(iscell(run_script_name))
                   run_script_name_aux = run_script_name;
                   run_script_name = cell(0);
                   run_script_name{1} = run_script_name_aux;
                end

                for i = 1 : length(run_script_name)
                    eval(run_script_name{i});
                end

                option_counter = option_counter + 2;

            elseif ismember(varargin{option_counter},lower('exit_zeffiro'))
                zef_close_all;
                option_counter = option_counter + 1;
            elseif ismember(varargin{option_counter},lower('quit_matlab'))
                quit force;
                option_counter = option_counter + 1;
            else
                option_counter = option_counter + 1;
            end
            else
                 option_counter = option_counter + 1;
            end
        end


        if exist('zef','var')
        if isfield(zef,'h_zeffiro_window_main')
            if isvalid(zef.h_zeffiro_window_main)
                if exist('use_display','var')
                    start_mode = 'display';
                end
                if ismember(start_mode,'display')
                    zef.start_mode = start_mode;
                    zef.h_zeffiro.Visible = 1;
                    zef.h_zeffiro_window_main.Visible = 1;
                    zef.h_mesh_tool.Visible = 1;
                    zef.h_mesh_visualization_tool.Visible = 1;
                    zef.h_zeffiro_menu.Visible = 1;
                    zef.use_display = 1;


       if and(zef.gpu_count > 0, zef.use_gpu)
      gpuDevice(zef.gpu_num);
       end

                end
            end
        end


        if exist('use_github','var')
         zef.use_github = use_github;
          end
         if exist('use_gpu','var')
         zef.use_gpu = use_gpu;
          end
         if exist('use_gpu_graphic','var')
           zef.use_gpu_graphic = use_gpu_graphic;
          end
             if exist('gpu_num','var')
           zef.gpu_num = gpu_num;
       end
         if exist('parallel_processes','var')
         zef.parallel_processes = parallel_processes;
         end
         if exist('verbose_mode','var')
         zef.zeffiro_verbose_mode = verbose_mode;
         end
         if exist('use_log','var')
         zef.use_log = use_log;
        end
        if exist('log_file_name','var')
         zef.zeffiro_log_file_name = log_file_name;
        end
    if exist('use_waitbar','var')
         zef.use_waitbar = use_waitbar;
         end
           if exist('use_display','var')
               zef.use_display =  use_display;
           end

    zef.zeffiro_restart = 0;
        end

            else
        zef = zef_start(zef);
        zef = zef_start_log(zef);
    if isequal(zef.zeffiro_restart,0) && exist([zef.program_path filesep 'data' filesep 'default_project.mat'],'file')
        zef = zef_load(zef,'default_project.mat',[zef.program_path filesep 'data' filesep]);
    end
    end

        if nargout == 0
    if    exist('zef','var')
        assignin('base','zef',zef);
    end
    end


warning on;
end

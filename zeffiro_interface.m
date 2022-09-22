function zef = zeffiro_interface(args)

%
% zeffiro_interface
%
% This fuction starts Zeffiro Interface. It can be run with a variable number
% of key–value arguments, which can be called as a list of name-value pairs as
% follows:
%
%   zef = zeffiro_interface('name_1', value_1, 'name_2', value_2, …);
%
% Matlab versions newer or equal to r2021a may also use the syntax
%
%   zef = zeffiro_interafce(name_1=value_1, name_2=value_2, …)
%
% Either way, this enables running Zeffiro with or without a display and
% performing different operations. The list of names (and their values) is the
% following:
%
%   name                            value
%
%   'restart'                       none
%   'start_mode'                    'display' or 'nodisplay'
%   'open_project'                  <project file name>,
%   'import_to_new_project'         <file name>,
%   'import_to_existing_project'    <file name>,
%   'save_project'                  <file name>,
%   'export_fem_mesh'               <file name>,
%   'open_figure'                   <file name>,
%   'open_figure_folder'            <file name>,
%   'run_script'                    <file name>,
%   'exit_zeffiro'                  none
%   'quit_matlab'                   none
%   'use_github'                    1 (yes) or 0 (no)
%   'use_gpu'                       1 (yes) or 0 (no)
%   'use_gpu_graphic'               1 (yes) or 0 (no)
%   'gpu_num'                       <gpu device number>
%   'use_display'                   1 (yes) or 0 (no)
%   'parallel_processes'            <parallel pool size>
%   'verbose_mode'                  1 (yes) or 0 (no)
%   'use_waitbar'                   1 (yes) or 0 (no)
%   'use_log'                       1 (yes) or 0 (no)
%   'log_file_name'                 <log file name>
%
%   NOTE: the value behind the name run_script is run using the Matlab
%   function eval, meaning one should be absolutely sure that the script
%   contents come from a trusted source.
%

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

        args.use_waitbar (1,1) logical = false;

        args.use_log (1,1) logical = false;

        args.log_file_name (1,1) string = "";

    end

    % Prevent starting of Zeffiro, if there is an existing value of zef.

    if not(args.restart) && evalin("base","exist('zef', 'var');")

        error("It looks like that another instance of Zeffiro interface is already open. To enable this script, close Zeffiro Interface by command 'zef_close_all' or clear zef by command 'clear zef'.")

    end

    %% Set zef fields based on name–value arguments.

    zef = struct;

    zef.zeffiro_restart = args.restart;

    zef.start_mode = args.start_mode;

    zef.use_github = args.use_github;

    zef.use_gpu = args.use_gpu;

    zef.use_gpu_graphic = args.use_gpu_graphic;

    zef.gpu_num = args.gpu_num;

    zef.use_display = args.use_display;

    zef.parallel_processes = args.parallel_processes;

    zef.verbose_mode = args.verbose_mode;

    zef.use_waitbar = args.use_waitbar;

    zef.use_log = args.use_log;

    zef.log_file_name = args.log_file_name;

    %% Then do initial preparations like path building and additions.

    program_path = mfilename("fullpath");

    [program_path, ~] = fileparts(program_path);

    program_path = string(program_path);

    code_path = program_path + filesep + "m";

    % TODO: should this be run here?
    %
    % run(code_path + filesep + "zef_close_all.m");

    zef.program_path = program_path;

    zef.code_path = code_path;

    zef.data_path = zef.program_path + filesep + "data";

    zef.external_path = zef.program_path + filesep + "external";

    zef.zeffiro_task_id = 0;

    zef.zeffiro_restart_time = now;

    zef.cluster_path =  zef.program_path + filesep + "cluster";

    addpath(zef.program_path);
    addpath(zef.code_path);
    addpath(genpath(zef.code_path));
    addpath(genpath(zef.cluster_path));

    addpath(genpath(zef.program_path + filesep + "mlapp"));
    addpath(genpath(zef.program_path + filesep + "fig"));
    addpath(genpath(zef.program_path + filesep + "plugins"));
    addpath(genpath(zef.program_path + filesep + "profile"));
    addpath(genpath(zef.program_path + filesep + "scripts"));

    addpath(zef.external_path);

    if not(zef.zeffiro_restart)

        addpath(zef.external_path + filesep + "SDPT3/");
        addpath(zef.external_path + filesep + "SeDuMi/");
        addpath(zef.external_path + filesep + "CVX/");

        % TODO: does not work.
        %
        % evalc("cvx_startup");

    end

    zef = zef_start(zef);

    if not(zef.zeffiro_restart) && isfile(zef.data_path + filesep + "default_project.mat")

        zef = zef_load(zef, "default_project.mat", zef.data_path + filesep);

    end

    zef = zef_start_log(zef);

    if isfield(zef, "h_zeffiro_window_main") ...
    && isvalid(zef.h_zeffiro_window_main) ...
    && zef.start_mode == "display"

        zef.start_mode = start_mode;
        zef.h_zeffiro.Visible = 1;
        zef.h_zeffiro_window_main.Visible = 1;
        zef.h_mesh_tool.Visible = 1;
        zef.h_mesh_visualization_tool.Visible = 1;
        zef.h_zeffiro_menu.Visible = 1;
        zef.use_display = 1;

    end

    %% Finally, do the things specified by the input arguments.

    % Choose GPU device, if available.

    zef.gpu_count = gpuDeviceCount;

    if zef.gpu_count > 0 && zef.use_gpu

        try

            gpuDevice(zef.gpu_num);

        catch

            warning("Tried using GPU with index " + zef.gpu_num + " but no such device was found. Starting without GPU...");

        end

    end % if

    % Open new project if given.

    if not(args.open_project == "")

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

    % Import given file contents to a new project.

    if not(args.import_to_new_project == "")

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

        zef.file_path = file_path;

        zef.file = file_1 + file_2;

        zef = zef_import_segmentation(zef);

        zef = zef_build_compartment_table(zef);

    end % if

    % Import given file contents into an existing project.

    if not(args.import_to_existing_project == "")

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

    % Export FE mesh to a given path.

    if not(args.export_fem_mesh == "")

        export_fem_mesh_file = export_fem_mesh;

        [file_path, file_1, file_2] = fileparts(export_fem_mesh_file );

        file_path = file_path + filesep;

        if isempty(file_path)

            file_path = "./data/";

        end

        if isempty(file_2)

            file_2 = ".mat";

        end

        zef.file_path = file_path;

        zef.file = file_1 + file_2;

        zef.save_switch = 1;

        zef_export_fem_mesh_as(zef);

        option_counter = option_counter + 2;

    end % if

    % Open figure in a given path.

    if not(args.open_figure == "")

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
                file_path = "./fig/";
            end

            if isempty(file_2)
                file_2 = ".fig";
            end

            zef.file_path = file_path;

            zef.file = file_1 + file_2;

            zef.save_switch = 1;

            zef = zef_import_figure(zef);

        end % for

    end % if

    % Open all figures in a given folder, if given.

    if not(args.open_figure_folder == "")

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

    % Finally, before possibly quitting, run the script given as an argument.
    %
    % NOTE: using eval here is very unsafe. Allows for arbitrary code
    % execution. Make sure given script is from a trusted source.

    if not(args.run_script == "")

        eval(args.run_script);

    end % if

    % Exit zeffiro, if told to.

    if args.exit_zeffiro
        zef_close_all;
    end

    % Close Matlab, if told to.

    if args.quit_matlab
        quit force;
    end

    % Make sure zef exists as a varible in the base workspace.

    if nargout == 0

        assignin("base", "zef", zef);

    end

end % function

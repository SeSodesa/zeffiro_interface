function rec = zef_call_inverse_method( ...
    zef, ...
    inverse_method_handle ...
)

% Calls the scripts needed to initialize a given inverse method and then the
% inversion method itself. Also filters the reconstruction produced by the
% inverse method and returns it.
%
% The returned reconstruction will be empty, if no inverse method was called
% successfully.

    arguments

        zef struct

        inverse_method_handle (1,1) function_handle

    end

    % Initialize empty reconstruction.

    rec = [];

    % Get the name of the given inverse method.

    inverse_method_name = func2str(inverse_method_handle);

    % See if a given inverse method is recognised and if so, call it after
    % running the required initialization script.

    if strcmp(inverse_method_name, 'zef_find_mne_reconstruction')
        [rec, ~] = inverse_method_handle();
    end

    if strcmp(inverse_method_name, 'zef_beamformer')
        % evalin('base', 'zef_init_beamformer;');
        [rec, ~, ~] = inverse_method_handle();
    end

    if strcmp(inverse_method_name, 'zef_ramus_iteration')
        % evalin('base', 'zef_init_ramus_inversion_tool;');
        % evalin('base', 'zef_update_ramus_inversion_tool;');
        [rec, ~] = inverse_method_handle();
    end

    if strcmp(inverse_method_name, 'SESAME_inversion')
        % evalin('base', 'SESAME_App_run');
        % evalin('base', 'zef_SESAME_init;');
        % evalin('base', 'zef_update_SESAME;');
        [rec] = inverse_method_handle();
    end

end

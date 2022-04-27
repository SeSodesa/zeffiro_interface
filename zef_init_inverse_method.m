function zef_init_inverse_method(inverse_method_handle)

    % Get the name of the given inverse method.

    inverse_method_name = func2str(inverse_method_handle);

    % See if a given inverse method is recognised and if so, call it after
    % running the required initialization script.

    if strcmp(inverse_method_name, 'zef_find_mne_reconstruction')
        % No need to do anything for MNE.
    end

    if strcmp(inverse_method_name, 'zef_beamformer')
        evalin('base', 'zef_init_beamformer;');
    end

    if strcmp(inverse_method_name, 'zef_ramus_iteration')
        evalin('base', 'zef_init_ramus_inversion_tool;');
        % evalin('base', '[zef.ramus_multires_dec, zef.ramus_multires_ind, zef.ramus_multires_count] = make_multires_dec();');
        evalin('base', 'zef_update_ramus_inversion_tool;');
    end

    if strcmp(inverse_method_name, 'SESAME_inversion')
        evalin('base', 'SESAME_App_run');
        evalin('base', 'zef_SESAME_init;');
        evalin('base', 'zef_update_SESAME;');
    end
end

function zef = zef_init_inverse_method(zef, inverse_method_handle)

    % Get the name of the given inverse method.

    inverse_method_name = func2str(inverse_method_handle);

    % See if a given inverse method is recognised and if so, call it after
    % running the required initialization script.

    if strcmp(inverse_method_name, 'zef_find_mne_reconstruction')
        % No need to do anything for MNE.
    end

    if strcmp(inverse_method_name, 'zef_beamformer')
        zef_init_beamformer;
    end

    if strcmp(inverse_method_name, 'zef_ramus_iteration')
        zef_init_ramus_inversion_tool;
        zef_update_ramus_inversion_tool;
    end

    if strcmp(inverse_method_name, 'SESAME_inversion')
        SESAME_App_run;
        zef_SESAME_init;
        zef_update_SESAME;
    end
end

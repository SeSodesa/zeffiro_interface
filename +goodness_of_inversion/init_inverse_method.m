function zef = zef_init_inverse_method(zef, inverse_method_name)

    % See if a given inverse method is recognised and if so, call it after
    % running the required initialization script.

    if strcmp(inverse_method_name, "MNE")
        % No need to do anything for MNE.
    end

    if strcmp(inverse_method_name, "Beamformer")
        % zef_init_beamformer;
    end

    if strcmp(inverse_method_name, "RAMUS")
        % zef_init_ramus_inversion_tool;
        % zef_update_ramus_inversion_tool;
    end

    if strcmp(inverse_method_name, "SESAME")
        % SESAME_App_run;
        % zef_SESAME_init;
        % zef_update_SESAME;
    end
end

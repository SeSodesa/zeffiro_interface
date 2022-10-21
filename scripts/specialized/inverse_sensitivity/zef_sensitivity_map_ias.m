function sensitivity_map = zef_sensitivity_map_ias(project_struct, n_reconstructions, noise_level, diff_type)

    arguments

        project_struct (1,1) struct

        n_reconstructions (1,1) double { mustBeInteger, mustBePositive } = 10

        noise_level (1,1) double { mustBeNonpositive } = -30

        diff_type (1,1) string { mustBeMember(diff_type, ["L2", "minabs"]) } = "L2"

    end

    % Set up the beginnings of a return value.

    sensitivity_map = struct;

    sensitivity_map.inverse_method = "IAS";

    project_struct = zef_init_ias(project_struct);

    % Call inverse method for a given number of times for averaging.

    for i = 1 : n_reconstructions

        [sensitivity_map.dist_vec{i},sensitivity_map.angle_vec{i}, sensitivity_map.mag_vec{i}] = zef_rec_diff( ...
            project_struct, ...
            @zef_ias_iteration, ...
            noise_level, ...
            diff_type ...
        );

    end

end % function

function self = preprocess_mesh(self)

    % process_meshes
    %
    % Processes a contained finite element mesh.

    arguments
        self zef_as_class.Zef
    end

    self.mesh_generation_phase = "pre-processing";

    output_mode = 'compact';

    explode_param = self.explode_everything;

    i = 0;
    sensors = [];
    reuna_submesh_ind = cell(0);

    % Adjusting existing compartment positions and directions.

    for k = 1 : length(self.compartments)

        compartment = self.compartments(k);

        on_val = compartment.is_on;

        if on_val

            i = i + 1;

            self.compartments(i).points_inf = compartment.points_inf;
            self.compartments(i).points = compartment.points;
            self.compartments(i).triangles = compartment.triangles;

            reuna_submesh_ind{i} = compartment.submesh_ind;

            mean_vec = repmat(mean(self.compartments(i).points,1),size(self.compartments(i).points,1),1);

            for t_ind = 1 : length(compartment.scaling)

                scaling_val = compartment.scaling(t_ind);

                translation_vec(1) = compartment.x_correction(t_ind);
                translation_vec(2) = compartment.y_correction(t_ind);
                translation_vec(3) = compartment.z_correction(t_ind);

                theta_angle_vec(1) = compartment.xy_rotation(t_ind);
                theta_angle_vec(2) = compartment.yz_rotation(t_ind);
                theta_angle_vec(3) = compartment.zx_rotation(t_ind);

                if ~ isempty(compartment.affine_transform)

                    if length(compartment.affine_transform) >= t_ind
                        affine_transform =  compartment.affine_transform(t_ind);
                    else
                        affine_transform = eye(4);
                    end

                    reuna_aux = [self.compartments(i).points ones(size(self.compartments(i).points,1),1)];

                    if not(isempty(reuna_aux)) && not(isempty(affine_transform))
                        reuna_aux = reuna_aux*affine_transform';
                        self.compartments(i).points = reuna_aux(:,1:3);
                    else
                        self.compartments(i).points = [];
                    end

                    if scaling_val ~= 1
                        self.compartments(i).points = scaling_val*self.compartments(i).points;
                        self.compartments(i).points_inf = scaling_val*self.compartments(i).points_inf;
                    end

                    for j = 1 : 3

                        switch j
                            case 1
                                axes_ind = [1 2];
                            case 2
                                axes_ind = [2 3];
                            case 3
                                axes_ind = [3 1];
                        end % switch

                        if theta_angle_vec(j) ~= 0

                            theta_angle = theta_angle_vec(j)*pi/180;
                            R_mat = [cos(theta_angle) -sin(theta_angle); sin(theta_angle) cos(theta_angle)];
                            self.compartments(i).points(:,axes_ind) = (self.compartments(i).points(:,axes_ind)-mean_vec(:,axes_ind))*R_mat' + mean_vec(:,axes_ind);
                            self.compartments(i).points_inf(:,axes_ind) = (self.compartments(i).points_inf(:,axes_ind)-mean_vec(:,axes_ind))*R_mat' + mean_vec(:,axes_ind);

                        end % if

                    end % for

                    for j = 1 : 3

                        if translation_vec(j) ~= 0
                            self.compartments(i).points(:,j) = self.compartments(i).points(:,j) + translation_vec(j);
                            self.compartments(i).points_inf(:,j) = self.compartments(i).points_inf(:,j) + translation_vec(j);
                        end

                    end % for

                end % if

                if explode_param ~= 1

                    for s_ind = 1 : length(reuna_submesh_ind{i})

                        if s_ind == 1
                            t_ind_1 = 1;
                        else
                            t_ind_1 = reuna_submesh_ind{i}(s_ind-1)+1;
                        end

                        t_ind_2 = reuna_submesh_ind{i}(s_ind);
                        p_ind = unique(self.compartments(i).triangles(t_ind_1:t_ind_2,:));
                        mean_aux = mean(self.compartments(i).points(p_ind,:),1);
                        self.compartments(i).points(p_ind,:) = self.compartments(i).points(p_ind,:) + (explode_param-1)*repmat(mean_aux,length(p_ind),1);

                        if not(isempty(self.compartments(i).points_inf))
                            mean_aux = mean(self.compartments(i).points_inf(p_ind,:),1);
                            self.compartments(i).points_inf(p_ind,:) = self.compartments(i).points_inf(p_ind,:) + (explode_param-1)*repmat(mean_aux,length(p_ind),1);
                        end

                    end % for

                end % if

            end % for

        end % if

        % Adjusting sensor positions and orientations.

        sensor_tag = self.data.current_sensors;

        sensor_tag_pts_name = [sensor_tag '_points'];

        sensor_tag_dirs_name = [sensor_tag '_directions'];

        s_points = self.data.(sensor_tag_pts_name);

        s_data_aux = [];

        if ismember(self.data.imaging_method, 1)

            f_handle = self.create_sensor_fn;

            if not(isempty(f_handle))
                % s_points = f_handle(s_points);
            end

        end % if

        if ismember(self.data.imaging_method, [2 3])

            s_directions = self.data.(sensor_tag_dirs_name)(:,1:3);

            s_directions_g = [];

            if size(self.data.(sensor_tag_dirs_name), 2) == 6
                s_directions_g = self.data.(sensor_tag_dirs_name)(:,4:6);
            end

        else

            if size(s_points,2)==6
                s_data_aux = s_points(:,4:6);
                s_points = s_points(:,1:3);
            end

            s_directions = [];
            s_directions_g = [];

        end % if

        sensor_tag_scaling_name = [sensor_tag '_scaling'];
        sensor_tag_x_corr_name = [sensor_tag '_x_correction'];
        sensor_tag_y_corr_name = [sensor_tag '_y_correction'];
        sensor_tag_z_corr_name = [sensor_tag '_z_correction'];
        sensor_tag_xy_rot_name = [sensor_tag '_xy_rotation'];
        sensor_tag_yz_rot_name = [sensor_tag '_yz_rotation'];
        sensor_tag_zx_rot_name = [sensor_tag '_zx_rotation'];
        sensor_tag_aff_trans_name = [sensor_tag '_affine_transform'];

        s_scaling = self.data.(sensor_tag_scaling_name);
        s_x_correction = self.data.(sensor_tag_x_corr_name);
        s_y_correction = self.data.(sensor_tag_y_corr_name);
        s_z_correction = self.data.(sensor_tag_z_corr_name);
        s_xy_rotation = self.data.(sensor_tag_xy_rot_name);
        s_yz_rotation = self.data.(sensor_tag_yz_rot_name);
        s_zx_rotation = self.data.(sensor_tag_zx_rot_name);

        if isfield(self.data, sensor_tag_aff_trans_name)
            s_affine_transform = self.data.(sensor_tag_aff_trans_name);
        else
            s_affine_transform = cell(0);
        end

        if isempty(s_directions)

            sensors = [s_points];

        else

            if isempty(s_directions_g)
                sensors = [s_points s_directions./repmat(sqrt(sum(s_directions.^2,2)),1,3)];
            else
                sensors = [s_points s_directions./repmat(sqrt(sum(s_directions.^2,2)),1,3) s_directions_g./repmat(sqrt(sum(s_directions_g.^2,2)),1,3)];
            end

        end % if

        use_pem = 0;

        if ismember(self.data.imaging_method, [1 5])
            use_pem = self.data.use_pem;
        end

        for t_ind = 1 : length(s_scaling)

            scaling_val = s_scaling(t_ind);
            translation_vec = [s_x_correction(t_ind) s_y_correction(t_ind) s_z_correction(t_ind)];
            theta_angle_vec = [s_xy_rotation(t_ind) s_yz_rotation(t_ind) s_zx_rotation(t_ind)];

            if not(isempty(sensors))

                mean_vec = repmat(mean(sensors(:,1:3),1),size(sensors(:,1:3),1),1);

                if length(s_affine_transform) >= t_ind

                    sensors_aux = [sensors(:,1:3) ones(size(sensors,1),1)];
                    sensors_aux = sensors_aux*s_affine_transform{t_ind}';
                    sensors(:,1: 3) = sensors_aux(:,1:3);

                    if size(sensors_aux,2) >= 6 && ismember(self.data.imaging_method, [2 3])
                        sensors(:,4:6) = sensors(:,4:6)*s_affine_transform{t_ind}(1:3,1:3)';
                    end

                    if size(sensors_aux,2) >= 9 && ismember(self.data.imaging_method, [3])
                        sensors(:,7:9) = sensors(:,7:9)*s_affine_transform{t_ind}(1:3,1:3)';
                    end

                end % if

                if scaling_val ~= 1
                    sensors(:,1:3) = scaling_val * sensors(:,1:3);
                end

                for j = 1 : 3

                    switch j
                        case 1
                            axes_ind_1 = [1 2];
                            axes_ind_2 = [4 5];
                            axes_ind_3 = [7 8];
                        case 2
                            axes_ind_1 = [2 3];
                            axes_ind_2 = [5 6];
                            axes_ind_3 = [8 9];
                        case 3
                            axes_ind_1 = [3 1];
                            axes_ind_2 = [6 4];
                            axes_ind_3 = [9 7];
                    end % switch

                    if theta_angle_vec(j) ~= 0

                        theta_angle = theta_angle_vec(j)*pi/180;

                        R_mat = [cos(theta_angle) -sin(theta_angle); sin(theta_angle) cos(theta_angle)];

                        sensors(:,axes_ind_1) = (sensors(:,axes_ind_1) - mean_vec(:,axes_ind_1))*R_mat' + mean_vec(:,axes_ind_1);

                        if not(isempty(s_directions))
                            sensors(:,axes_ind_2) = sensors(:,axes_ind_2)*R_mat';
                        end

                        if not(isempty(s_directions_g))
                            sensors(:,axes_ind_3) = sensors(:,axes_ind_3)*R_mat';
                        end

                    end % if

                end % for

                for j = 1 : 3

                    if translation_vec(j) ~= 0

                        sensors(:,j) = sensors(:,j) + translation_vec(j);

                    end

                end % for

            end % if

        end % for

        if not(isempty(sensors))

            for j = 1 : 3
                sensors(:,j) = sensors(:,j)  +  (explode_param-1)*(sensors(:,j) - mean(sensors(:,j),1));
            end

            if not(isempty(s_data_aux))
                sensors = [sensors s_data_aux];
            end

        else

            sensors = [NaN NaN NaN];

        end % if

        if use_pem
            sensors = sensors(:,1:3);
        end

        % More compartment adjustments.

        max_val = 0;
        box_ind = 0;

        for i_aux = 1 : length(self.compartments)

            compartment = self.compartments(i_aux);

            if not(isequal(compartment.sources, -1))

                max_val = max([abs(compartment.points(:)) ; max_val]);

            elseif isequal(compartment.sources, -1)

                box_ind = i_aux;

            end

        end % for

        if not(isequal(box_ind,0))

            if self.pml_outer_radius_unit == 1
                box_outer_radius = self.pml_outer_radius * max_val;
            elseif self.pml_outer_radius_unit == 2
                box_outer_radius = self.pml_outer_radius;
            else
                error("Unknown PML outer radius unit.");
            end

            self.compartments(box_ind).points = [
                -box_outer_radius   -box_outer_radius   -box_outer_radius ;
                 box_outer_radius   -box_outer_radius   -box_outer_radius ;
                 box_outer_radius    box_outer_radius   -box_outer_radius ;
                -box_outer_radius    box_outer_radius   -box_outer_radius ;
                -box_outer_radius   -box_outer_radius    box_outer_radius ;
                 box_outer_radius   -box_outer_radius    box_outer_radius ;
                 box_outer_radius    box_outer_radius    box_outer_radius ;
                -box_outer_radius    box_outer_radius    box_outer_radius ;
            ];

            self.compartments(box_ind).triangles = [
                1  2  6;
                6  5  1;
                3  4  7;
                8  7  4;
                2  3  7;
                2  7  6;
                5  1  8;
                4  1  8;
                2  1  4;
                4  3  2;
                5  6  8;
                7  8  6;
            ];

        end % if

    end % for

    % Place modified sensors into self.

    if size(sensors, 2) == 3

        self.sensors = zef_as_class.Sensor.from_arrays( ...
            sensors(:,1:3), ...
            s_directions, [], [], [] ...
        );

    elseif size(sensors, 2) == 6

        self.sensors = zef_as_class.Sensor.from_arrays( ...
            sensors(:,1:3), ...
            s_directions, ...
            sensors(:,4), ...
            sensors(:,5), ...
            sensors(:,6) ...
        );

    else

        self.sensors = [];

    end

end % function

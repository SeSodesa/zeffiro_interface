function self = create_fem_mesh(self)

    % create_fem_mesh
    %
    % Generates a finite element mesh from a surface segmentation.

    arguments
        self zef_as_class.Zef
    end

    void = [];

    self = self.segmentation_counter_step();

    % Determine mesh limits in a Cartesian system.

    x_lim = [0 0];
    y_lim = [0 0];
    z_lim = [0 0];

    n_compartments = 0;

    for k = 1 : length(self.compartments)

        compartment = self.compartments(k);

        n_compartments = n_compartments + max(1,length(compartment.submesh_ind));

        reuna_p = compartment.points;

        if not(isequal(compartment.sources,-1))
            x_lim = [min(x_lim(1),min(reuna_p(:,1))) max(x_lim(2),max(reuna_p(:,1)))];
            y_lim = [min(y_lim(1),min(reuna_p(:,2))) max(y_lim(2),max(reuna_p(:,2)))];
            z_lim = [min(z_lim(1),min(reuna_p(:,3))) max(z_lim(2),max(reuna_p(:,3)))];
        end

    end % for

    % Generate a grid based on the limits and the mesh resolution, either as
    % cubes or perfectly matched layers (PML).

    mesh_res = self.mesh_resolution;

    if isempty(self.data.pml_ind_aux)

        x_vec = [x_lim(1):mesh_res:x_lim(2)];
        y_vec = [y_lim(1):mesh_res:y_lim(2)];
        z_vec = [z_lim(1):mesh_res:z_lim(2)];

        [X, Y, Z] = meshgrid(x_vec,y_vec,z_vec);

        n_cubes = (length(x_vec)-1)*(length(y_vec)-1)*(length(z_vec)-1);

    else

        pml_inner_radius = max(abs([x_lim(:); y_lim(:); z_lim(:)]));
        pml_outer_radius_unit = self.pml_outer_radius_unit;
        pml_outer_radius = self.pml_outer_radius;
        pml_max_size_unit = self.pml_max_size_unit;
        pml_max_size = self.pml_max_size;

        if isequal(pml_outer_radius_unit,1)
            pml_outer_radius = pml_inner_radius*pml_outer_radius;
        end

        if isequal(pml_max_size_unit,1)
            pml_max_size = mesh_res*pml_max_size;
        end

        [X, Y, Z, pml_ind] = pml_mesh(pml_inner_radius,pml_outer_radius,mesh_res,pml_max_size);

        n_cubes = prod(size(X)-1);

    end % if

    size_xyz = size(X);

    wb = waitbar(0,'Initial mesh.');

    cleanup_fn = @(h) close(h);

    cleanup_obj = onCleanup(@() cleanup_fn(wb));

    % Construct initial mesh.

     [self.nodes, self.tetra, self.label_ind] = initial_nodes_tetra_and_label_ind(self, X, Y, Z, n_cubes);

    % Start labeling generated tetra.

    self.labeling_mode = "initial";

    self = mesh_labeling_step(self);

    % Then perform mesh refinement, either surface and/or volumetric, in
    % adapted and/or non-adapted mode.

    refinement_compartments_aux = self.data.refinement_surface_compartments;

    refinement_compartments = [];

    if ismember(1,refinement_compartments_aux)
        refinement_compartments = aux_brain_ind(:);
    end

    refinement_compartments_aux = setdiff(refinement_compartments_aux,1)-1;

    refinement_compartments = [refinement_compartments ; refinement_compartments_aux(:)];

    refinement_flag = 1;

    n_surface_refinement = self.data.refinement_surface_number;

    if self.refinement_on

        if self.surface_refinement_on

            if length(n_surface_refinement) == 1

                for i_surface_refinement = 1 : n_surface_refinement

                    self = self.mesh_refinement_step();

                    if self.data.mesh_relabeling

                        pml_ind = [];
                        label_ind = uint32(tetra);
                        self.labeling_mode = "repeated";
                        self = mesh_labeling_step(self);

                    end % if

                end % for

            else

                for j_surface_refinement = 1 : length(n_surface_refinement)

                    for i_surface_refinement = 1 : n_surface_refinement(j_surface_refinement)

                        self = self.mesh_refinement_step();

                        if self.data.mesh_relabeling

                            pml_ind = [];
                            label_ind = uint32(tetra);
                            self.labeling_mode = "repeated";
                            self = mesh_labeling_step(self);

                        end % if

                    end % for

                end % for

            end % if

        end % if

    else

        pml_ind = [];
        label_ind = uint32(tetra);
        self.labeling_mode = "repeated";
        self = mesh_labeling_step(self);

    end % if

    if self.volume_refinement_on

        n_refinement = self.data.refinement_volume_number;
        refinement_compartments_aux = sort(self.data.refinement_volume_compartments);

        refinement_compartments = [];

        if ismember(1,refinement_compartments_aux)
            refinement_compartments = aux_brain_ind(:);
        end

        refinement_compartments_aux = setdiff(refinement_compartments_aux,1)-1;
        refinement_compartments = [refinement_compartments ; refinement_compartments_aux(:)];

        if length(n_refinement) == 1

            waitbar(0,wb,'Volume refinement.');

            for i = 1 : n_refinement

                [nodes,tetra,self.domain_labels] = zef_mesh_refinement(nodes,tetra,self.domain_labels,refinement_compartments);

                waitbar(i/n_refinement,wb,'Volume refinement.');

                if self.data.mesh_relabeling

                    pml_ind = [];
                    label_ind = uint32(tetra);
                    self.labeling_mode = "repeated";
                    self = mesh_labeling_step(self);

                end % if

            end % for

        else

            waitbar(0/length(n_refinement),wb,'Volume refinement.');

            for j = 1 : length(n_refinement)

                for i = 1 : n_refinement(j)

                    [nodes,tetra,self.domain_labels] = zef_mesh_refinement(nodes,tetra,self.domain_labels,refinement_compartments(j));

                    if self.data.mesh_relabeling

                        pml_ind = [];
                        label_ind = uint32(tetra);
                        self.labeling_mode = 2;
                        self = mesh_labeling_step(self);

                    end % if

                    waitbar(i/length(n_refinement(j)),wb,'Volume refinement.');

                end % for

            end % if

        end % if

    end % if

    if self.adaptive_refinement_on

        n_refinement = self.data.adaptive_refinement_number;

        refinement_compartments_aux = sort(self.data.adaptive_refinement_compartments);

        refinement_compartments = [];

        if ismember(1,refinement_compartments_aux)
            refinement_compartments = aux_brain_ind(:);
        end

        refinement_compartments_aux = setdiff(refinement_compartments_aux,1)-1;
        refinement_compartments = [refinement_compartments ; refinement_compartments_aux(:)];

        if length(n_refinement) == 1

            waitbar(0,wb,'Adaptive volume refinement.');

            for i = 1 : n_refinement

                k_param = self.data.adaptive_refinement_k_param;

                thresh_val  = self.data.adaptive_refinement_thresh_val;

                tetra_refine_ind = zef_get_tetra_to_refine(refinement_compartments, thresh_val, k_param, nodes, tetra,self.domain_labels,reuna_p,reuna_t);

                [nodes,tetra,self.domain_labels,tetra_interp_vec] = zef_mesh_refinement(nodes,tetra,self.domain_labels,refinement_compartments, tetra_refine_ind);

                tetra_refine_ind = find(ismember(tetra_interp_vec,tetra_refine_ind));

                waitbar(i/n_refinement,wb,'Adaptive volume refinement.');

                if self.data.mesh_relabeling

                    pml_ind = [];
                    label_ind = uint32(tetra);
                    self.labeling_mode = "adaptive-repeated";
                    self = mesh_labeling_step(self);

                end % if

            end % for

        else

            waitbar(0/length(n_refinement),wb,'Adaptive volume refinement.');

            for j = 1 : length(n_refinement)

                for i = 1 : n_refinement(j)

                    k_param = self.adaptive_refinement_k_param;

                    thresh_val  = self.adaptive_refinement_thresh_val;

                    tetra_refine_ind = zef_get_tetra_to_refine(refinement_compartments(j), thresh_val, k_param, nodes, tetra,self.domain_labels,reuna_p,reuna_t);

                    [nodes,tetra,self.domain_labels,tetra_interp_vec] = mesh_refinement(nodes,tetra,self.domain_labels,refinement_compartments(j),tetra_refine_ind);

                    tetra_refine_ind = find(ismember(tetra_interp_vec,tetra_refine_ind));

                    if self.data.mesh_relabeling

                        pml_ind = [];
                        label_ind = uint32(tetra);
                        self.labeling_mode = "adaptive-repeated";
                        self = mesh_labeling_step(self);

                    end % if

                    waitbar(i/length(n_refinement(j)),wb,'Adaptive volume refinement.');

                end % for

            end % for

        end % if

    end % if

% if nargout == 0
%     assignin('base','zef_data',struct('nodes',nodes,'nodes_raw',nodes,'tetra',tetra,'tetra_raw',tetra,'self.domain_labels',self.domain_labels,'domain_labels_raw',self.domain_labels,'name_tags',name_tags));
% end

end % function

%% Local helper functions

function [nodes, tetra, label_ind] = initial_nodes_tetra_and_label_ind( ...
    zef, ...
    X, ...
    Y, ...
    Z, ...
    n_cubes ...
)

    % initial_nodes_tetra_and_label_ind
    %
    % Constructs an initial mesh in mode 1, whatever that means.
    %
    % Input:
    %
    % - zef
    %
    % - X
    %
    % - Y
    %
    % - Z
    %
    % - n_cubes
    %
    % Output:
    %
    % - nodes
    %
    % - tetra
    %
    % - label_ind
    %

    arguments
        zef zef_as_class.Zef
        X (:,:,:) double
        Y (:,:,:) double
        Z (:,:,:) double
        n_cubes (1,1) double { mustBeInteger, mustBePositive }
    end

    % Initialize waitbar and cleanup object.

    if zef.use_gui

        wb = waitbar(0, "Initial mesh.");

        cu_fn = @(h) close(h);

        cu_obj = onCleanup(@() cu_fn(wb));

    else

        wb = zef_as_class.TerminalWaitbar("Initial mesh", size(X,2) - 1);

    end

    % Routine start.

    if zef.initial_mesh_mode == 1

        ind_mat_1{1}{2}{1} = [ 2 5 6 7 ; 7 5 4 2 ; 2 3 4 7 ; 1 2 4 5 ; 4 7 8 5 ];
        ind_mat_1{1}{2}{2} = [ 6 2 1 3 ; 1 3 8 6 ; 8 7 6 3 ; 5 8 6 1 ; 3 8 4 1 ];
        ind_mat_1{2}{2}{2} = [ 5 2 1 4 ; 4 2 7 5 ; 5 8 7 4 ; 5 7 6 2 ; 3 7 4 2 ];
        ind_mat_1{2}{2}{1} = [ 1 5 6 8 ; 6 8 3 1 ; 3 4 1 8 ; 2 3 1 6 ; 3 7 8 6 ];
        ind_mat_1{1}{1}{2} = [ 4 3 7 2 ; 2 7 4 5 ; 5 7 6 2 ; 1 5 2 4 ; 8 7 5 4 ];
        ind_mat_1{2}{1}{2} = [ 3 6 8 1 ; 1 3 4 8 ; 5 8 6 1 ; 1 6 2 3 ; 8 7 6 3 ];
        ind_mat_1{1}{1}{1} = [ 7 8 3 6 ; 8 1 3 6 ; 2 3 1 6 ; 1 5 6 8 ; 1 3 4 8 ];
        ind_mat_1{2}{1}{1} = [ 7 8 4 5 ; 5 4 7 2 ; 2 4 1 5 ; 2 5 6 7 ; 2 3 4 7 ];

        tetra = zeros(5*n_cubes,4);

        if isequal(zef.mesh_labeling_approach,1)

            label_ind = zeros(5*n_cubes,8);

        elseif isequal(zef.mesh_labeling_approach, 2)

            label_ind = zeros(5*n_cubes,4);

        else

            error("Unknown mesh labeling approach.");

        end

    elseif zef.initial_mesh_mode == 2

        ind_mat_1 = [
            3     4     1     7 ;
            2     3     1     7 ;
            1     2     7     6 ;
            7     1     6     5 ;
            7     4     1     8 ;
            7     8     1     5
        ];

        tetra = zeros(6*n_cubes,4);

        if isequal(zef.mesh_labeling_approach,1)

            label_ind = zeros(6*n_cubes,8);

        elseif isequal(zef.mesh_labeling_approach, 2)

            label_ind = zeros(6*n_cubes,4);

        else

            error("Unknown mesh labeling approach.");

        end

    else

        error("Unknown initial mesh mode.")

    end

    nodes = [X(:) Y(:) Z(:)];

    i = 1;

    size_xyz = size(X);

    for i_x = 1 : size(X,2) - 1

        if zef.use_gui

            waitbar(i_x/(size(X,2)-1),wb,'Initial mesh.');

        else

            wb = wb.progress();

        end

        for i_y = 1 : size(X,1) - 1

            for i_z = 1 : size(X,3) - 1

                x_ind = [i_x   i_x+1  i_x+1  i_x    i_x    i_x+1  i_x+1  i_x]';
                y_ind = [i_y   i_y    i_y+1  i_y+1  i_y    i_y    i_y+1  i_y+1]';
                z_ind = [i_z   i_z    i_z    i_z    i_z+1  i_z+1  i_z+1  i_z+1]';

                ind_mat_2 = sub2ind(size_xyz,y_ind,x_ind,z_ind);

                if zef.initial_mesh_mode == 1

                    tetra(i:i+4,:) = ind_mat_2(ind_mat_1{2-mod(i_x,2)}{2-mod(i_y,2)}{2-mod(i_z,2)});

                    if isequal(zef.mesh_labeling_approach, 1)

                        label_ind(i:i+4,:) = ind_mat_2(:,ones(5,1))';

                    elseif isequal(zef.mesh_labeling_approach, 2)

                        label_ind(i:i+4,:) = ind_mat_2(ind_mat_1{2-mod(i_x,2)}{2-mod(i_y,2)}{2-mod(i_z,2)});

                    else

                        error("Unknown mesh labeling approach.");

                    end

                    i = i + 5;

                elseif zef.initial_mesh_mode == 2

                    tetra(i:i+5,:) = ind_mat_2(ind_mat_1);

                    if isequal(zef.mesh_labeling_approach, 1)

                        label_ind(i:i+5,:) = ind_mat_2(:,ones(6,1))';

                    elseif isequal(zef.mesh_labeling_approach, 2)

                        label_ind(i:i+5,:) = ind_mat_2(ind_mat_1);

                    else

                        error("Unknown mesh labeling approach.");

                    end;

                    i = i + 6;

                else

                    error("Unknown initial mesh mode.");

                end

            end % for

        end % for

    end % for

end % function

function [X, Y, Z, pml_ind] = pml_mesh(inner_radius,outer_radius,lattice_size,max_size)

    % pml_mesh
    %
    % TODO
    %
    % Inputs:
    %
    % - inner_radius
    %
    % - outer_radius
    %
    % - lattice_size
    %
    % - max_size
    %
    % Outputs:
    %
    % - X
    %
    % - Y
    %
    % - Z
    %
    % - pml_ind

    arguments
        inner_radius (1,1) double { mustBePositive }
        outer_radius (1,1) double { mustBePositive }
        lattice_size (1,1) double { mustBePositive }
        max_size     (1,1) double { mustBePositive }
    end

    growth_param = 1 + 1e-15;

    convergence_value = Inf;

    while convergence_value > 1e-5

        extra_layers = round(log(((outer_radius-inner_radius)*(growth_param - 1))/(lattice_size*growth_param) + 1)/log(growth_param));

        growth_param_new = exp(log(max_size/lattice_size)/extra_layers);

        convergence_value = abs(max_size - lattice_size*growth_param^extra_layers)/max_size;

        growth_param = growth_param_new;

    end

    extra_layers = round(log(((outer_radius-inner_radius)*(growth_param - 1))/(lattice_size*growth_param) + 1)/log(growth_param));

    intermediate_radius = extra_layers*lattice_size + inner_radius;

    [X,Y,Z] = meshgrid([-intermediate_radius:2*intermediate_radius/(round(2*intermediate_radius/lattice_size)):intermediate_radius]);

    I_X = find(abs(X) > inner_radius);
    I_Y = find(abs(Y) > inner_radius);
    I_Z = find(abs(Z) > inner_radius);

    pml_ind = unique([I_X(:); I_Y(:) ; I_Z(:)]);

    R_aux = max(abs([X(pml_ind) Y(pml_ind) Z(pml_ind)]),[],2);

    N = round((R_aux - inner_radius)/lattice_size);

    X_inner = inner_radius.*X(pml_ind)./R_aux;
    Y_inner = inner_radius.*Y(pml_ind)./R_aux;
    Z_inner = inner_radius.*Z(pml_ind)./R_aux;

    R_new = inner_radius + lattice_size.*growth_param*(1-growth_param.^N)./(1-growth_param);

    R_old = inner_radius + N*lattice_size;

    X(pml_ind) = R_new.*X(pml_ind)./R_old;
    Y(pml_ind) = R_new.*Y(pml_ind)./R_old;
    Z(pml_ind) = R_new.*Z(pml_ind)./R_old;

    s = max(abs([X(:) ; Y(:); Z(:)]));
    X = outer_radius*X/s;
    Y = outer_radius*Y/s;
    Z = outer_radius*Z/s;

end

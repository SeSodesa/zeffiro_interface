function self = segmentation_counter_step(self)

    arguments
        self zef_as_class.Zef
    end

    n_of_iterations = length(self.compartments);

    fn_title = "Segmentation counter";

    if self.use_gui

        wb = waitbar(0, fn_title);

        cleanup_fn = @(h) close(h);

        cleanup_obj = onCleanup(@() cleanup_fn(wb));

    else

        wb = zef_as_class.TerminalWaitbar(fn_title, n_of_iterations);

    end

    pml_ind_aux = [];
    pml_ind = [];

    mesh_res = self.mesh_resolution;
    sensors = self.sensors;

    i = 0;

    sigma_vec = [];
    priority_vec = [];
    submesh_cell = cell(0);
    aux_brain_ind = [];

    for k = 1 : length(self.compartments)

        if self.use_gui

            waitbar(i_x/(size(X,2)-1),wb, fn_title);

        else

            wb = wb.progress();

        end

        compartment = self.compartments(k);

        on_val = compartment.is_on;
        sigma_val = compartment.sigma;
        priority_val = compartment.priority;

        if on_val

            i = i + 1;

            sigma_vec(i,1) = sigma_val;
            priority_vec(i,1) = priority_val;
            self.compartments(i).submesh_ind = compartment.submesh_ind_original_surface_mesh;
            self.compartments(i).name = compartment.name;

            if isequal(compartment.sources,-1)
                pml_ind_aux = i;
            end

            if ismember(compartment.sources,[1 2])
                aux_brain_ind = [aux_brain_ind i];
            end

        end % if

    end % for

    n_compartments = 0;

    for k = 1 : length(self.compartments)

        compartment = self.compartments(k);

        n_compartments = n_compartments + max(1,length(compartment.submesh_ind));

    end

    priority_vec_aux = zeros(n_compartments,1);
    compartment_counter = 0;
    submesh_ind_1 = ones(n_compartments,1);
    submesh_ind_2 = ones(n_compartments,1);

    % Re-organize compartment priorities.

    old_priorities = arrayfun(@(c) c.priority, self.compartments);

    for i = 1 : length(self.compartments)

        compartment = self.compartments(i);

        for k = 1 : max(1,length(compartment.submesh_ind))

            compartment_counter = compartment_counter + 1;
            self.compartments(compartment_counter).priority = old_priorities(i);
            submesh_ind_1(compartment_counter) = i;
            submesh_ind_2(compartment_counter) = k;

        end % for

    end % for

    self.data.submesh_ind_1 = submesh_ind_1;

    self.data.submesh_ind_2 = submesh_ind_2;

    self.data.submesh_cell = submesh_cell;

    self.data.pml_ind_aux = pml_ind_aux;

    self.data.pml_ind = pml_ind;

end % function

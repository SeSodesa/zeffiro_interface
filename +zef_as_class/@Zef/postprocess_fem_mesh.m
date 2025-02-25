function self = postprocess_fem_mesh(self)

    % postprocess_fem_mesh
    %
    % Post-processes a given finite element mesh.

    arguments
        self zef_as_class.Zef
    end

    self.mesh_generation_phase = "post-processing";

end

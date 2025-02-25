classdef VolumeCompartment < zef_as_class.AffinelyTransformable

    % zef_as_class.VolumeCompartment
    %
    % A class that describes how volumetric compartments are structured and
    % what operations are defined on them.
    %
    % Properties:
    %
    % - affine_transform
    %
    %   A transformation matrix used during the mesh construction of this
    %   compartment.
    %
    % - color
    %
    %   The RGB color of this compartment. Should probably not be stored here, but
    %   in the GUI.
    %
    % - condition_number
    %
    %   The condition number of this compartment.
    %
    % - delta
    %
    %   A parameter that might be used for something.
    %
    % - epsilon
    %
    %   A parameter that might be used for something.
    %
    % - is_on
    %
    %   A boolean flag for denoting whether this volumetric compartment is
    %   active.
    %
    % - name
    %
    %   The textual label of this compartment.
    %
    % - scaling
    %
    %   A scalar which denotes whether this compartment is scaled.
    %
    % - points
    %
    %   The set of points that belong to this volume compartment.
    %
    % - points_inf
    %
    %   Some set of points.
    %
    % - points_original_surface_mesh
    %
    %   The set of points that belong to this volume compartment, before being
    %   modified by the mesh creation process.
    %
    % - priority
    %
    %   Describes the priority of this compartment.
    %
    % - sigma
    %
    %   The condutivity of this compartment.
    %
    % - sources
    %
    %   A number indicating the number of sources in this volume.
    %
    % - sources_old
    %
    %   A number indicating the number of sources in this volume, before being
    %   modified by the mesh creation routine.
    %
    % - submesh_ind
    %
    %   The indices of the submeshes that this compartment is a member of?
    %
    % - submesh_ind_original_surface_mesh
    %
    %   The indices of the submeshes that this compartment is a member of,
    %   before being modified by the mesh generation routine?
    %
    % - transform_name
    %
    %   The name of the affine transform applied to this compartment.
    %
    % - triangles
    %
    %   The index set signifying which triangles or triples of nodes belong to
    %   this compartment.
    %
    % - triangles_original_surface_mesh
    %
    %   The index set signifying which triangles or triples of nodes belong to
    %   this compartment, before being modified by the mesh creation process.
    %
    % - x_correction
    %
    %   A scalar which indicates whether the compartment has been corrected in
    %   the x-direction in a Cartesian coordinate system.
    %
    % - xy_rotation
    %
    %   A scalar which indicates how much this compartment has been rotated on
    %   relation to the xy-plane in a Cartesian system.
    %
    % - y_correction
    %
    %   A scalar which indicates whether the compartment has been corrected in
    %   the y-direction in a Cartesian coordinate system.
    %
    % - yz_rotation
    %
    %   A scalar which indicates how much this compartment has been rotated on
    %   relation to the yz-plane in a Cartesian system.
    %
    % - z_correction
    %
    %   A scalar which indicates whether the compartment has been corrected in
    %   the z-direction in a Cartesian coordinate system.
    %
    % - zx_rotation
    %
    %   A scalar which indicates how much this compartment has been rotated on
    %   relation to the zx-plane in a Cartesian system.
    %

    properties

        color (1,3) double {} = [0 0 0];

        condition_number (1,1) double = 1;

        delta (1,1) double = 0;

        epsilon (1,1) double = 0;

        is_on (1,1) logical = true;

        name (1,1) string = "";

        points_inf (:,3) double = [];

        points (:,3) double = [];

        points_original_surface_mesh (:,3) double = [];

        priority (1,1) double { mustBeInteger, mustBePositive } = 1;

        triangles (:,3) double { mustBeInteger, mustBePositive } = [];

        triangles_original_surface_mesh (:,3) double { mustBeInteger, mustBePositive } = [];

        sigma (1,1) double { mustBeNonnegative } = 1;

        submesh_ind double { mustBeInteger, mustBePositive } = [];

        submesh_ind_original_surface_mesh double { mustBeInteger, mustBePositive } = [];

        sources (1,1) double { mustBeInteger } = 0;

        sources_old (1,1) double { mustBeInteger }  = 0;

    end % properties

    properties (Constant)

        % The names of the volumentric compartments that can be found in some
        % older Zeffiro project files. The point of this is to help with
        % backwards-compatibility.

        VOLUME_COMPARTMENT_FIELD_NAMES = [
            "affine_transform",
            "color",
            "condition_number",
            "delta",
            "epsilon",
            "invert",
            "merge",
            "mu",
            "name",
            "on",
            "points",
            "points_original_surface_mesh",
            "points_inf",
            "priority",
            "rho",
            "scaling",
            "sigma",
            "sources",
            "sources_old",
            "submesh_ind",
            "submesh_ind_original_surface_mesh",
            "tetra",
            "transform_name",
            "triangles",
            "triangles_original_surface_mesh",
            "visible",
            "x_correction",
            "xy_rotation",
            "y_correction",
            "yz_rotation",
            "z_correction",
            "zx_rotation"
        ];

    end % properties (Constant)

    methods

        function self = VolumeCompartment(args)

            % zef_as_class.VolumeCompartment.VolumeCompartment
            %
            % Constructor for this class. Takes as an input a struct whose
            % contents will be converted to the properties of this
            % compartment.

            arguments

                args struct

            end

            % Initialize superclass fields.

            self = self@zef_as_class.AffinelyTransformable(args);

            % Initialize own fields.

            if isfield(args, "color")

                try

                    self.color = args.color;

                catch

                    warning("Could not initialize self.is_on from legacy file. Using default value.")

                end

            end

            if isfield(args, "condition_number")

                self.condition_number = args.condition_number;

            end

            if isfield(args, "delta")

                self.delta = args.delta;

            end

            if isfield(args, "epsilon")

                self.epsilon = args.epsilon;

            end

            if isfield(args, "on")

                try

                    self.is_on = args.on;

                catch

                    warning("Could not initialize self.is_on from legacy file. Using default value.")

                end

            end

            if isfield(args, "name")

                self.name = args.name;

            end

            if isfield(args, "points")

                self.points = args.points;

            end

            if isfield(args, "points_inf")

                self.points_inf = args.points_inf;

            end

            if isfield(args, "points_original_surface_mesh")

                self.points_original_surface_mesh = args.points_original_surface_mesh;

            end

            if isfield(args, "priority")

                self.priority = args.priority;

            end

            if isfield(args, "sigma")

                self.sigma = args.sigma;

            end

            if isfield(args, "sources")

                self.sources = args.sources;

            end

            if isfield(args, "sources_old")

                self.sources_old = args.sources_old;

            end

            if isfield(args, "submesh_ind")

                self.submesh_ind = args.submesh_ind;

            end

            if isfield(args, "submesh_ind_original_surface_mesh")

                self.submesh_ind_original_surface_mesh = args.submesh_ind_original_surface_mesh;

            end

            if isfield(args, "triangles")

                self.triangles = args.triangles;

            end

            if isfield(args, "triangles_original_surface_mesh")

                self.triangles_original_surface_mesh = args.triangles_original_surface_mesh;

            end

        end % function

    end % methods

    methods (Static)

        function compartments = from_legacy_zef(zef)

            %
            % from_legacy_zef
            %
            % Builds an instace of VolumeCompartment from a legacy version of
            % the zef struct.
            %
            % Input:
            %
            % - zef
            %
            %   A struct that contains the necessary compartment info:
            %
            %   - compartment_tags
            %
            %   - The fields names listed in VolumeCompartment.VOLUME_COMPARTMENT_FIELD_NAMES,
            %     prefixed with the tags in compartment_tags.
            %
            % Output
            %
            % - compartments
            %
            %   An array of VolumeCompartment objects.
            %

            arguments

                zef (1,1) struct

            end

            % Check for missing fields and other errors.

            if ~ isfield(zef, "compartment_tags")

                error("A legacy zef struct needs to contain a field compartment_tags. Aborting construction of VolumeCompartment...")

            end

            if isempty(zef.compartment_tags)

                error("Cannot construct VolumeCompartments from an empty zef.compartment_tags field.")

            end

            % Loop over the fields to construct the compartment tables, after
            % compartment tags are known.

            field_names = fieldnames(zef);

            compartment_table = cell(0,2);

            for fi = 1 : length(field_names)

                field_name = field_names{fi};

                field_value = zef.(field_name);

                % Parse the field name to see if it is a known compartment
                % property, and store them in a set of compartment properties.

                compartment_table_len = size(compartment_table, 1);

                if contains(field_name, "_")

                    % Strip the compartment tag prefix from the possible
                    % compartment field.

                    split_fi_name = string(strsplit(field_name, "_"));

                    prefix = split_fi_name(1);

                    suffix_parts = split_fi_name(2:end);

                    suffix = string(join(suffix_parts, "_"));

                    if ismember(prefix, zef.compartment_tags) ...
                    && ismember(suffix, zef_as_class.VolumeCompartment.VOLUME_COMPARTMENT_FIELD_NAMES)

                        % Check if compartment name is already in the
                        % compartment table.

                        name_col = [compartment_table{:,1}];

                        if isempty(name_col)

                            compartment_ind = compartment_table_len + 1;

                            compartment_exists = false;

                        else

                            compartment_ind = find(ismember(name_col, prefix));

                            compartment_exists = any(compartment_ind);

                            if ~ compartment_exists

                                compartment_ind = compartment_table_len + 1;

                            end

                        end

                        if compartment_exists

                            compartment_table{compartment_ind, 2}.(suffix) = field_value;

                        else

                            compartment_table{compartment_ind, 1} = prefix;

                            compartment_table{compartment_ind, 2}.(suffix) = struct;

                        end

                    end % if

                end % if

            end % for

            % Then iterate over the accumulated compartment info to generate
            % the compartment objects.

            n_of_compartment_structs = size(compartment_table, 1);

            compartment_cells = cell(0);

            for ii = 1 : n_of_compartment_structs

                a_struct = compartment_table{ii,2};

                compartment_cells{ii} = zef_as_class.VolumeCompartment(a_struct);

            end

            % Turn cell array into vector.

            compartments = [compartment_cells{:}];

        end % function

    end % methods (Static)

end % classdef

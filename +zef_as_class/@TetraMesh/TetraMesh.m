classdef TetraMesh < zef_as_class.FEMesh

    %
    % TetraMesh < FEMesh
    %
    % Describes a finite element mesh, whose elements are tetrahedral.
    %

    properties

        nodes = []

        elements = []

    end

    methods

        function self = TetraMesh(nodes, tetra)

            %
            % TetraMesh
            %
            % The constructor of this class.
            %
            % Input
            %
            % - nodes
            %
            %   Nodes with 3D Cartesian coordinates.
            %
            % - tetra
            %
            %   An M × 4 array of node indices.
            %
            % Output
            %
            % - self
            %
            %   An instance of this class.
            %

            self.nodes = nodes;

            self.elements = tetra;

        end

        function self = set.elements(self, tetra)

            %
            % set.elements(self, tetra)
            %
            % A setter method to enforce that the elements of this class are
            % indeed tetrahedral, that is their rows contain exactly 4 node
            % indices.
            %

            element_size = size(tetra, 2);

            if size(tetra, 2) ~= 4

                error("Tetrahedra have 4 node indices. Got " + element_size + " instead.")

            end

            self.elements = tetra;

        end % function

    end

end % classdef

function compartments = active_compartments(self)

    %
    % active_compartments
    %
    % Returns an array of compartments that are active.
    %
    % Input:
    %
    % - self
    %
    %   The instance of Zef that called this method.
    %
    % Output:
    %
    % - compartments
    %
    %   The active compartments in a vector.
    %

    arguments

        self zef_as_class.Zef

    end

    active_inds = self.active_compartment_inds();

    compartments = self.compartments(active_inds);

end

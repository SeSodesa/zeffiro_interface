function zef = zef_struct(args)

    %
    % zef_struct
    %
    % A function which defines the structure of the core zef struct of Zeffiro
    % Interface. The structure is enforced via Matlab's function argument
    % validation:
    %
    % - The names, types, validation functions and default values of the
    %   fields of zef are declared in the arguments block.
    %
    % - This function then feeds either the default or given fields into the
    %   returned zef.
    %
    % To modify a single field called f with a value of fval, call this function with
    %
    %   zef = zef_struct("f", fval)
    %
    % or starting from Matlab r2021a with
    %
    %   zef = zef_struct(f=fval)
    %
    % This removes the need of setting every field there is, if the defaults
    % are good enough.
    %

    arguments

        args.nodes (:,3) double = []

        args.tetra (:,4) double { mustBeInteger, mustBePositive } = []

    end

    zef.nodes = args.nodes;

    zef.tetra = args.tetra;

end

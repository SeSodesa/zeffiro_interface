function out_reduced_interpolation_matrix = zef_kron_reduction( ...
    in_interpolation_matrix, ...
    in_schur_complement, ...
    in_electrode_model, ...
    in_source_model ...
)

    % Documentation
    %
    % Performs a Kron reduction on a given interpolation matrix G. In
    % practical terms, this might mean simplifying the structure of the
    % interpolation matrix, such that the remaining graph relation is still
    % equivalent to the original one.
    %
    % Input:
    %
    % - in_interpolation_matrix: the matrix being reduced.
    %
    % - in_schur_compement: this is applied to in_interpolation_matrix to
    %
    % - in_electrode_model: if this is not 'CEM' but 'PEM', the matrix is
    %   returned as-is. Has to be one of these options.
    %
    % Output:
    %
    % - out_reduced_interpolation_matrix: the reduced (Complete Electrode
    %   Model) or unreduced (Partial Elecrode Model) interpolation matrix.

    arguments
        in_interpolation_matrix
        in_schur_complement
        in_electrode_model { mustBeText, mustBeMember(in_electrode_model, {'CEM', 'PEM'}) }
    end

    out_reduced_interpolation_matrix = in_interpolation_matrix;

    schur_size = size(in_schur_complement);

    if isequal(in_electrode_model,'CEM')

        inv_schur_complement = in_schur_complement \ eye(schur_size);

        out_reduced_interpolation_matrix = ...
            inv_schur_complement ...
            * ...
            in_interpolation_matrix ...
        ;

    end
end

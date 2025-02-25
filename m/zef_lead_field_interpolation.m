function [G, dipole_locations] =  zef_lead_field_interpolation( ...
    p_nodes, ...
    p_tetrahedra, ...
    p_brain_inds, ...
    p_source_model, ...
    p_intended_source_inds, ...
    p_nearest_neighbour_inds, ...
    p_optimization_system_type, ...
    p_regparam ...
)

    % Chooses the interpolation method based on given source_model and
    % constructs the interpolation matrix G and interpolation (dipole)
    % positions. Errors in case of unknown source model.
    %
    % Input (must be given even if not used, 'cause nobody wants to spend time
    % parsing the varargin tuple):
    %
    % - p_nodes (common)
    % - p_tetrahedra (common)
    % - p_brain_inds (common)
    % - p_source_model (common).
    % - p_intended_source_inds (Whitney and Hdiv)
    % - p_nearest_neighbour_inds (all continuous variants of the source models)
    % - p_optimization_system_type
    % - p_regparam (St Venant)
    %
    % Output:
    %
    % - Interpolation matrix G
    % - source_positions

    arguments
        p_nodes (:,3) double {mustBeNonNan}
        p_tetrahedra (:,4) double {mustBeInteger, mustBePositive}
        p_brain_inds (:,1) double {mustBeInteger, mustBePositive}
        p_source_model
        p_intended_source_inds (:,1) double {mustBeInteger, mustBePositive}
        p_nearest_neighbour_inds (:,1) double {mustBeInteger, mustBePositive}
        p_optimization_system_type { ...
            mustBeText, ...
            mustBeMember(p_optimization_system_type,{'pbo','mpo','none'}) ...
        }
        p_regparam (1,1) double
    end

    switch ZefSourceModel.from(p_source_model)

        case { ZefSourceModel.Whitney, ZefSourceModel.ContinuousWhitney }

            [G, dipole_locations] = zef_whitney_interpolation( ...
                p_nodes, ...
                p_tetrahedra, ...
                p_brain_inds, ...
                p_intended_source_inds, ...
                p_nearest_neighbour_inds, ...
                p_optimization_system_type ...
            );

        case { ZefSourceModel.Hdiv, ZefSourceModel.ContinuousHdiv }

            [G, dipole_locations] = zef_hdiv_interpolation( ...
                p_nodes, ...
                p_tetrahedra, ...
                p_brain_inds, ...
                p_intended_source_inds, ...
                p_nearest_neighbour_inds, ...
                p_optimization_system_type ...
            );

        case { ZefSourceModel.StVenant, ZefSourceModel.ContinuousStVenant }

            [G, dipole_locations] = zef_st_venant_interpolation( ...
                p_nodes, ...
                p_tetrahedra, ...
                p_brain_inds, ...
                p_intended_source_inds, ...
                p_nearest_neighbour_inds, ...
                p_regparam ...
            );

        otherwise

            error("Unknown source model.");

    end

end

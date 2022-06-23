function nearest_list = zef_nearest_points( ...
    points, ...
    neighbour_points, ...
    quantity, ...
    quantity_interpretation ...
)

    % Documentation
    %
    % Find the nearest neighbours either within a single set of points, or
    % between two point clouds. In the first case, simply pass in the same set
    % of points twice. The notion of "nearness" must be specified as a
    % quantity–interpretation pair.
    %
    % Input:
    %
    % - points: one set of points from which we are looking nearest
    %   neighbours from.
    %
    % - neighbour_points: the point cloud we are looking the neighbours of
    %   points from.
    %
    % - quantity: a real number whose interpretation is defined by the below
    %   argument. Might be restricted to integers, for example.
    %
    % - quantity_interpretation: an interpretation for the quantity
    %   parameter. Can be either 'single', 'count' or 'range'.
    %
    % Output
    %
    % - nearest_list: a list of indices of the nearest neighbours between the
    %   two given point clouds.

    arguments
        points (:, 3) double
        neighbour_points (:, 3) double
        quantity (1,1) double { mustBeReal, mustBeNonnegative }
        quantity_interpretation { mustBeText, mustBeMember(quantity_interpretation, ['single', 'count', 'range']) }
    end

    if strcmp(quantity_interpretation, 'count')

        MdlKDT = KDTreeSearcher(neighbour_points);

        if ~ isinteger(quantity) | quantity < 0
            error('Given quantity must be a positive integer when its interpretation is the number of neighbours');
        end

        nearest_list = knnsearch(MdlKDT, points, 'K', quantity);

    elseif strcmp(quantity_interpretation, 'range')

        nearest_list = rangesearch(points, neighbour_points, quantity);

    elseif strcmp(quantity_interpretation, 'single')

        MdlKDT = KDTreeSearcher(neighbour_points);

        nearest_list = knnsearch(MdlKDT, points);

    else

        error('Undefined quantity interpretation. Must be either ''single'', ''count'' or ''range''.');

    end % if

end % function

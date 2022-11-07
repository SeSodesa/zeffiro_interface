function [sr, sc] = zef_ES_objective_function(varargin)
%% Check for zef
if nargin == 0
    zef = evalin('base','zef');
end
if nargin == 1
    zef = varargin{1};
end
%% pre-allocation of the objective functions and threshold parameters
vec     = zef_ES_table(zef);
obj1    = zef.ES_obj_fun;
obj2    = zef.ES_obj_fun_2;
AT      = zef.ES_acceptable_threshold;
TT      = zef.ES_threshold_condition;

metacriteria_names  = vec.Properties.VariableNames;
metacriteria_minmax = vec.Properties.VariableDescriptions;
metacriteria_objfun = metacriteria_names(~strcmpi(metacriteria_minmax, 'none'));

obj_fun_1 = vec.(metacriteria_objfun{obj1}){1};
obj_fun_2 = vec.(metacriteria_objfun{obj2}){1};
if isempty(obj_fun_1) || isempty(obj_fun_2)
    error('ZI error: No data has been calculated yet.')
end
%% 'sweet spot' indexing based on objective
if isequal(obj1, obj2)
    if     strcmpi(metacriteria_minmax(obj1), 'minimum')
        [~, Idx] = min(abs(obj_fun_1),[],'all');
    elseif strcmpi(metacriteria_minmax(obj1), 'maximum')
        [~, Idx] = max(abs(obj_fun_1),[],'all');
    end
    [sr, sc] = ind2sub(size(obj_fun_1), Idx);

else

    if isequal(TT, 1)
        AT = AT*max((obj_fun_1(:)));
    end

    if     strcmpi(metacriteria_minmax(obj1), 'minimum')
        [Idx] = find(abs(obj_fun_1(:)) <= AT);
        if isempty(Idx)
            [~, Idx] = min(abs(obj_fun_1(:)));
        end
    elseif strcmpi(metacriteria_minmax(obj1), 'maximum')
        [Idx] = find(abs(obj_fun_1(:)) >= AT);
        if isempty(Idx)
            [~, Idx] = max(abs(obj_fun_1(:)));
        end
    end

    if     strcmpi(metacriteria_minmax(obj2), 'minimum')
        [~, Idx_2] = min(abs(obj_fun_2(Idx)));
    elseif strcmpi(metacriteria_minmax(obj2), 'maximum')
        [~, Idx_2] = max(abs(obj_fun_2(Idx)));
    end
    [sr, sc] = ind2sub(size(obj_fun_2), Idx(Idx_2));
end
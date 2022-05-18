function source_model = zef_source_model_from(input)

    % A function that creates an instance of the enumeration ZefSourceModel
    % based on given input.

    % Initialize value as invalid

    source_model = ZefSourceModel.Error;

    % Check that we receive input.

    if ~ nargin == 1
        warning('I accept exactly one input.')
        return
    end

    % Function is idempotent with respect to given ZefSourceModel variants.

    if isenum(input)
        switch input
            case ZefSourceModel.Error
                source_model = input;
            case ZefSourceModel.Whitney
                source_model = input;
            case ZefSourceModel.Hdiv
                source_model = input;
            otherwise
                warning("I received an enum that is not a ZefSourceModel. Setting erraneous return value.")
                source_model = ZefSourceModel.Error;
        end

        return

    end

    KNOWN_INTEGERS = [1, 2];

    % Check for valid inputs.

    if isreal(input)

        if input == 1
            source_model = ZefSourceModel.Whitney;
        elseif input == 2
            source_model = ZefSourceModel.Hdiv;
        else
            % Do nothing
        end

    end

    if ischar(input) || isstring(input)

        if strcmp(input, '1')
            source_model = ZefSourceModel.Whitney;
        elseif strcmp(input, '2')
            source_model = ZefSourceModel.Hdiv;
        else
            % Do nothing
        end

    end

    if source_model == ZefSourceModel.Error
        warning(strcat( ...
            "The source model was not set properly. Needs to be initialized with one of the known values ", ...
            "{" , num2str(KNOWN_INTEGERS), "}." ...
        ));
    end

end

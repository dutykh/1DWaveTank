function val = get_param(params_struct, field_name, default_value)

    %GET_PARAM Retrieves a parameter value from a structure, using a default if not found.
    %   val = GET_PARAM(params_struct, field_name, default_value) checks if
    %   'field_name' exists in the structure 'params_struct'. If it exists,
    %   its value is returned. Otherwise, 'default_value' is returned.
    %
    %   Inputs:
    %       params_struct - Structure containing parameters.
    %       field_name    - String name of the parameter field to retrieve.
    %       default_value - Value to return if the field does not exist.
    %
    %   Outputs:
    %       val           - The retrieved parameter value or the default value.

    if isfield(params_struct, field_name)
        val = params_struct.(field_name);
    else
        val = default_value;
        % Optional: Display a message when using a default value
        % fprintf('Info: Parameter ''%s'' not found in config structure. Using default value: %s\n', ...
        %         field_name, mat2str(default_value));
    end

end
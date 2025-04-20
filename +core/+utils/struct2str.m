function s = struct2str(st)

    %STRUCT2STR Convert a simple scalar structure to a comma-separated string.
    %
    %   S = STRUCT2STR(ST) converts the scalar structure ST into a string S
    %   where fields and values are represented as 'fieldname1=value1, fieldname2=value2'.
    %   Values are converted to strings using mat2str.
    %
    %   Example:
    %       params.a = 0.05;
    %       params.T = 2.0;
    %       str_rep = core.utils.struct2str(params); % Returns 'a=0.05, T=2'

    f = fieldnames(st);
    vals = cellfun(@(x) mat2str(st.(x)), f, 'UniformOutput', false);
    s_parts = cellfun(@(name, val) [name '=' val], f, vals, 'UniformOutput', false);
    s = strjoin(s_parts, ', ');

end
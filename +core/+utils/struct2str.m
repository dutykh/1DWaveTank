function s = struct2str(st)
    f = fieldnames(st);
    vals = cellfun(@(x) mat2str(st.(x)), f, 'UniformOutput', false);
    s_parts = cellfun(@(name, val) [name '=' val], f, vals, 'UniformOutput', false);
    s = strjoin(s_parts, ', ');
end
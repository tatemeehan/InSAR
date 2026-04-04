function data = load_export_struct(filename, fieldName)
S = load(filename);
if ~isfield(S, fieldName)
    error('Field "%s" not found in %s', fieldName, filename);
end
data = S.(fieldName);
end
function CR = read_corner_reflectors(csvPath, nameLength)
% READ_CORNER_REFLECTORS Load CR table and standardize names
if nargin < 2, nameLength = 4; end
CR = readtable(csvPath);
CR.Name = cellfun(@(s) s(1:min(nameLength, numel(s))), CR.Name, 'UniformOutput', false);
end

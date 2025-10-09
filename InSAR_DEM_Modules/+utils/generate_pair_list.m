function pairList = generate_pair_list(N, mode, customPairs)
%GENERATE_PAIR_LIST Generate list of index pairs for interferometric processing
%
%   pairList = generate_pair_list(N, mode)
%   pairList = generate_pair_list(N, 'custom', customPairs)
%
%   Inputs:
%       N           - Number of trajectories
%       mode        - Pairing mode: 'sequential', 'all', or 'custom'
%       customPairs - (optional) Nx2 array of custom pairs, required if mode='custom'
%
%   Output:
%       pairList    - Mx2 array of index pairs

if nargin < 2
    mode = 'sequential';
end

switch lower(mode)
    case 'sequential'
        % (1,2), (2,3), ..., (N-1,N)
        pairList = [(1:N-1)', (2:N)'];

    case 'all'
        % All unique combinations i < j
        pairList = nchoosek(1:N, 2);

    case 'custom'
        if nargin < 3 || isempty(customPairs)
            error('Custom mode requires a customPairs input (Nx2 array).');
        end
        if size(customPairs,2) ~= 2
            error('customPairs must be an Nx2 array.');
        end
        pairList = customPairs;

    otherwise
        error('Unknown pairing mode: %s. Valid options are: sequential, all, custom.', mode);
end
end

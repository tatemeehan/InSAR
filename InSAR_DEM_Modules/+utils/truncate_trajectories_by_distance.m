function [trajOut, commonLength] = truncate_trajectories_by_distance(trajIn, maxDist)
%TRUNCATE_TRAJECTORIES_BY_DISTANCE Truncate multiple trajectories to overlapping region.
%
%   Inputs:
%       trajIn   - cell array of [N x 3] trajectories (UTM coordinates)
%       maxDist  - allowed Euclidean separation [m] at start and end
%
%   Outputs:
%       trajOut      - cell array of truncated trajectories
%       commonLength - number of points in the truncated trajectories

% Number of trajectories
nTraj = numel(trajIn);

% Extract start and end points for all trajectories
starts = cellfun(@(t) t(1,:), trajIn, 'UniformOutput', false);
ends   = cellfun(@(t) t(end,:), trajIn, 'UniformOutput', false);

% Compute Euclidean distances to the reference trajectory (first one)
startDists = cellfun(@(pt) norm(pt - starts{1}), starts);
endDists   = cellfun(@(pt) norm(pt - ends{1}), ends);

% Check which trajectories are within threshold
valid = (startDists <= maxDist) & (endDists <= maxDist);
if any(~valid)
    warning('Some trajectories exceed maxDist: they will be truncated to the shortest valid trajectory.');
end

% Find minimum length of valid trajectories
minLen = min(cellfun(@(t) size(t, 1), trajIn(valid)));

% Truncate all trajectories to the same length
trajOut = cellfun(@(t) t(1:minLen, :), trajIn, 'UniformOutput', false);
commonLength = minLen;
end

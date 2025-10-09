function [trajOut, N] = truncate_trajectories_by_overlap(trajIn)
%TRUNCATE_TRAJECTORIES_BY_OVERLAP Find and truncate trajectories to best common spatial overlap
%
% Inputs:
%   trajIn - cell array of [N x 3] position arrays
%
% Outputs:
%   trajOut - truncated trajectory segments (same length)
%   N       - number of points in the common truncated region

nTraj = numel(trajIn);
trajOut = cell(size(trajIn));

% Use first trajectory as spatial reference
refTraj = trajIn{1};

% Store each trajectory’s best start and end indices based on spatial proximity
startIdx = zeros(nTraj,1);
endIdx   = zeros(nTraj,1);

% Reference start and end points
refStart = refTraj(1,:);
refEnd   = refTraj(end,:);

for i = 1:nTraj
    thisTraj = trajIn{i};
    dStart = vecnorm(thisTraj - refStart, 2, 2);
    dEnd   = vecnorm(thisTraj - refEnd, 2, 2);
    [~, startIdx(i)] = min(dStart);
    [~, endIdx(i)]   = min(dEnd);

    if endIdx(i) <= startIdx(i)
        error('Trajectory %d: end index precedes start index.', i);
    end
end

% Now determine minimal length possible
maxStart = max(startIdx);
minEnd   = min(endIdx);
N = minEnd - maxStart + 1;

if N < 2
    error('Insufficient overlapping segment found across trajectories.');
end

% Truncate each trajectory
for i = 1:nTraj
    trajOut{i} = trajIn{i}(maxStart:minEnd, :);
end
end

% function [trajOut, N] = truncate_trajectories_by_overlap(trajIn)
% %TRUNCATE_TRAJECTORIES_BY_OVERLAP Truncate all trajectories to best-aligned overlapping region.
% %
% %   [trajOut, N] = truncate_trajectories_by_overlap(trajIn)
% %
% %   Uses nearest-point logic to find the best overlapping subsegment.
% 
% nTraj = numel(trajIn);
% trajLens = cellfun(@(T) size(T,1), trajIn);
% 
% % Use the first trajectory as a reference
% refTraj = trajIn{1};
% refStart = refTraj(1,:);
% refEnd   = refTraj(end,:);
% 
% % Initialize index arrays
% startIndices = zeros(nTraj, 1);
% endIndices   = zeros(nTraj, 1);
% 
% for i = 1:nTraj
%     thisTraj = trajIn{i};
% 
%     % Distance to reference start
%     dStart = vecnorm(thisTraj - refStart, 2, 2);
%     [~, iStart] = min(dStart);
% 
%     % Distance to reference end
%     dEnd = vecnorm(thisTraj - refEnd, 2, 2);
%     [~, iEnd] = min(dEnd);
% 
%     % Ensure ordering
%     if iEnd <= iStart
%         error('Trajectory %d: nearest end index occurs before start index.', i);
%     end
% 
%     startIndices(i) = iStart;
%     endIndices(i)   = iEnd;
% end
% 
% % Common overlapping index range
% iStartCommon = max(startIndices);
% iEndCommon   = min(endIndices);
% 
% if iEndCommon <= iStartCommon
%     error('No valid overlapping segment found among all trajectories.');
% end
% 
% N = iEndCommon - iStartCommon + 1;
% trajOut = cell(size(trajIn));
% for i = 1:nTraj
%     trajOut{i} = trajIn{i}(iStartCommon:iEndCommon, :);
% end
% end
% 
% % function [trajOut, N] = truncate_trajectories_by_overlap(trajIn, tol)
% % %TRUNCATE_TRAJECTORIES_BY_OVERLAP Truncate all trajectories to shared spatial segment
% % %
% % %   [trajOut, N] = truncate_trajectories_by_overlap(trajIn, tol)
% % %   Ensures each trajectory begins and ends within `tol` meters of shared start/end.
% % 
% % nTraj = numel(trajIn);
% % startPts = zeros(nTraj, 3);
% % endPts   = zeros(nTraj, 3);
% % lengths  = zeros(nTraj, 1);
% % 
% % for i = 1:nTraj
% %     startPts(i,:) = trajIn{i}(1,:);
% %     endPts(i,:)   = trajIn{i}(end,:);
% %     lengths(i)    = size(trajIn{i},1);
% % end
% % 
% % % Use mean location as reference (you could also try min/max combinations)
% % refStart = mean(startPts,1);
% % refEnd   = mean(endPts,1);
% % 
% % % Find start and end indices within tol
% % startIndices = zeros(nTraj,1);
% % endIndices   = zeros(nTraj,1);
% % for i = 1:nTraj
% %     dStart = vecnorm(trajIn{i} - refStart, 2, 2);
% %     dEnd   = vecnorm(trajIn{i} - refEnd, 2, 2);
% % 
% %     iStart = find(dStart < tol, 1, 'first');
% %     iEnd   = find(dEnd   < tol, 1, 'last');
% % 
% %     if isempty(iStart)
% %         error('Trajectory %d has no starting point within %.2f m of average start.', i, tol);
% %     end
% %     if isempty(iEnd)
% %         error('Trajectory %d has no ending point within %.2f m of average end.', i, tol);
% %     end
% % 
% %     startIndices(i) = iStart;
% %     endIndices(i)   = iEnd;
% % end
% % 
% % % Define common truncation range
% % iStart = max(startIndices);
% % iEnd   = min(endIndices);
% % 
% % if iEnd <= iStart
% %     error('No valid overlapping segment found among all trajectories within %.2f m.', tol);
% % end
% % 
% % % Truncate all
% % trajOut = cell(size(trajIn));
% % N = iEnd - iStart + 1;
% % for i = 1:nTraj
% %     trajOut{i} = trajIn{i}(iStart:iEnd, :);
% % end
% % end
% % 
% % % function [trajOut, N] = truncate_trajectories_by_overlap(trajIn, tol)
% % % %TRUNCATE_TRAJECTORIES_BY_OVERLAP Truncate all trajectories so they
% % % % start and end within `tol` meters of a common point.
% % % %
% % % % Inputs:
% % % %   trajIn - cell array of [N x 3] matrices (X,Y,Z)
% % % %   tol    - max Euclidean distance for start/end point similarity (meters)
% % % %
% % % % Outputs:
% % % %   trajOut - truncated trajectories, same length
% % % %   N       - new common length
% % % 
% % % nTraj = numel(trajIn);
% % % startPts = zeros(nTraj, 3);
% % % endPts   = zeros(nTraj, 3);
% % % lengths  = zeros(nTraj, 1);
% % % 
% % % for i = 1:nTraj
% % %     startPts(i,:) = trajIn{i}(1,:);
% % %     endPts(i,:)   = trajIn{i}(end,:);
% % %     lengths(i)    = size(trajIn{i},1);
% % % end
% % % 
% % % % Find reference start and end (mean for now)
% % % refStart = mean(startPts,1);
% % % refEnd   = mean(endPts,1);
% % % 
% % % % Find the earliest index in each traj where it is within tol of refStart
% % % startIndices = zeros(nTraj,1);
% % % endIndices   = zeros(nTraj,1);
% % % for i = 1:nTraj
% % %     dStart = vecnorm(trajIn{i} - refStart, 2, 2);
% % %     dEnd   = vecnorm(trajIn{i} - refEnd, 2, 2);
% % % 
% % %     startIndices(i) = find(dStart < tol, 1, 'first');
% % %     endIndices(i)   = find(dEnd   < tol, 1, 'last');
% % % end
% % % 
% % % % Compute the common overlapping index range
% % % iStart = max(startIndices);
% % % iEnd   = min(endIndices);
% % % 
% % % if iEnd <= iStart
% % %     error('Trajectories do not have overlapping start/end segments within %.2f m', tol);
% % % end
% % % 
% % % % Truncate all trajectories to the overlap segment
% % % trajOut = cell(size(trajIn));
% % % N = iEnd - iStart + 1;
% % % 
% % % for i = 1:nTraj
% % %     trajOut{i} = trajIn{i}(iStart:iEnd, :);
% % % end
% % % end

function [trajOut, startCoord, endCoord, truncateIndices] = truncate_trajectories_global(trajIn)
%TRUNCATE_TRAJECTORIES_GLOBAL Find common segment for all trajectories
%
% This version compares all start/end candidates against **all points**
% in every trajectory to find a globally consistent overlapping segment.

    N = numel(trajIn);
    allPoints = vertcat(trajIn{:});  % All XYZ points from all trajectories

    % STEP 1: Build distance matrix from each trajectory start point to ALL points
    candidateStart = [];
    candidateEnd = [];

    for i = 1:N
        candidateStart = [candidateStart; trajIn{i}(1,:)];      % first point
        candidateEnd   = [candidateEnd; trajIn{i}(end,:)];      % last point
    end

    % STEP 2: Find best global start and end points (the ones that are closest on average to all other trajectories)
    % Evaluate candidateStart against all trajIn
    minStartDist = inf(size(candidateStart,1),1);
    minEndDist   = inf(size(candidateEnd,1),1);

    for i = 1:size(candidateStart,1)
        totalDist = 0;
        for j = 1:N
            d = vecnorm(trajIn{j} - candidateStart(i,:), 2, 2);
            totalDist = totalDist + min(d);
        end
        minStartDist(i) = totalDist;
    end

    for i = 1:size(candidateEnd,1)
        totalDist = 0;
        for j = 1:N
            d = vecnorm(trajIn{j} - candidateEnd(i,:), 2, 2);
            totalDist = totalDist + min(d);
        end
        minEndDist(i) = totalDist;
    end

    % STEP 3: Choose global start and end with minimum total distance to all trajs
    [~, bestStartIdx] = min(minStartDist);
    [~, bestEndIdx]   = min(minEndDist);

    startCoord = candidateStart(bestStartIdx, :);
    endCoord   = candidateEnd(bestEndIdx, :);

    % STEP 4: Find truncation bounds for all trajectories
    trajOut = cell(size(trajIn));
    truncateIndices = zeros(N, 2);
    segmentLengths = zeros(N,1);

    % First pass: find nearest indices to global start/end for each trajectory
    for i = 1:N
        traj = trajIn{i};
        dStart = vecnorm(traj - startCoord, 2, 2);
        dEnd   = vecnorm(traj - endCoord,   2, 2);

        [~, iStart] = min(dStart);
        [~, iEnd]   = min(dEnd);

        if iEnd < iStart
            warning('Trajectory %d has reversed direction. Skipping.', i);
            trajOut{i} = [];
            truncateIndices(i,:) = [NaN NaN];
            segmentLengths(i) = NaN;
        else
            truncateIndices(i,:) = [iStart, iEnd];
            segmentLengths(i) = iEnd - iStart + 1;
        end
    end

    % Get minimum common length
    minLen = min(segmentLengths(~isnan(segmentLengths)));

    % Second pass: enforce equal length truncation
    for i = 1:N
        if isnan(truncateIndices(i,1))
            continue;  % Skip invalid trajectories
        end

        iStart = truncateIndices(i,1);
        iEnd   = truncateIndices(i,2);
        len = iEnd - iStart + 1;

        % Adjust to match minLen
        if len > minLen
            excess = len - minLen;
            padBefore = floor(excess / 2);
            padAfter  = ceil(excess / 2);

            iStart = iStart + padBefore;
            iEnd   = iEnd   - padAfter;
        end

        trajOut{i} = trajIn{i}(iStart:iEnd, :);
        truncateIndices(i,:) = [iStart, iEnd];  % update actual used indices
    end

end
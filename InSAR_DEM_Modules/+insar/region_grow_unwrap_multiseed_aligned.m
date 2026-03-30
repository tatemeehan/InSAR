function [unwrapped, owner, basinShift] = region_grow_unwrap_multiseed_aligned( ...
    wrapped, quality, quality_seed_thresh, quality_valid_thresh)
%REGION_GROW_UNWRAP_MULTISEED_ALIGNED
% Fast component-wise multiseed phase unwrapping with post-hoc basin
% alignment by integer multiples of 2*pi.
%
% INPUTS:
%   wrapped              - wrapped phase [rad]
%   quality              - quality/coherence-like map
%   quality_seed_thresh  - threshold used to define seed regions
%   quality_valid_thresh - minimum quality for propagation mask
%
% OUTPUTS:
%   unwrapped   - globally corrected unwrapped phase
%   owner       - basin/seed label for each pixel
%   basinShift  - integer 2*pi shift applied to each basin label
%
% NOTES:
%   1) Local unwrapping is done independently inside connected valid regions.
%   2) Multiple seed basins are allowed within a connected component.
%   3) After local unwrap, adjacent basins are aligned using robust integer
%      offset estimates across basin boundaries.
%   4) Truly isolated basins cannot be aligned unless some boundary relation
%      exists.

if nargin < 3 || isempty(quality_seed_thresh)
    quality_seed_thresh = 0.75;
end
if nargin < 4 || isempty(quality_valid_thresh)
    quality_valid_thresh = 0.0;
end

[rows, cols] = size(wrapped);
unwrapped = nan(rows, cols);
owner     = zeros(rows, cols, 'uint32');

% Valid propagation domain
valid = ~isnan(wrapped) & ~isnan(quality) & (quality >= quality_valid_thresh);

% Connected valid components
ccValid = bwconncomp(valid, 4);

nextBasinID = uint32(0);

for c = 1:ccValid.NumObjects
    compIdx = ccValid.PixelIdxList{c};
    compMask = false(rows, cols);
    compMask(compIdx) = true;

    % Seed regions inside this component
    seedMask = compMask & (quality >= quality_seed_thresh);
    ccSeed = bwconncomp(seedMask, 4);

    % Fallback: if no seed regions, choose best quality pixel in component
    if ccSeed.NumObjects == 0
        [~, kBest] = max(quality(compIdx));
        seedList = compIdx(kBest);
        seedIDs  = nextBasinID + 1;
        nextBasinID = seedIDs;
    else
        seedList = zeros(ccSeed.NumObjects,1);
        seedIDs  = zeros(ccSeed.NumObjects,1,'uint32');
        for s = 1:ccSeed.NumObjects
            region = ccSeed.PixelIdxList{s};
            [~, kBest] = max(quality(region));
            seedList(s) = region(kBest);
            nextBasinID = nextBasinID + 1;
            seedIDs(s) = nextBasinID;
        end
    end

    visited = false(rows, cols);

    % Queue stores [linearIndex, priority, basinID]
    % priority = -quality so smaller is better
    queue = zeros(numel(compIdx), 3);
    qsize = 0;

    % Initialize seeds
    for s = 1:numel(seedList)
        idx = seedList(s);
        if visited(idx), continue; end
        visited(idx)   = true;
        unwrapped(idx) = wrapped(idx);
        owner(idx)     = seedIDs(s);

        qsize = qsize + 1;
        queue(qsize,:) = [idx, -quality(idx), double(seedIDs(s))];
    end

    % Local multiseed growth within this connected component
    while qsize > 0
        [~, minIdx] = min(queue(1:qsize, 2));
        currentIdx  = queue(minIdx, 1);
        currentID   = uint32(queue(minIdx, 3));

        queue(minIdx,:) = queue(qsize,:);
        qsize = qsize - 1;

        [i, j] = ind2sub([rows, cols], currentIdx);
        currentWrapped   = wrapped(currentIdx);
        currentUnwrapped = unwrapped(currentIdx);

        % --- up
        if i > 1
            nidx = currentIdx - 1;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(wrapped(nidx) - currentWrapped);
                unwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx)     = currentID;
                visited(nidx)   = true;

                qsize = qsize + 1;
                queue(qsize,:) = [nidx, -quality(nidx), double(currentID)];
            end
        end

        % --- down
        if i < rows
            nidx = currentIdx + 1;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(wrapped(nidx) - currentWrapped);
                unwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx)     = currentID;
                visited(nidx)   = true;

                qsize = qsize + 1;
                queue(qsize,:) = [nidx, -quality(nidx), double(currentID)];
            end
        end

        % --- left
        if j > 1
            nidx = currentIdx - rows;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(wrapped(nidx) - currentWrapped);
                unwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx)     = currentID;
                visited(nidx)   = true;

                qsize = qsize + 1;
                queue(qsize,:) = [nidx, -quality(nidx), double(currentID)];
            end
        end

        % --- right
        if j < cols
            nidx = currentIdx + rows;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(wrapped(nidx) - currentWrapped);
                unwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx)     = currentID;
                visited(nidx)   = true;

                qsize = qsize + 1;
                queue(qsize,:) = [nidx, -quality(nidx), double(currentID)];
            end
        end
    end
end

% -------------------------------------------------------------------------
% Post-hoc basin alignment: solve integer 2*pi offsets between adjacent basins
% -------------------------------------------------------------------------
numBasins = double(max(owner(:)));
basinShift = zeros(numBasins,1);   % integer shifts in units of 2*pi

if numBasins > 1
    edges = build_basin_edges(owner, wrapped, unwrapped, quality);

    if ~isempty(edges.a)
        % Choose anchor = basin with largest area (or highest quality could also work)
        basinCounts = accumarray(double(owner(owner>0)), 1, [numBasins 1]);
        [~, anchor] = max(basinCounts);

        % BFS / graph traversal
        solved = false(numBasins,1);
        solved(anchor) = true;
        basinShift(anchor) = 0;

        changed = true;
        while changed
            changed = false;

            for e = 1:numel(edges.a)
                a = edges.a(e);
                b = edges.b(e);
                d = edges.delta(e);  % shift(b) - shift(a) = d

                if solved(a) && ~solved(b)
                    basinShift(b) = basinShift(a) + d;
                    solved(b) = true;
                    changed = true;
                elseif solved(b) && ~solved(a)
                    basinShift(a) = basinShift(b) - d;
                    solved(a) = true;
                    changed = true;
                end
            end
        end

        % Apply solved shifts
        for k = 1:numBasins
            if solved(k)
                mask = (owner == k);
                unwrapped(mask) = unwrapped(mask) + 2*pi*basinShift(k);
            end
        end
    end
end

end

% =========================================================================
function edges = build_basin_edges(owner, wrapped, unwrapped, quality)
% Build robust integer offset constraints between adjacent seed basins.

[rows, cols] = size(owner);

pairA = [];
pairB = [];
pairK = [];
pairW = [];

% --- Vertical neighbor pairs: (i,j) <-> (i+1,j)
oa = owner(1:end-1, :);
ob = owner(2:end,   :);

mask = oa > 0 & ob > 0 & oa ~= ob & ...
       ~isnan(wrapped(1:end-1,:)) & ~isnan(wrapped(2:end,:)) & ...
       ~isnan(unwrapped(1:end-1,:)) & ~isnan(unwrapped(2:end,:));

if any(mask(:))
    wa = wrapped(1:end-1,:);
    wb = wrapped(2:end,:);
    ua = unwrapped(1:end-1,:);
    ub = unwrapped(2:end,:);
    qa = quality(1:end-1,:);
    qb = quality(2:end,:);

    dphi = wrap_to_pi_local(wb(mask) - wa(mask));
    kraw = (ua(mask) + dphi - ub(mask)) / (2*pi);
    kest = round(kraw);
    wgt  = 0.5*(qa(mask) + qb(mask));

    a = double(oa(mask));
    b = double(ob(mask));

    % canonical ordering
    flipMask = a > b;
    tmp = a(flipMask); a(flipMask) = b(flipMask); b(flipMask) = tmp;
    kest(flipMask) = -kest(flipMask);

    pairA = [pairA; a];
    pairB = [pairB; b];
    pairK = [pairK; kest];
    pairW = [pairW; wgt];
end

% --- Horizontal neighbor pairs: (i,j) <-> (i,j+1)
oa = owner(:, 1:end-1);
ob = owner(:, 2:end);

mask = oa > 0 & ob > 0 & oa ~= ob & ...
       ~isnan(wrapped(:,1:end-1)) & ~isnan(wrapped(:,2:end)) & ...
       ~isnan(unwrapped(:,1:end-1)) & ~isnan(unwrapped(:,2:end));

if any(mask(:))
    wa = wrapped(:,1:end-1);
    wb = wrapped(:,2:end);
    ua = unwrapped(:,1:end-1);
    ub = unwrapped(:,2:end);
    qa = quality(:,1:end-1);
    qb = quality(:,2:end);

    dphi = wrap_to_pi_local(wb(mask) - wa(mask));
    kraw = (ua(mask) + dphi - ub(mask)) / (2*pi);
    kest = round(kraw);
    wgt  = 0.5*(qa(mask) + qb(mask));

    a = double(oa(mask));
    b = double(ob(mask));

    flipMask = a > b;
    tmp = a(flipMask); a(flipMask) = b(flipMask); b(flipMask) = tmp;
    kest(flipMask) = -kest(flipMask);

    pairA = [pairA; a];
    pairB = [pairB; b];
    pairK = [pairK; kest];
    pairW = [pairW; wgt];
end

% Consolidate repeated basin-pair observations using weighted mode
edges.a = [];
edges.b = [];
edges.delta = [];

if isempty(pairA)
    return;
end

pairs = [pairA pairB];
[uniquePairs, ~, ic] = unique(pairs, 'rows');

for p = 1:size(uniquePairs,1)
    sel = (ic == p);
    kvals = pairK(sel);
    wvals = pairW(sel);

    % reject very weak or tiny supports
    if numel(kvals) < 3
        continue;
    end

    d = weighted_mode_integer(kvals, wvals);

    edges.a(end+1,1) = uniquePairs(p,1); %#ok<AGROW>
    edges.b(end+1,1) = uniquePairs(p,2); %#ok<AGROW>
    edges.delta(end+1,1) = d;            %#ok<AGROW>
end

end

% =========================================================================
function m = weighted_mode_integer(kvals, wvals)
% Weighted mode for integer-valued offsets.
uk = unique(kvals(:));
score = zeros(size(uk));
for i = 1:numel(uk)
    score(i) = sum(wvals(kvals == uk(i)));
end
[~, idx] = max(score);
m = uk(idx);
end

% =========================================================================
function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
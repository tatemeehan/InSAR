function [phzUnwrapped, unwrapMeta] = region_grow_unwrap_bucket( ...
    phzWrapped, quality, mask, unwrapOpts)
%REGION_GROW_UNWRAP_BUCKET
% Fast multiseed quality-guided phase unwrapping using a bucket queue.
%
% This function performs LOCAL unwrapping only. It does not try to resolve
% basin-to-basin 2*pi offsets globally. That is handled later by
% ALIGN_UNWRAP_BASINS_WITH_REFERENCE.
%
% INPUTS
%   phzWrapped  - wrapped phase [rad]
%   quality     - coherence / quality map [0..1]
%   mask        - valid processing mask
%   unwrapOpts  - options struct
%
% Required/optional fields used from unwrapOpts
%   .bucketBins         (default 256)
%   .seedQualityThresh  (default 0.75)
%   .validQualityThresh (default = qualityThresh if present, else 0.3)
%   .minSeedRegionSize  (default 1)
%
% OUTPUTS
%   phzUnwrapped - locally unwrapped phase
%   unwrapMeta   - metadata struct

if nargin < 4
    unwrapOpts = struct();
end

if ~isfield(unwrapOpts,'bucketBins') || isempty(unwrapOpts.bucketBins)
    unwrapOpts.bucketBins = 256;
end
if ~isfield(unwrapOpts,'seedQualityThresh') || isempty(unwrapOpts.seedQualityThresh)
    unwrapOpts.seedQualityThresh = 0.75;
end
if ~isfield(unwrapOpts,'validQualityThresh') || isempty(unwrapOpts.validQualityThresh)
    if isfield(unwrapOpts,'qualityThresh') && ~isempty(unwrapOpts.qualityThresh)
        unwrapOpts.validQualityThresh = unwrapOpts.qualityThresh;
    else
        unwrapOpts.validQualityThresh = 0.3;
    end
end
if ~isfield(unwrapOpts,'minSeedRegionSize') || isempty(unwrapOpts.minSeedRegionSize)
    unwrapOpts.minSeedRegionSize = 1;
end

numBins = double(unwrapOpts.bucketBins);

[rows, cols] = size(phzWrapped);

phzUnwrapped = nan(rows, cols, 'like', phzWrapped);
owner        = zeros(rows, cols, 'uint32');

validMask = logical(mask) & isfinite(phzWrapped) & isfinite(quality) & ...
            (quality >= unwrapOpts.validQualityThresh);

ccValid = bwconncomp(validMask, 4);

nextBasinID = uint32(0);
seedIdxAll = [];
seedBasinIDAll = [];

maxQueueLenGlobal = 0;
pushCountGlobal   = 0;
popCountGlobal    = 0;

for c = 1:ccValid.NumObjects
    compIdx = ccValid.PixelIdxList{c};
    if isempty(compIdx)
        continue;
    end

    compMask = false(rows, cols);
    compMask(compIdx) = true;

    seedMask = compMask & (quality >= unwrapOpts.seedQualityThresh);

    if unwrapOpts.minSeedRegionSize > 1
        seedMask = bwareaopen(seedMask, unwrapOpts.minSeedRegionSize, 4);
    end

    ccSeed = bwconncomp(seedMask, 4);

    if ccSeed.NumObjects == 0
        [~, kBest] = max(quality(compIdx));
        seedList = compIdx(kBest);
        nextBasinID = nextBasinID + 1;
        seedIDs = nextBasinID;
    else
        seedList = zeros(ccSeed.NumObjects, 1);
        seedIDs  = zeros(ccSeed.NumObjects, 1, 'uint32');

        for s = 1:ccSeed.NumObjects
            region = ccSeed.PixelIdxList{s};
            [~, kBest] = max(quality(region));
            seedList(s) = region(kBest);

            nextBasinID = nextBasinID + 1;
            seedIDs(s) = nextBasinID;
        end
    end

    seedIdxAll     = [seedIdxAll; seedList(:)]; %#ok<AGROW>
    seedBasinIDAll = [seedBasinIDAll; seedIDs(:)]; %#ok<AGROW>

    visited = false(rows, cols);

    % Queue state
    queue = init_bucket_queue(numBins);

    localPushCount = 0;
    localPopCount  = 0;
    localMaxQueue  = 0;
    localQueueLen  = 0;

    % Initialize seeds
    for s = 1:numel(seedList)
        idx = seedList(s);
        basinID = seedIDs(s);

        if visited(idx)
            continue;
        end

        visited(idx) = true;
        phzUnwrapped(idx) = phzWrapped(idx);
        owner(idx) = basinID;

        [queue, localPushCount, localQueueLen, localMaxQueue] = push_bucket_queue( ...
            queue, idx, basinID, quality(idx), localPushCount, localQueueLen, localMaxQueue);
    end

    % Grow
    while true
        [queue, currentIdx, currentBasinID, ok, localPopCount, localQueueLen] = ...
            pop_bucket_queue(queue, localPopCount, localQueueLen);

        if ~ok
            break;
        end

        [i, j] = ind2sub([rows, cols], double(currentIdx));
        currentWrapped   = phzWrapped(currentIdx);
        currentUnwrapped = phzUnwrapped(currentIdx);

        % up
        if i > 1
            nidx = currentIdx - 1;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(phzWrapped(nidx) - currentWrapped);
                phzUnwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx) = currentBasinID;
                visited(nidx) = true;

                [queue, localPushCount, localQueueLen, localMaxQueue] = push_bucket_queue( ...
                    queue, nidx, currentBasinID, quality(nidx), ...
                    localPushCount, localQueueLen, localMaxQueue);
            end
        end

        % down
        if i < rows
            nidx = currentIdx + 1;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(phzWrapped(nidx) - currentWrapped);
                phzUnwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx) = currentBasinID;
                visited(nidx) = true;

                [queue, localPushCount, localQueueLen, localMaxQueue] = push_bucket_queue( ...
                    queue, nidx, currentBasinID, quality(nidx), ...
                    localPushCount, localQueueLen, localMaxQueue);
            end
        end

        % left
        if j > 1
            nidx = currentIdx - rows;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(phzWrapped(nidx) - currentWrapped);
                phzUnwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx) = currentBasinID;
                visited(nidx) = true;

                [queue, localPushCount, localQueueLen, localMaxQueue] = push_bucket_queue( ...
                    queue, nidx, currentBasinID, quality(nidx), ...
                    localPushCount, localQueueLen, localMaxQueue);
            end
        end

        % right
        if j < cols
            nidx = currentIdx + rows;
            if compMask(nidx) && ~visited(nidx)
                dphi = wrap_to_pi_local(phzWrapped(nidx) - currentWrapped);
                phzUnwrapped(nidx) = currentUnwrapped + dphi;
                owner(nidx) = currentBasinID;
                visited(nidx) = true;

                [queue, localPushCount, localQueueLen, localMaxQueue] = push_bucket_queue( ...
                    queue, nidx, currentBasinID, quality(nidx), ...
                    localPushCount, localQueueLen, localMaxQueue);
            end
        end
    end

    maxQueueLenGlobal = max(maxQueueLenGlobal, localMaxQueue);
    pushCountGlobal   = pushCountGlobal + localPushCount;
    popCountGlobal    = popCountGlobal + localPopCount;
end

numBasins = double(max(owner(:)));

basinArea = zeros(numBasins,1);
basinMeanQuality = nan(numBasins,1);

for b = 1:numBasins
    m = (owner == b);
    basinArea(b) = nnz(m);
    if basinArea(b) > 0
        basinMeanQuality(b) = mean(double(quality(m)), 'omitnan');
    end
end

unwrapMeta = struct();
unwrapMeta.owner = owner;
unwrapMeta.seedIdx = seedIdxAll;
unwrapMeta.seedBasinID = seedBasinIDAll;
unwrapMeta.basinArea = basinArea;
unwrapMeta.basinMeanQuality = basinMeanQuality;
unwrapMeta.validMask = validMask;
unwrapMeta.phzLocal = phzUnwrapped;
unwrapMeta.queueStats = struct( ...
    'maxQueueLen', maxQueueLenGlobal, ...
    'pushCount',   pushCountGlobal, ...
    'popCount',    popCountGlobal);

end

% =========================================================================
function queue = init_bucket_queue(numBins)

queue = struct();
queue.numBins = numBins;
queue.bucketIdx = cell(numBins,1);
queue.bucketBID = cell(numBins,1);
queue.bucketCnt = zeros(numBins,1,'uint32');
queue.maxNonEmpty = uint32(0);

for b = 1:numBins
    queue.bucketIdx{b} = zeros(128,1,'uint32');
    queue.bucketBID{b} = zeros(128,1,'uint32');
end

end

% =========================================================================
function [queue, pushCount, queueLen, maxQueueLen] = push_bucket_queue( ...
    queue, idx, basinID, qval, pushCount, queueLen, maxQueueLen)

qval = double(qval);
if ~isfinite(qval)
    qval = 0;
end
qval = max(0, min(1, qval));

bin = 1 + floor(qval * (queue.numBins - 1));
bin = uint32(bin);

queue.bucketCnt(bin) = queue.bucketCnt(bin) + 1;
k = queue.bucketCnt(bin);

if k > numel(queue.bucketIdx{bin})
    oldN = numel(queue.bucketIdx{bin});
    queue.bucketIdx{bin}(oldN+1:2*oldN) = uint32(0);
    queue.bucketBID{bin}(oldN+1:2*oldN) = uint32(0);
end

queue.bucketIdx{bin}(k) = uint32(idx);
queue.bucketBID{bin}(k) = uint32(basinID);

if bin > queue.maxNonEmpty
    queue.maxNonEmpty = bin;
end

pushCount = pushCount + 1;
queueLen = queueLen + 1;
if queueLen > maxQueueLen
    maxQueueLen = queueLen;
end

end

% =========================================================================
function [queue, idx, basinID, ok, popCount, queueLen] = ...
    pop_bucket_queue(queue, popCount, queueLen)

idx = uint32(0);
basinID = uint32(0);
ok = false;

while queue.maxNonEmpty >= 1
    if queue.bucketCnt(queue.maxNonEmpty) > 0
        k = queue.bucketCnt(queue.maxNonEmpty);

        idx = queue.bucketIdx{queue.maxNonEmpty}(k);
        basinID = queue.bucketBID{queue.maxNonEmpty}(k);
        queue.bucketCnt(queue.maxNonEmpty) = k - 1;

        popCount = popCount + 1;
        queueLen = max(queueLen - 1, 0);
        ok = true;
        return;
    else
        queue.maxNonEmpty = queue.maxNonEmpty - 1;
    end
end

end

% =========================================================================
function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
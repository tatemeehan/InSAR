function unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% REGION_GROW_UNWRAP_MULTISEED Fast region-growing phase unwrapping with multiple seeds
%
%   unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
%
% INPUTS:
%   wrapped         - 2D wrapped phase image (radians)
%   quality         - 2D phase quality map (0 to 1)
%   quality_thresh  - Threshold for region seeding (e.g., 0.3)
%
% OUTPUT:
%   unwrapped       - 2D unwrapped phase image
%

if nargin < 3
    quality_thresh = 0.75; % phase quality threshold for seeding
end

[rows, cols] = size(wrapped);
N = rows * cols;

unwrapped = nan(rows, cols);
visited = false(rows, cols);

% Find connected components above threshold
mask = quality >= quality_thresh & ~isnan(wrapped);
cc = bwconncomp(mask, 4);
numSeeds = cc.NumObjects;

% Use linear indexing for queue (preallocate)
queue = zeros(N, 2);
offsets = [-1, 0; 1, 0; 0, -1; 0, 1];
qsize = 0;

for s = 1:numSeeds
    region = cc.PixelIdxList{s};
    [~, bestIdx] = max(quality(region));
    seedIdx = region(bestIdx);

    unwrapped(seedIdx) = wrapped(seedIdx);
    visited(seedIdx) = true;

    qsize = qsize + 1;
    queue(qsize,:) = [seedIdx, -quality(seedIdx)];
end

% Begin region-growing from all seeds
while qsize > 0
    [~, minIdx] = min(queue(1:qsize, 2));
    currentIdx = queue(minIdx, 1);
    [i, j] = ind2sub([rows, cols], currentIdx);
    queue(minIdx,:) = queue(qsize,:);
    qsize = qsize - 1;

    for n = 1:4
        ni = i + offsets(n,1);
        nj = j + offsets(n,2);
        if ni < 1 || nj < 1 || ni > rows || nj > cols
            continue;
        end
        if visited(ni, nj) || isnan(wrapped(ni, nj))
            continue;
        end
        neighborIdx = sub2ind([rows, cols], ni, nj);
        dphi = wrapToPi(wrapped(ni, nj) - wrapped(i, j));
        unwrapped(ni, nj) = unwrapped(i, j) + dphi;
        visited(ni, nj) = true;

        qsize = qsize + 1;
        queue(qsize,:) = [neighborIdx, -quality(ni, nj)];
    end
end
end
% % --- Post-process: attempt regression-based correction in untrusted regions ---
% trusted = quality >= 0.3 & isfinite(wrapped) & bwdist(isnan(wrapped)) > 3;
% 
% if any(trusted(:))
%     [X, Y] = meshgrid(1:cols, 1:rows);
%     xt = X(trusted); yt = Y(trusted); zt = unwrapped(trusted);
%     xf = X(~trusted); yf = Y(~trusted);
% 
%     % Fit surface using polyfitn or griddata
%     try
%         F = scatteredInterpolant(xt, yt, zt, 'linear', 'none');
%         unwrapped(~trusted) = F(xf, yf);
%     catch
%         warning('Surface regression failed');
%     end
% end
% end

% function unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% % REGION_GROW_UNWRAP_MULTISEED Fast region-growing phase unwrapping with multiple seeds
% %
% %   unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% %
% % INPUTS:
% %   wrapped         - 2D wrapped phase image (radians)
% %   quality         - 2D phase quality map (0 to 1)
% %   quality_thresh  - Threshold for region seeding (e.g., 0.3)
% %
% % OUTPUT:
% %   unwrapped       - 2D unwrapped phase image
% %
% 
% if nargin < 3
%     quality_thresh = 0.75; % phase quality threshold for seeding
% end
% 
% [rows, cols] = size(wrapped);
% N = rows * cols;
% 
% unwrapped = nan(rows, cols);
% visited = false(rows, cols);
% 
% % Find connected components above threshold
% mask = quality >= quality_thresh & ~isnan(wrapped);
% cc = bwconncomp(mask, 4);
% numSeeds = cc.NumObjects;
% 
% % Use linear indexing for queue (preallocate)
% queue = zeros(N, 2);
% offsets = [-1, 0; 1, 0; 0, -1; 0, 1];
% qsize = 0;
% 
% for s = 1:numSeeds
%     region = cc.PixelIdxList{s};
%     [~, bestIdx] = max(quality(region));
%     seedIdx = region(bestIdx);
% 
%     unwrapped(seedIdx) = wrapped(seedIdx);
%     visited(seedIdx) = true;
% 
%     qsize = qsize + 1;
%     queue(qsize,:) = [seedIdx, -quality(seedIdx)];
% end
% 
% % Begin region-growing from all seeds
% while qsize > 0
%     [~, minIdx] = min(queue(1:qsize, 2));
%     currentIdx = queue(minIdx, 1);
%     [i, j] = ind2sub([rows, cols], currentIdx);
%     queue(minIdx,:) = queue(qsize,:);
%     qsize = qsize - 1;
% 
%     for n = 1:4
%         ni = i + offsets(n,1);
%         nj = j + offsets(n,2);
%         if ni < 1 || nj < 1 || ni > rows || nj > cols
%             continue;
%         end
%         if visited(ni, nj) || isnan(wrapped(ni, nj))
%             continue;
%         end
%         neighborIdx = sub2ind([rows, cols], ni, nj);
%         dphi = wrapToPi(wrapped(ni, nj) - wrapped(i, j));
%         unwrapped(ni, nj) = unwrapped(i, j) + dphi;
%         visited(ni, nj) = true;
% 
%         qsize = qsize + 1;
%         queue(qsize,:) = [neighborIdx, -quality(ni, nj)];
%     end
% end
% % --- Intra-region gradient ramp correction ---
% [L, numRegions] = bwlabel(~isnan(unwrapped), 4);
% 
% for r = 1:numRegions
%     mask = (L == r);
%     [iy, ix] = find(mask);
%     regionVals = unwrapped(mask);
% 
%     % Fit a 2D plane to unwrapped phase: z = ax + by + c
%     A = [ix, iy, ones(size(ix))];
%     coeffs_unwrapped = A \ regionVals;
% 
%     % Fit a 2D plane to wrapped phase
%     regionWrapped = wrapped(mask);
%     coeffs_wrapped = A \ regionWrapped;
% 
%     % Compute gradient difference
%     gradDiff = coeffs_unwrapped(1:2) - coeffs_wrapped(1:2);
%     gradMag = norm(gradDiff);
% 
%     % If discrepancy exceeds threshold, subtract fitted bias
%     if gradMag > 0.5  % radians/pixel
%         % Rebuild ramp
%         ramp = A * [gradDiff; 0];
%         corrected = regionVals - ramp;
%         unwrapped(mask) = corrected;
%     end
% end
% end
% 
% % function unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% % % REGION_GROW_UNWRAP_MULTISEED MST-enhanced region-growing phase unwrapping with multiple seeds
% % %
% % %   unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% % %
% % % INPUTS:
% % %   wrapped         - 2D wrapped phase image (radians)
% % %   quality         - 2D phase quality map (0 to 1)
% % %   quality_thresh  - Threshold for region seeding (e.g., 0.3)
% % %
% % % OUTPUT:
% % %   unwrapped       - 2D unwrapped phase image
% % %
% % 
% % if nargin < 3
% %     quality_thresh = 0.75; % phase quality threshold for seeding
% % end
% % 
% % [rows, cols] = size(wrapped);
% % unwrapped = nan(rows, cols);
% % 
% % % Find connected components above threshold
% % mask = quality >= quality_thresh & ~isnan(wrapped);
% % cc = bwconncomp(mask, 4);
% % 
% % % Offsets for 4-connectivity
% % offsets = [-1, 0; 1, 0; 0, -1; 0, 1];
% % 
% % for s = 1:cc.NumObjects
% %     regionIdx = cc.PixelIdxList{s};
% %     [ii, jj] = ind2sub([rows, cols], regionIdx);
% %     numPix = numel(regionIdx);
% % 
% %     % Build sparse graph: nodes = region pixels
% %     G = sparse(numPix, numPix);
% %     idxMap = zeros(rows, cols);
% %     idxMap(regionIdx) = 1:numPix;
% % 
% %     for k = 1:numPix
% %         i = ii(k); j = jj(k);
% %         for n = 1:4
% %             ni = i + offsets(n,1);
% %             nj = j + offsets(n,2);
% %             if ni < 1 || nj < 1 || ni > rows || nj > cols
% %                 continue;
% %             end
% %             if idxMap(ni,nj) == 0
% %                 continue;
% %             end
% %             a = idxMap(i,j);
% %             b = idxMap(ni,nj);
% %             dphi = wrapToPi(wrapped(ni,nj) - wrapped(i,j));
% %             weight = 1 - 0.5*(quality(i,j) + quality(ni,nj));
% %             G(a,b) = weight;
% %             G(b,a) = weight;
% %         end
% %     end
% % 
% %     % Use best-quality pixel as root
% %     [~, bestIdx] = max(quality(regionIdx));
% %     Ggraph = graph(G);
% %     T = minspantree(Ggraph);
% % 
% %     % Traverse tree and unwrap
% %     unwrappedRegion = nan(numPix,1);
% %     unwrappedRegion(bestIdx) = wrapped(regionIdx(bestIdx));
% %     visited = false(numPix,1);
% %     visited(bestIdx) = true;
% %     stack = bestIdx;
% % 
% %     while ~isempty(stack)
% %         current = stack(end);
% %         stack(end) = [];
% %         nbrs = neighbors(T, current)';
% %         for nbr = nbrs
% %             if visited(nbr), continue; end
% %             dphi = wrapToPi(wrapped(regionIdx(nbr)) - wrapped(regionIdx(current)));
% %             unwrappedRegion(nbr) = unwrappedRegion(current) + dphi;
% %             visited(nbr) = true;
% %             stack(end+1) = nbr;
% %         end
% %     end
% % 
% %     unwrapped(regionIdx) = unwrappedRegion;
% % end
% % 
% % % --- Intra-region gradient ramp correction ---
% % [L, numRegions] = bwlabel(~isnan(unwrapped), 4);
% % 
% % for r = 1:numRegions
% %     mask = (L == r);
% %     [iy, ix] = find(mask);
% %     regionVals = unwrapped(mask);
% % 
% %     % Fit a 2D plane to unwrapped phase: z = ax + by + c
% %     A = [ix, iy, ones(size(ix))];
% %     coeffs_unwrapped = A \ regionVals;
% % 
% %     % Fit a 2D plane to wrapped phase
% %     regionWrapped = wrapped(mask);
% %     coeffs_wrapped = A \ regionWrapped;
% % 
% %     % Compute gradient difference
% %     gradDiff = coeffs_unwrapped(1:2) - coeffs_wrapped(1:2);
% %     gradMag = norm(gradDiff);
% % 
% %     % If discrepancy exceeds threshold, subtract fitted bias
% %     if gradMag > 0.5  % radians/pixel
% %         % Rebuild ramp
% %         ramp = A * [gradDiff; 0];
% %         corrected = regionVals - ramp;
% %         unwrapped(mask) = corrected;
% %     end
% % end
% % end
% % 
% % % function unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% % % % REGION_GROW_UNWRAP_MULTISEED MST-enhanced region-growing phase unwrapping with multiple seeds
% % % %
% % % %   unwrapped = region_grow_unwrap_multiseed(wrapped, quality, quality_thresh)
% % % %
% % % % INPUTS:
% % % %   wrapped         - 2D wrapped phase image (radians)
% % % %   quality         - 2D phase quality map (0 to 1)
% % % %   quality_thresh  - Threshold for region seeding (e.g., 0.3)
% % % %
% % % % OUTPUT:
% % % %   unwrapped       - 2D unwrapped phase image
% % % %
% % % 
% % % if nargin < 3
% % %     quality_thresh = 0.75; % phase quality threshold for seeding
% % % end
% % % 
% % % [rows, cols] = size(wrapped);
% % % unwrapped = nan(rows, cols);
% % % 
% % % % Find connected components above threshold
% % % mask = quality >= quality_thresh & ~isnan(wrapped);
% % % cc = bwconncomp(mask, 4);
% % % 
% % % % Offsets for 4-connectivity
% % % offsets = [-1, 0; 1, 0; 0, -1; 0, 1];
% % % 
% % % for s = 1:cc.NumObjects
% % %     regionIdx = cc.PixelIdxList{s};
% % %     [ii, jj] = ind2sub([rows, cols], regionIdx);
% % %     numPix = numel(regionIdx);
% % % 
% % %     % Build sparse graph: nodes = region pixels
% % %     G = sparse(numPix, numPix);
% % %     idxMap = zeros(rows, cols);
% % %     idxMap(regionIdx) = 1:numPix;
% % % 
% % %     for k = 1:numPix
% % %         i = ii(k); j = jj(k);
% % %         for n = 1:4
% % %             ni = i + offsets(n,1);
% % %             nj = j + offsets(n,2);
% % %             if ni < 1 || nj < 1 || ni > rows || nj > cols
% % %                 continue;
% % %             end
% % %             if idxMap(ni,nj) == 0
% % %                 continue;
% % %             end
% % %             a = idxMap(i,j);
% % %             b = idxMap(ni,nj);
% % %             dphi = wrapToPi(wrapped(ni,nj) - wrapped(i,j));
% % %             weight = 1 - 0.5*(quality(i,j) + quality(ni,nj));
% % %             G(a,b) = weight;
% % %             G(b,a) = weight;
% % %         end
% % %     end
% % % 
% % %     % Use best-quality pixel as root
% % %     [~, bestIdx] = max(quality(regionIdx));
% % %     % T = graphminspantree(G, bestIdx);
% % %     % tree = graph(T);
% % %     Ggraph = graph(G);
% % %     T = minspantree(Ggraph);%, 'StartVertex', bestIdx);
% % %     tree = digraph(T);  % Or use `graph(T)` if you want an undirected MST
% % % 
% % %     % Traverse tree and unwrap
% % %     unwrappedRegion = nan(numPix,1);
% % %     unwrappedRegion(bestIdx) = wrapped(regionIdx(bestIdx));
% % %     visited = false(numPix,1);
% % %     visited(bestIdx) = true;
% % %     stack = bestIdx;
% % % 
% % %     while ~isempty(stack)
% % %         current = stack(end);
% % %         stack(end) = [];
% % %         nbrs = neighbors(tree, current)';
% % %         for nbr = nbrs
% % %             if visited(nbr), continue; end
% % %             dphi = wrapToPi(wrapped(regionIdx(nbr)) - wrapped(regionIdx(current)));
% % %             unwrappedRegion(nbr) = unwrappedRegion(current) + dphi;
% % %             visited(nbr) = true;
% % %             stack(end+1) = nbr;
% % %         end
% % %     end
% % % 
% % %     unwrapped(regionIdx) = unwrappedRegion;
% % % end
% % % 
% % % % --- Intra-region gradient ramp correction ---
% % % [L, numRegions] = bwlabel(~isnan(unwrapped), 4);
% % % 
% % % for r = 1:numRegions
% % %     mask = (L == r);
% % %     [iy, ix] = find(mask);
% % %     regionVals = unwrapped(mask);
% % % 
% % %     % Fit a 2D plane to unwrapped phase: z = ax + by + c
% % %     A = [ix, iy, ones(size(ix))];
% % %     coeffs_unwrapped = A \ regionVals;
% % % 
% % %     % Fit a 2D plane to wrapped phase
% % %     regionWrapped = wrapped(mask);
% % %     coeffs_wrapped = A \ regionWrapped;
% % % 
% % %     % Compute gradient difference
% % %     gradDiff = coeffs_unwrapped(1:2) - coeffs_wrapped(1:2);
% % %     gradMag = norm(gradDiff);
% % % 
% % %     % If discrepancy exceeds threshold, subtract fitted bias
% % %     if gradMag > 0.5  % radians/pixel
% % %         % Rebuild ramp
% % %         ramp = A * [gradDiff; 0];
% % %         corrected = regionVals - ramp;
% % %         unwrapped(mask) = corrected;
% % %     end
% % % end
% % % end
% 

function unwrapped = postprocess_unwrap_regions(unwrapped, quality, varargin)
%POSTPROCESS_UNWRAP_REGIONS Enhance multi-seed unwrapped phase continuity.
%   Supports stitching across disjoint regions using:
%       - MST-based region linking (default)
%       - Phase gradient blending
%       - Region bridging (neighborhood extrapolation)
%       - Small region removal or fallback
%
%   Usage:
%       unwrapped = postprocess_unwrap_regions(unwrapped, quality, 'method', 'mst', ...)
%
%   Options (passed as name-value pairs):
%       'method'     : 'mst' (default), 'blend', 'bridge', or 'none'
%       'tinySize'   : Minimum region size to retain (default = 30)
%       'alpha'      : Weighting balance between distance and quality (default = 0.5)

% opts = struct('method', 'mst', 'tinySize', 100, 'alpha', 0.5);
if nargin > 2 && ~isempty(varargin)
    % opts = parse_inputs(opts, varargin{:});
    opts = varargin{:};
end

% Identify connected unwrap regions
mask = ~isnan(unwrapped);
cc = bwconncomp(mask, 4);

% Remove or flag tiny regions
L = labelmatrix(cc);
regionSizes = cellfun(@numel, cc.PixelIdxList);
isTiny = regionSizes < opts.tinyRegionSize;
for r = find(isTiny)
    unwrapped(cc.PixelIdxList{r}) = NaN;
end

% Recompute region labels after pruning
mask = ~isnan(unwrapped);
[L, numRegions] = bwlabel(mask, 4);

switch lower(opts.stitchMethod)
    case 'none'
        return

    case 'mst'
        % Compute quality-weighted mean phase and centroid per region
        regMeans = zeros(numRegions, 1);
        regCentroids = zeros(numRegions, 2);
        regWeights = zeros(numRegions, 1);
        for r = 1:numRegions
            idx = find(L == r);
            q = quality(idx);
            p = unwrapped(idx);
            regMeans(r) = sum(p .* q, 'omitnan') / sum(q, 'omitnan');
            regWeights(r) = mean(q, 'omitnan');
            [i, j] = ind2sub(size(unwrapped), idx);
            regCentroids(r,:) = [mean(i), mean(j)];
        end

        % Build MST using combined distance and quality weighting
        D = squareform(pdist(regCentroids));
        Q = regWeights(:);
        Q(Q == 0) = eps;  % avoid division by zero
        D = D / max(D(:));
        Q = Q / max(Q);
        W = (1 - opts.stitchAlpha) * D + opts.stitchAlpha * (1 ./ (Q + Q'));
        W(isnan(W) | isinf(W)) = max(W(:));

        G = graph(W);
        T = minspantree(G);
        edges = table2array(T.Edges);

        for e = 1:size(edges,1)
            r1 = edges(e,1);
            r2 = edges(e,2);
            dp = regMeans(r2) - regMeans(r1);
            offset = round(dp / (2*pi)) * 2*pi;
            unwrapped(L == r2) = unwrapped(L == r2) - offset;
        end
    case 'mstblend'
        % MST-blending using phase differences over boundary regions
        se = strel('disk', 1);
        smoothU = unwrapped;
        regCentroids = zeros(numRegions, 2);
        for r = 1:numRegions
            [i, j] = find(L == r);
            regCentroids(r,:) = [mean(i), mean(j)];
        end
        D = squareform(pdist(regCentroids));
        G = graph(D);
        T = minspantree(G);
        edges = table2array(T.Edges);

        for e = 1:size(edges,1)
            r1 = edges(e,1);
            r2 = edges(e,2);
            % Get phase boundary between r1 and r2
            b1 = imdilate(L == r1, se) & (L == r2);
            b2 = imdilate(L == r2, se) & (L == r1);
            p1 = unwrapped(b1);
            p2 = unwrapped(b2);
            if ~isempty(p1) && ~isempty(p2)
                offset = median(p2(:)) - median(p1(:));
                smoothU(L == r2) = smoothU(L == r2) - offset;
            end
        end
        unwrapped = smoothU;

    case 'blend'
        % Gradient blending across region boundaries
        se = strel('disk', 1);
        smoothU = unwrapped;
        for r = 1:numRegions
            maskR = (L == r);
            dilated = imdilate(maskR, se);
            boundary = dilated & ~maskR;
            neighborPhases = unwrapped(boundary);
            phaseR = unwrapped(maskR);
            neighborPhases = neighborPhases(~isnan(neighborPhases));
            phaseR = phaseR(~isnan(phaseR));
            if isempty(neighborPhases) || isempty(phaseR)
                offset = 0;
            else
                % Only compute offset if there are enough samples
                if numel(neighborPhases) > 1 && numel(phaseR) > 1
                    offset = median(neighborPhases(:)) - median(phaseR(:));
                else
                    offset = 0;
                end
            end
            smoothU(maskR) = smoothU(maskR) + offset;
        end
        unwrapped = smoothU;

    case 'bridge'
        % Simple bridging across small gaps (nearest extrapolation)
        se = strel('disk', 1);
        filled = unwrapped;
        for iter = 1:5
            nanMask = isnan(filled);
            grow = imdilate(~nanMask, se) & nanMask;
            [i, j] = find(grow);
            for k = 1:numel(i)
                ni = i(k); nj = j(k);
                local = filled(max(ni-1,1):min(ni+1,end), max(nj-1,1):min(nj+1,end));
                filled(ni,nj) = mean(local(~isnan(local)), 'omitnan');
            end
        end
        unwrapped = filled;

    otherwise
        error('Unknown postprocess method: %s', opts.stitchMethod);
end
end

function opts = parse_inputs(opts, varargin)
% Helper function for name-value pair parsing
if mod(numel(varargin),2) ~= 0
    error('Arguments must be name-value pairs.');
end
for i = 1:2:numel(varargin)
    key = (varargin{i});
    val = varargin{i+1};
    if isfield(opts, key)
        opts.(key) = val;
    else
        warning('Unknown option: %s', key);
    end
end
end

% function unwrapped = postprocess_unwrap_regions(unwrapped, quality, varargin)
% %POSTPROCESS_UNWRAP_REGIONS Enhance multi-seed unwrapped phase continuity.
% %   Supports stitching across disjoint regions using:
% %       - MST-based region linking (default)
% %       - Phase gradient blending
% %       - Region bridging (neighborhood extrapolation)
% %       - Small region removal or fallback
% %
% %   Usage:
% %       unwrapped = postprocess_unwrap_regions(unwrapped, quality, 'method', 'mst', ...)
% %
% %   Options (passed as name-value pairs):
% %       'method'     : 'mst' (default), 'blend', 'bridge', or 'none'
% %       'tinySize'   : Minimum region size to retain (default = 30)
% %       'alpha'      : Weighting balance between distance and quality (default = 0.5)
% 
% opts = struct('method', 'mst', 'tinySize', 30, 'alpha', 0.5);
% if nargin > 2 && ~isempty(varargin)
%     opts = parse_inputs(opts, varargin{:});
% end
% 
% % Identify connected unwrap regions
% mask = ~isnan(unwrapped);
% cc = bwconncomp(mask, 4);
% 
% % Remove or flag tiny regions
% L = labelmatrix(cc);
% regionSizes = cellfun(@numel, cc.PixelIdxList);
% isTiny = regionSizes < opts.tinySize;
% for r = find(isTiny)
%     unwrapped(cc.PixelIdxList{r}) = NaN;
% end
% 
% % Recompute region labels after pruning
% mask = ~isnan(unwrapped);
% [L, numRegions] = bwlabel(mask, 4);
% 
% switch lower(opts.method)
%     case 'none'
%         return
% 
%     case 'mst'
%         % Compute mean phase and location per region
%         regMeans = zeros(numRegions, 1);
%         regCentroids = zeros(numRegions, 2);
%         regWeights = zeros(numRegions, 1);
%         for r = 1:numRegions
%             idx = find(L == r);
%             regMeans(r) = mean(unwrapped(idx), 'omitnan');
%             regWeights(r) = mean(quality(idx), 'omitnan');
%             [i, j] = ind2sub(size(unwrapped), idx);
%             regCentroids(r,:) = [mean(i), mean(j)];
%         end
% 
%         % Build MST using combined distance and quality weighting
%         D = squareform(pdist(regCentroids));
%         Q = regWeights(:);
%         Q(Q == 0) = eps;  % avoid division by zero
%         D = D / max(D(:));
%         Q = Q / max(Q);
%         W = (1 - opts.alpha) * D + opts.alpha * (1 ./ (Q + Q'));
%         W(isnan(W) | isinf(W)) = max(W(:));
% 
%         G = graph(W);
%         T = minspantree(G);
%         edges = table2array(T.Edges);
% 
%         for e = 1:size(edges,1)
%             r1 = edges(e,1);
%             r2 = edges(e,2);
%             dp = regMeans(r2) - regMeans(r1);
%             offset = round(dp / (2*pi)) * 2*pi;
%             unwrapped(L == r2) = unwrapped(L == r2) - offset;
%         end
% 
%     case 'blend'
%         % Gradient blending across region boundaries
%         se = strel('disk', 1);
%         smoothU = unwrapped;
%         for r = 1:numRegions
%             maskR = (L == r);
%             dilated = imdilate(maskR, se);
%             boundary = dilated & ~maskR;
%             neighborPhases = unwrapped(boundary);
%             offset = median(neighborPhases - unwrapped(maskR), 'omitnan');
%             smoothU(maskR) = smoothU(maskR) + offset;
%         end
%         unwrapped = smoothU;
% 
%     case 'bridge'
%         % Simple bridging across small gaps (nearest extrapolation)
%         se = strel('disk', 1);
%         filled = unwrapped;
%         for iter = 1:5
%             nanMask = isnan(filled);
%             grow = imdilate(~nanMask, se) & nanMask;
%             [i, j] = find(grow);
%             for k = 1:numel(i)
%                 ni = i(k); nj = j(k);
%                 local = filled(max(ni-1,1):min(ni+1,end), max(nj-1,1):min(nj+1,end));
%                 filled(ni,nj) = mean(local(~isnan(local)), 'omitnan');
%             end
%         end
%         unwrapped = filled;
% 
%     otherwise
%         error('Unknown postprocess method: %s', opts.method);
% end
% end
% 
% function opts = parse_inputs(opts, varargin)
% % Helper function for name-value pair parsing
% if mod(numel(varargin),2) ~= 0
%     error('Arguments must be name-value pairs.');
% end
% for i = 1:2:numel(varargin)
%     key = lower(varargin{i});
%     val = varargin{i+1};
%     if isfield(opts, key)
%         opts.(key) = val;
%     else
%         warning('Unknown option: %s', key);
%     end
% end
% end
% 
% % function unwrapped = postprocess_unwrap_regions(unwrapped, quality, varargin)
% % %POSTPROCESS_UNWRAP_REGIONS Enhance multi-seed unwrapped phase continuity.
% % %   Supports stitching across disjoint regions using:
% % %       - MST-based region linking (default)
% % %       - Phase gradient blending
% % %       - Region bridging (neighborhood extrapolation)
% % %       - Small region removal or fallback
% % %
% % %   Usage:
% % %       unwrapped = postprocess_unwrap_regions(unwrapped, quality, 'method', 'mst', ...)
% % %
% % %   Options (passed as name-value pairs):
% % %       'method'     : 'mst' (default), 'blend', 'bridge', or 'none'
% % %       'tinySize'   : Minimum region size to retain (default = 30)
% % 
% % opts = struct('method', 'mst', 'tinySize', 30);
% % if nargin > 2 && ~isempty(varargin)
% %     opts = parse_inputs(opts, varargin{:});
% % end
% % 
% % % Identify connected unwrap regions
% % mask = ~isnan(unwrapped);
% % cc = bwconncomp(mask, 4);
% % 
% % % Remove or flag tiny regions
% % L = labelmatrix(cc);
% % regionSizes = cellfun(@numel, cc.PixelIdxList);
% % isTiny = regionSizes < opts.tinySize;
% % for r = find(isTiny)
% %     unwrapped(cc.PixelIdxList{r}) = NaN;
% % end
% % 
% % % Recompute region labels after pruning
% % mask = ~isnan(unwrapped);
% % [L, numRegions] = bwlabel(mask, 4);
% % 
% % switch lower(opts.method)
% %     case 'none'
% %         return
% % 
% %     case 'mst'
% %         % Compute mean phase and location per region
% %         regMeans = zeros(numRegions, 1);
% %         regCentroids = zeros(numRegions, 2);
% %         regWeights = zeros(numRegions, 1);
% %         for r = 1:numRegions
% %             idx = find(L == r);
% %             regMeans(r) = mean(unwrapped(idx), 'omitnan');
% %             regWeights(r) = mean(quality(idx), 'omitnan');
% %             [i, j] = ind2sub(size(unwrapped), idx);
% %             regCentroids(r,:) = [mean(i), mean(j)];
% %         end
% % 
% %         % Build MST weighted by centroid distance and inverse quality
% %         D = squareform(pdist(regCentroids));
% %         Q = regWeights(:);
% %         Q(Q == 0) = eps;  % avoid division by zero
% %         W = D ./ (Q + Q');
% %         W(isnan(W) | isinf(W)) = max(W(:));
% %         G = graph(W);
% %         T = minspantree(G);
% %         edges = table2array(T.Edges);
% % 
% %         for e = 1:size(edges,1)
% %             r1 = edges(e,1);
% %             r2 = edges(e,2);
% %             dp = regMeans(r2) - regMeans(r1);
% %             offset = round(dp / (2*pi)) * 2*pi;
% %             unwrapped(L == r2) = unwrapped(L == r2) - offset;
% %         end
% % 
% %     case 'blend'
% %         % Gradient blending across region boundaries
% %         se = strel('disk', 1);
% %         smoothU = unwrapped;
% %         for r = 1:numRegions
% %             maskR = (L == r);
% %             dilated = imdilate(maskR, se);
% %             boundary = dilated & ~maskR;
% %             neighborPhases = unwrapped(boundary);
% %             offset = median(neighborPhases - unwrapped(maskR), 'omitnan');
% %             smoothU(maskR) = smoothU(maskR) + offset;
% %         end
% %         unwrapped = smoothU;
% % 
% %     case 'bridge'
% %         % Simple bridging across small gaps (nearest extrapolation)
% %         se = strel('disk', 1);
% %         filled = unwrapped;
% %         for iter = 1:5
% %             nanMask = isnan(filled);
% %             grow = imdilate(~nanMask, se) & nanMask;
% %             [i, j] = find(grow);
% %             for k = 1:numel(i)
% %                 ni = i(k); nj = j(k);
% %                 local = filled(max(ni-1,1):min(ni+1,end), max(nj-1,1):min(nj+1,end));
% %                 filled(ni,nj) = mean(local(~isnan(local)), 'omitnan');
% %             end
% %         end
% %         unwrapped = filled;
% % 
% %     otherwise
% %         error('Unknown postprocess method: %s', opts.method);
% % end
% % end
% % 
% % function opts = parse_inputs(opts, varargin)
% % % Helper function for name-value pair parsing
% % if mod(numel(varargin),2) ~= 0
% %     error('Arguments must be name-value pairs.');
% % end
% % for i = 1:2:numel(varargin)
% %     key = (varargin{i});
% %     val = varargin{i+1};
% %     if isfield(opts, key)
% %         opts.(key) = val;
% %     else
% %         warning('Unknown option: %s', key);
% %     end
% % end
% % end
% % 
% % % function unwrapped = postprocess_unwrap_regions(unwrapped, quality, varargin)
% % % %POSTPROCESS_UNWRAP_REGIONS Enhance multi-seed unwrapped phase continuity.
% % % %   Supports stitching across disjoint regions using:
% % % %       - MST-based region linking (default)
% % % %       - Phase gradient blending
% % % %       - Region bridging (neighborhood extrapolation)
% % % %       - Small region removal or fallback
% % % %
% % % %   Usage:
% % % %       unwrapped = postprocess_unwrap_regions(unwrapped, quality, ...)
% % % %
% % % %   Options (passed as name-value pairs):
% % % %       'method'     : 'mst' (default), 'blend', 'bridge', or 'none'
% % % %       'tinySize'   : Minimum region size to retain (default = 30)
% % % 
% % % opts = struct('method', 'mst', 'tinySize', 30);
% % % opts = parse_inputs(opts, varargin{:});
% % % 
% % % % Identify connected unwrap regions
% % % mask = ~isnan(unwrapped);
% % % cc = bwconncomp(mask, 4);
% % % 
% % % % Remove or flag tiny regions
% % % L = labelmatrix(cc);
% % % regionSizes = cellfun(@numel, cc.PixelIdxList);
% % % isTiny = regionSizes < opts.tinySize;
% % % for r = find(isTiny)
% % %     unwrapped(cc.PixelIdxList{r}) = NaN;
% % % end
% % % 
% % % % Recompute region labels after pruning
% % % mask = ~isnan(unwrapped);
% % % [L, numRegions] = bwlabel(mask, 4);
% % % 
% % % switch lower(opts.method)
% % %     case 'none'
% % %         return
% % % 
% % %     case 'mst'
% % %         % Compute mean phase and location per region
% % %         regMeans = zeros(numRegions, 1);
% % %         regCentroids = zeros(numRegions, 2);
% % %         regWeights = zeros(numRegions, 1);
% % %         for r = 1:numRegions
% % %             idx = find(L == r);
% % %             regMeans(r) = mean(unwrapped(idx), 'omitnan');
% % %             regWeights(r) = mean(quality(idx), 'omitnan');
% % %             [i, j] = ind2sub(size(unwrapped), idx);
% % %             regCentroids(r,:) = [mean(i), mean(j)];
% % %         end
% % % 
% % %         % Build MST weighted by centroid distance and inverse quality
% % %         D = squareform(pdist(regCentroids));
% % %         Q = regWeights(:);
% % %         Q(Q == 0) = eps;  % avoid division by zero
% % %         W = D ./ (Q + Q');
% % %         W(isnan(W) | isinf(W)) = max(W(:));
% % %         G = graph(W);
% % %         T = minspantree(G);
% % %         edges = table2array(T.Edges);
% % % 
% % %         for e = 1:size(edges,1)
% % %             r1 = edges(e,1);
% % %             r2 = edges(e,2);
% % %             dp = regMeans(r2) - regMeans(r1);
% % %             offset = round(dp / (2*pi)) * 2*pi;
% % %             unwrapped(L == r2) = unwrapped(L == r2) - offset;
% % %         end
% % % 
% % %     case 'blend'
% % %         % Gradient blending across region boundaries
% % %         se = strel('disk', 1);
% % %         smoothU = unwrapped;
% % %         for r = 1:numRegions
% % %             maskR = (L == r);
% % %             dilated = imdilate(maskR, se);
% % %             boundary = dilated & ~maskR;
% % %             neighborPhases = unwrapped(boundary);
% % %             offset = median(neighborPhases - unwrapped(maskR), 'omitnan');
% % %             smoothU(maskR) = smoothU(maskR) + offset;
% % %         end
% % %         unwrapped = smoothU;
% % % 
% % %     case 'bridge'
% % %         % Simple bridging across small gaps (nearest extrapolation)
% % %         se = strel('disk', 1);
% % %         filled = unwrapped;
% % %         for iter = 1:5
% % %             nanMask = isnan(filled);
% % %             grow = imdilate(~nanMask, se) & nanMask;
% % %             [i, j] = find(grow);
% % %             for k = 1:numel(i)
% % %                 ni = i(k); nj = j(k);
% % %                 local = filled(max(ni-1,1):min(ni+1,end), max(nj-1,1):min(nj+1,end));
% % %                 filled(ni,nj) = mean(local(~isnan(local)), 'omitnan');
% % %             end
% % %         end
% % %         unwrapped = filled;
% % % 
% % %     otherwise
% % %         error('Unknown postprocess method: %s', opts.method);
% % % end
% % % end
% % % 
% % % function opts = parse_inputs(opts, varargin)
% % % % Helper function for name-value pair parsing
% % % if mod(numel(varargin),2) ~= 0
% % %     error('Arguments must be name-value pairs.');
% % % end
% % % for i = 1:2:numel(varargin)
% % %     key = lower(varargin{i});
% % %     val = varargin{i+1};
% % %     if isfield(opts, key)
% % %         opts.(key) = val;
% % %     else
% % %         warning('Unknown option: %s', key);
% % %     end
% % % end
% % % end
% % % 
% % % % function unwrapped = postprocess_unwrap_regions(unwrapped, varargin)
% % % % %POSTPROCESS_UNWRAP_REGIONS Enhance multi-seed unwrapped phase continuity.
% % % % %   Supports stitching across disjoint regions using:
% % % % %       - MST-based region linking (default)
% % % % %       - Phase gradient blending
% % % % %       - Region bridging (neighborhood extrapolation)
% % % % %       - Small region removal or fallback
% % % % %
% % % % %   Usage:
% % % % %       unwrapped = postprocess_unwrap_regions(unwrapped, ...)
% % % % %
% % % % %   Options (passed as name-value pairs):
% % % % %       'method'     : 'mst' (default), 'blend', 'bridge', or 'none'
% % % % %       'tinySize'   : Minimum region size to retain (default = 30)
% % % % 
% % % % opts = struct('method', 'mst', 'tinySize', 30);
% % % % opts = parse_inputs(opts, varargin{:});
% % % % 
% % % % % Identify connected unwrap regions
% % % % mask = ~isnan(unwrapped);
% % % % cc = bwconncomp(mask, 4);
% % % % 
% % % % % Remove or flag tiny regions
% % % % L = labelmatrix(cc);
% % % % regionSizes = cellfun(@numel, cc.PixelIdxList);
% % % % isTiny = regionSizes < opts.tinySize;
% % % % for r = find(isTiny)
% % % %     unwrapped(cc.PixelIdxList{r}) = NaN;
% % % % end
% % % % 
% % % % % Recompute region labels after pruning
% % % % mask = ~isnan(unwrapped);
% % % % [L, numRegions] = bwlabel(mask, 4);
% % % % 
% % % % switch lower(opts.method)
% % % %     case 'none'
% % % %         return
% % % % 
% % % %     case 'mst'
% % % %         % Compute mean phase and location per region
% % % %         regMeans = zeros(numRegions, 1);
% % % %         regCentroids = zeros(numRegions, 2);
% % % %         for r = 1:numRegions
% % % %             idx = find(L == r);
% % % %             regMeans(r) = mean(unwrapped(idx), 'omitnan');
% % % %             [i, j] = ind2sub(size(unwrapped), idx);
% % % %             regCentroids(r,:) = [mean(i), mean(j)];
% % % %         end
% % % % 
% % % %         % Build MST based on distance between centroids
% % % %         D = squareform(pdist(regCentroids));
% % % %         G = graph(D);
% % % %         T = minspantree(G);
% % % %         edges = table2array(T.Edges);
% % % % 
% % % %         for e = 1:size(edges,1)
% % % %             r1 = edges(e,1);
% % % %             r2 = edges(e,2);
% % % %             dp = regMeans(r2) - regMeans(r1);
% % % %             offset = round(dp / (2*pi)) * 2*pi;
% % % %             unwrapped(L == r2) = unwrapped(L == r2) - offset;
% % % %         end
% % % % 
% % % %     case 'blend'
% % % %         % Gradient blending across region boundaries
% % % %         se = strel('disk', 1);
% % % %         smoothU = unwrapped;
% % % %         for r = 1:numRegions
% % % %             maskR = (L == r);
% % % %             dilated = imdilate(maskR, se);
% % % %             boundary = dilated & ~maskR;
% % % %             neighborPhases = unwrapped(boundary);
% % % %             offset = median(neighborPhases - unwrapped(maskR), 'omitnan');
% % % %             smoothU(maskR) = smoothU(maskR) + offset;
% % % %         end
% % % %         unwrapped = smoothU;
% % % % 
% % % %     case 'bridge'
% % % %         % Simple bridging across small gaps (nearest extrapolation)
% % % %         se = strel('disk', 1);
% % % %         filled = unwrapped;
% % % %         for iter = 1:5
% % % %             nanMask = isnan(filled);
% % % %             grow = imdilate(~nanMask, se) & nanMask;
% % % %             [i, j] = find(grow);
% % % %             for k = 1:numel(i)
% % % %                 ni = i(k); nj = j(k);
% % % %                 local = filled(max(ni-1,1):min(ni+1,end), max(nj-1,1):min(nj+1,end));
% % % %                 filled(ni,nj) = mean(local(~isnan(local)), 'omitnan');
% % % %             end
% % % %         end
% % % %         unwrapped = filled;
% % % % 
% % % %     otherwise
% % % %         error('Unknown postprocess method: %s', opts.method);
% % % % end
% % % % end
% % % % 
% % % % function opts = parse_inputs(opts, varargin)
% % % % % Helper function for name-value pair parsing
% % % % if mod(numel(varargin),2) ~= 0
% % % %     error('Arguments must be name-value pairs.');
% % % % end
% % % % for i = 1:2:numel(varargin)
% % % %     key = varargin{i};
% % % %     val = varargin{i+1};
% % % %     if isfield(opts, key)
% % % %         opts.(key) = val;
% % % %     else
% % % %         warning('Unknown option: %s', key);
% % % %     end
% % % % end
% % % % end

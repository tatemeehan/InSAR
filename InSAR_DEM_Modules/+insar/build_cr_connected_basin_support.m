function [crSupportMask, supportBasinIDs, seedBasinIDs] = ...
    build_cr_connected_basin_support(owner, crResult, basinMeta, opts)
%BUILD_CR_CONNECTED_BASIN_SUPPORT
% Build a CR-supported mask from:
%   1) CR index neighborhoods
%   2) basins touched by those neighborhoods
%   3) graph expansion through contiguous basin adjacencies
%
% INPUTS
%   owner      - basin owner raster
%   crResult   - struct array with Index field
%   basinMeta  - basin metadata containing directEdges (preferred)
%   opts       - options struct
%
% OPTIONS
%   .crRadiusPx      radius around CR indices in pixels (default 30)
%   .maxGraphHops    graph expansion hops from CR basins (default 2)
%   .includeGapEdges use gap edges too (default false)
%
% OUTPUTS
%   crSupportMask  - logical raster support mask
%   supportBasinIDs - all basin IDs included in support
%   seedBasinIDs    - basin IDs directly touched by CR neighborhoods

if nargin < 4
    opts = struct();
end
if ~isfield(opts,'crRadiusPx') || isempty(opts.crRadiusPx)
    opts.crRadiusPx = 30;
end
if ~isfield(opts,'maxGraphHops') || isempty(opts.maxGraphHops)
    opts.maxGraphHops = 2;
end
if ~isfield(opts,'includeGapEdges') || isempty(opts.includeGapEdges)
    opts.includeGapEdges = false;
end

[nr, nc] = size(owner);
crSupportMask = false(nr, nc);

% ------------------------------------------------------------
% 1) Build CR neighborhood mask from CR indices
% ------------------------------------------------------------
crPts = false(nr, nc);
for ii = 1:numel(crResult)
    if ~isfield(crResult(ii),'Index') || isempty(crResult(ii).Index)
        continue;
    end

    idx = double(crResult(ii).Index);
    idx = idx(isfinite(idx));
    idx = round(idx);
    idx = idx(idx >= 1 & idx <= nr*nc);

    if ~isempty(idx)
        crPts(idx) = true;
    end
end

if ~any(crPts(:))
    supportBasinIDs = [];
    seedBasinIDs = [];
    return;
end

crNbrMask = bwdist(crPts) <= opts.crRadiusPx;

% ------------------------------------------------------------
% 2) Seed basins = basins touched by CR neighborhoods
% ------------------------------------------------------------
seedBasinIDs = unique(owner(crNbrMask & owner > 0));
seedBasinIDs = double(seedBasinIDs(:));

if isempty(seedBasinIDs)
    supportBasinIDs = [];
    crSupportMask = crNbrMask;
    return;
end

% ------------------------------------------------------------
% 3) Build basin adjacency graph
% ------------------------------------------------------------
numBasins = double(max(owner(:)));
adj = cell(numBasins,1);

if isfield(basinMeta,'directEdges') && ~isempty(basinMeta.directEdges) && ...
        isfield(basinMeta.directEdges,'a') && ~isempty(basinMeta.directEdges.a)

    a = double(basinMeta.directEdges.a(:));
    b = double(basinMeta.directEdges.b(:));

    for k = 1:numel(a)
        if a(k) >= 1 && a(k) <= numBasins && b(k) >= 1 && b(k) <= numBasins
            adj{a(k)}(end+1) = b(k); %#ok<AGROW>
            adj{b(k)}(end+1) = a(k); %#ok<AGROW>
        end
    end
end

if opts.includeGapEdges && isfield(basinMeta,'gapEdges') && ~isempty(basinMeta.gapEdges) && ...
        isfield(basinMeta.gapEdges,'a') && ~isempty(basinMeta.gapEdges.a)

    a = double(basinMeta.gapEdges.a(:));
    b = double(basinMeta.gapEdges.b(:));

    for k = 1:numel(a)
        if a(k) >= 1 && a(k) <= numBasins && b(k) >= 1 && b(k) <= numBasins
            adj{a(k)}(end+1) = b(k); %#ok<AGROW>
            adj{b(k)}(end+1) = a(k); %#ok<AGROW>
        end
    end
end

for k = 1:numBasins
    if ~isempty(adj{k})
        adj{k} = unique(adj{k});
    end
end

% ------------------------------------------------------------
% 4) BFS expansion from CR basins
% ------------------------------------------------------------
visited = false(numBasins,1);
hopDist = inf(numBasins,1);

queue = seedBasinIDs(:);
visited(queue) = true;
hopDist(queue) = 0;

head = 1;
while head <= numel(queue)
    u = queue(head);
    head = head + 1;

    if hopDist(u) >= opts.maxGraphHops
        continue;
    end

    nbrs = adj{u};
    for jj = 1:numel(nbrs)
        v = nbrs(jj);
        if ~visited(v)
            visited(v) = true;
            hopDist(v) = hopDist(u) + 1;
            queue(end+1,1) = v; %#ok<AGROW>
        end
    end
end

supportBasinIDs = find(visited);

% ------------------------------------------------------------
% 5) Rasterize support basins
% ------------------------------------------------------------
validOwner = owner > 0;
crSupportMask(validOwner) = visited(double(owner(validOwner)));

% Always include the immediate CR neighborhoods too
crSupportMask = crSupportMask | crNbrMask;

end
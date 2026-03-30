function crSupportMask = build_cr_support_mask_from_index(insarPair, crResult, basinOwner, basinMeta, opts)
%BUILD_CR_SUPPORT_MASK_FROM_INDEX
% Build a CR/reference support mask from:
%   1) anchor basin
%   2) CR pixel indices in image coordinates
%
% INPUTS
%   insarPair   - struct containing phzReferenced (defines raster size)
%   crResult    - struct array with Index field
%   basinOwner  - basin owner raster
%   basinMeta   - basin metadata with anchorBasin
%   opts        - options
%
% OPTIONS
%   .crRadiusPx        radius around CR pixels in image pixels (default 40)
%   .includeAnchorBasin include anchor basin (default true)

if nargin < 5
    opts = struct();
end
if ~isfield(opts,'crRadiusPx') || isempty(opts.crRadiusPx)
    opts.crRadiusPx = 40;
end
if ~isfield(opts,'includeAnchorBasin') || isempty(opts.includeAnchorBasin)
    opts.includeAnchorBasin = true;
end

[nr, nc] = size(insarPair.phzReferenced);
crSupportMask = false(nr, nc);

% --- 1) Anchor basin
if opts.includeAnchorBasin && ~isempty(basinOwner) && ~isempty(basinMeta) && isfield(basinMeta,'anchorBasin')
    anchorBasin = basinMeta.anchorBasin;
    if ~isempty(anchorBasin) && isfinite(anchorBasin) && anchorBasin > 0
        crSupportMask = crSupportMask | (basinOwner == anchorBasin);
    end
end

% --- 2) CR index neighborhoods
if ~isempty(crResult)
    crPts = false(nr, nc);

    for ii = 1:numel(crResult)
        if ~isfield(crResult(ii),'Index') || isempty(crResult(ii).Index)
            continue;
        end

        idx = double(crResult(ii).Index);

        % support scalar or short vector of indices
        idx = idx(isfinite(idx));
        idx = round(idx);
        idx = idx(idx >= 1 & idx <= nr*nc);

        if isempty(idx)
            continue;
        end

        crPts(idx) = true;
    end

    if any(crPts(:))
        crSupportMask = crSupportMask | (bwdist(crPts) <= opts.crRadiusPx);
    end
end

end
function unwrapped_corrected = correct_nan_isolated_regions(unwrapped, windowSize)
%CORRECT_NAN_ISOLATED_REGIONS Fix overwrapped islands surrounded by NaNs.
%
%   unwrapped_corrected = correct_nan_isolated_regions(unwrapped, wrapped, quality, windowSize)
%
%   INPUTS:
%       unwrapped   - 2D unwrapped phase map
%       wrapped     - 2D wrapped phase map (for reference)
%       quality     - 2D coherence or quality map
%       windowSize  - Size of search window around border (default = 3)
%
%   OUTPUT:
%       unwrapped_corrected - Phase map after local offset correction

if nargin < 4, windowSize = 3; end

unwrapped_corrected = unwrapped;

% Identify connected regions within overwrapped (bad) areas
overwrapMask = abs(unwrapped) > 2*pi & ~isnan(unwrapped);
cc = bwconncomp(overwrapMask, 8);

% Valid (good) phase data for reference
validMask = ~isnan(unwrapped);

for r = 1:cc.NumObjects
    regionIdx = cc.PixelIdxList{r};
    regionMask = false(size(unwrapped));
    regionMask(regionIdx) = true;

    % Dilate region to find border with valid area
    se = strel('disk', windowSize);
    dilated = imdilate(regionMask, se);
    borderMask = dilated & validMask & ~regionMask;

    if nnz(borderMask) < 5
        continue; % skip small or fully isolated regions
    end

    % Compute phase offset
    internalPhase = unwrapped(regionMask);
    borderPhase = unwrapped(borderMask);
    offset = round(median(borderPhase(:)) - median(internalPhase(:)));

    % Apply integer multiple of 2pi correction if offset large
    if abs(offset) > pi
        delta = round(offset / (2*pi)) * 2*pi;
        unwrapped_corrected(regionIdx) = unwrapped(regionIdx) - delta;
    end
end
end

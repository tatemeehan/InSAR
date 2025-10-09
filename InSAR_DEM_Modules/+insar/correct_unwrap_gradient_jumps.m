function unwrapped = correct_unwrap_gradient_jumps(unwrapped, quality, jumpThresh)
%CORRECT_UNWRAP_GRADIENT_JUMPS Detect and correct local 2pi phase jumps
%   Attempts to fix wrapped blobs in otherwise connected regions
%
%   unwrapped = correct_unwrap_gradient_jumps(unwrapped, quality, jumpThresh)
%
%   INPUTS:
%       unwrapped   - initial unwrapped phase (2D matrix)
%       quality      - quality map (same size)
%       jumpThresh   - threshold for detecting strong gradient jumps
%                      (e.g., 3 or 4 * median gradient magnitude)
%
%   OUTPUT:
%       unwrapped   - corrected phase

if nargin < 3
    jumpThresh = 10 * median(abs(gradient(unwrapped(:))), 'omitnan');
end

[rows, cols] = size(unwrapped);
[gradX, gradY] = gradient(unwrapped);
gradMag = sqrt(gradX.^2 + gradY.^2);

% Detect strong local phase jumps
% jumpMask = gradMag > jumpThresh;
jumpMask = unwrapped > 2.*pi | unwrapped < -2.*pi;
cc = bwconncomp(jumpMask, 8);

for k = 1:cc.NumObjects
    regionMask = false(rows, cols);
    regionMask(cc.PixelIdxList{k}) = true;

    % Find bordering pixels on either side of the jump
    dilated = imdilate(regionMask, strel('disk', 1));
    outerEdge = dilated & ~regionMask & ~isnan(unwrapped);
    innerRegion = regionMask & ~isnan(unwrapped);

    if nnz(outerEdge) < 5 || nnz(innerRegion) < 5
        continue; % skip small patches
    end

    outerVal = median(unwrapped(outerEdge), 'omitnan');
    innerVal = median(unwrapped(innerRegion), 'omitnan');
    dphi = innerVal - outerVal;

    % Snap to nearest multiple of 2pi
    dphi_snap = round(dphi / (2*pi)) * 2*pi;

    if abs(dphi_snap) > pi % Only apply if it's significant
        unwrapped(innerRegion) = unwrapped(innerRegion) - dphi_snap;
    end
end
end

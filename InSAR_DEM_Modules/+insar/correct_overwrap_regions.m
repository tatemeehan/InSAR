function unwrapped = correct_overwrap_regions(unwrapped, mask)
%CORRECT_OVERWRAP_REGIONS Fixes overwrapped blobs in unwrapped phase image.
%   Identifies regions exceeding +/-2pi and corrects their mean phase to
%   match surrounding pixels using boundary referencing.
%
%   INPUTS:
%       unwrapped - unwrapped phase image (2D)
%       mask      - logical mask of regions suspected to be overwrapped
%
%   OUTPUT:
%       unwrapped - corrected unwrapped phase
% function unwrapped = correct_overwrap_gradient_patch(unwrapped, quality, grad_thresh)

if nargin < 3, grad_thresh = pi; end

% Compute gradient magnitude
[dx, dy] = gradient(unwrapped);
gradMag = sqrt(dx.^2 + dy.^2);

% Threshold to detect steep gradients
mask = gradMag > grad_thresh & ~isnan(unwrapped);

% Label regions
cc = bwconncomp(mask, 8);
L = labelmatrix(cc);

for r = 1:cc.NumObjects
    regionMask = (L == r);
    
    % Get boundary of region
    se = strel('disk', 1);
    dilated = imdilate(regionMask, se);
    boundary = dilated & ~regionMask & ~isnan(unwrapped);
    
    % Skip tiny or disconnected regions
    if nnz(regionMask) < 10 || nnz(boundary) < 10, continue; end

    % Estimate offset: median boundary - median region
    offset = median(unwrapped(boundary), 'omitnan') - ...
             median(unwrapped(regionMask), 'omitnan');

    % Apply flat correction
    unwrapped(regionMask) = unwrapped(regionMask) + offset;
end
end

% % Label each disconnected blob in the overwrap mask
% cc = bwconncomp(mask, 4);
% 
% % Structuring element for boundary expansion
% se = strel('disk', 1);
% 
% for i = 1:cc.NumObjects
%     blobIdx = cc.PixelIdxList{i};
%     blobMask = false(size(unwrapped));
%     blobMask(blobIdx) = true;
% 
%     % Create a boundary mask around the blob
%     dilated = imdilate(blobMask, se);
%     boundaryMask = dilated & ~blobMask;
% 
%     % Ensure boundary is within valid phase region
%     validBoundary = boundaryMask & isfinite(unwrapped);
%     validBlob = blobMask & isfinite(unwrapped);
% 
%     if ~any(validBoundary(:)) || ~any(validBlob(:))
%         continue;  % Skip if no valid pixels
%     end
% 
%     % Compute median phase difference
%     dphi = median(unwrapped(validBlob)) - median(unwrapped(validBoundary));
% 
%     % Apply correction
%     unwrapped(blobIdx) = unwrapped(blobIdx) - dphi;
% end
% end

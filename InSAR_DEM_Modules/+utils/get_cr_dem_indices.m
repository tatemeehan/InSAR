function [crIdxs, weights, xCR, yCR, zCR] = get_cr_dem_indices(CR, crName, demData, k, method, sigma, amp)
% GET_CR_DEM_INDICES Find DEM indices and weights for a CR name.
% Includes option to weight based on max amplitude within local window.

if nargin < 4, k = 5; end
if nargin < 5, method = 'gaussian'; end
if nargin < 6, sigma = 1; end

ix = find(strcmp(CR.Name, crName));
xCR = CR.Easting(ix); 
yCR = CR.Northing(ix);
zCR = CR.Elevation(ix);

% Define larger window for amplitude-based filtering
kLarge = 2 * k + 1;
dist = sqrt((demData.X(:) - xCR).^2 + (demData.Y(:) - yCR).^2);
[~, wideIdx] = mink(dist, kLarge);

% Use amplitude to find best-matching CR pixels
if nargin >= 7 && ~isempty(amp)
    ampVals = abs(amp(wideIdx));
    [~, ampOrder] = maxk(ampVals, k);
    crIdxs = wideIdx(ampOrder);
else
    [~, crIdxs] = mink(dist, k);
end

switch lower(method)
    case 'inverse'
        epsilon = 1e-6;
        weights = 1 ./ (dist(crIdxs) + epsilon);
    case 'gaussian'
        weights = exp(-dist(crIdxs).^2 / (2 * sigma^2));
    otherwise
        error('Unknown weighting method: %s', method)
end

weights = weights / sum(weights);  % normalize
end

% function [crIdxs, weights, xCR, yCR,zCR] = get_cr_dem_indices(CR, crName, demData, k, method, sigma)
% % GET_CR_DEM_INDICES Find DEM indices and weights for a CR name
% if nargin < 4, k = 5; end
% if nargin < 5, method = 'gaussian'; end
% if nargin < 6, sigma = 1; end
% 
% ix = find(strcmp(CR.Name, crName));
% xCR = CR.Easting(ix); yCR = CR.Northing(ix);zCR = CR.EllipsoidHeight_meters_(ix);
% 
% dist = sqrt((demData.X(:) - xCR).^2 + (demData.Y(:) - yCR).^2);
% [~, crIdxs] = mink(dist, k);
% 
% switch lower(method)
%     case 'inverse'
%         epsilon = 1e-6;
%         weights = 1 ./ (dist(crIdxs) + epsilon);
%     case 'gaussian'
%         weights = exp(-dist(crIdxs).^2 / (2 * sigma^2));
%     otherwise
%         error('Unknown weighting method: %s', method)
% end
% 
% weights = weights / sum(weights);  % normalize
% end

function feat = build_phase_fill_features(demData, geomData, insarPair, opts)
%BUILD_PHASE_FILL_FEATURES Build feature rasters for phase completion.
%
% INPUTS
%   demData    - DEM struct
%   geomData   - geometry struct for the pair
%   insarPair  - current insarData(k) entry after referencing
%   opts       - options struct
%
% OUTPUT
%   feat       - struct containing feature rasters and masks
%
% Required fields expected in insarPair:
%   .phzReferenced
%   .cor
%   .power   (optional)
%
% Expected geomData fields when available:
%   .slant
%   .incidence
%   .Bperp
%
% Expected demData fields when available:
%   .X / .x
%   .Y / .y
%   .Z / .dem

if ~isfield(opts, 'trainCohThresh') || isempty(opts.trainCohThresh)
    opts.trainCohThresh = 0.75;
end
if ~isfield(opts, 'crSupportMask') || isempty(opts.crSupportMask)
    opts.crSupportMask = false(size(insarPair.phzReferenced));
end
if ~isfield(opts, 'directSolvedMask') || isempty(opts.directSolvedMask)
    opts.directSolvedMask = true(size(insarPair.phzReferenced));
end
if ~isfield(opts, 'maxTrainDistToCR') || isempty(opts.maxTrainDistToCR)
    opts.maxTrainDistToCR = inf;
end

phz = double(insarPair.phzReferenced);
cor = double(insarPair.cor);

[nr, nc] = size(phz);

% --- Coordinates
[X, Y] = get_xy_from_dem(demData, nr, nc);

% --- DEM elevation
Z = get_dem_z(demData, nr, nc);

% --- Terrain derivatives
[dZdx, dZdy] = gradient(Z);
slope = atan(sqrt(dZdx.^2 + dZdy.^2));
aspect = atan2(dZdy, dZdx);
aspect_sin = sin(aspect);
aspect_cos = cos(aspect);

[d2Zdx2, ~] = gradient(dZdx);
[~, d2Zdy2] = gradient(dZdy);
curvature = d2Zdx2 + d2Zdy2;

% --- Geometry
slant = get_field_or_nan(geomData, 'slant', nr, nc);
incidence = get_field_or_nan(geomData, 'incidence', nr, nc);
Bperp = get_field_or_nan(geomData, 'Bperp', nr, nc);

% --- Power/backscatter if available
if isfield(insarPair, 'power') && ~isempty(insarPair.power)
    power = double(insarPair.power);
else
    power = nan(nr, nc);
end

% --- Normalized image coordinates
[rr, cc] = ndgrid(1:nr, 1:nc);
rnorm = (rr - 1) ./ max(nr - 1, 1);
cnorm = (cc - 1) ./ max(nc - 1, 1);

% --- Distance to trusted phase pixels
if ~isfield(opts, 'trainCohThresh') || isempty(opts.trainCohThresh)
    opts.trainCohThresh = 0.65;
end
% 
% trustedMask = isfinite(phz) & isfinite(cor) & cor >= opts.trainCohThresh;
% distToTrusted = bwdist(trustedMask);

crSupportMask = logical(opts.crSupportMask);
directSolvedMask = logical(opts.directSolvedMask);

if ~isequal(size(crSupportMask), [nr nc])
    crSupportMask = false(nr, nc);
end
if ~isequal(size(directSolvedMask), [nr nc])
    directSolvedMask = true(nr, nc);
end

distToTrustedSeed = bwdist(isfinite(phz) & isfinite(cor) & cor >= opts.trainCohThresh);
distToCRSupport = bwdist(crSupportMask);
distToDirectSolved = bwdist(directSolvedMask & isfinite(phz));
% 
% trustedMask = ...
%     isfinite(phz) & ...
%     isfinite(cor) & ...
%     cor >= opts.trainCohThresh & ...
%     directSolvedMask & ...
%     (distToCRSupport <= opts.maxTrainDistToCR);

useCRGate = any(crSupportMask(:)) && isfinite(opts.maxTrainDistToCR);

if useCRGate
    crGate = distToCRSupport <= opts.maxTrainDistToCR;
else
    crGate = true(size(phz));
end

trustedMask = ...
    isfinite(phz) & ...
    isfinite(cor) & ...
    cor >= opts.trainCohThresh & ...
    directSolvedMask & ...
    crGate;
% --- Optional local coherence mean
localCoh = local_mean_nan(cor, 7);

% --- Valid feature mask
validFeatureMask = isfinite(Z) & isfinite(cor);

feat = struct();
feat.phz = phz;
feat.cor = cor;
feat.X = X;
feat.Y = Y;
feat.Z = Z;
feat.slope = slope;
feat.aspect_sin = aspect_sin;
feat.aspect_cos = aspect_cos;
feat.curvature = curvature;
feat.slant = slant;
feat.incidence = incidence;
feat.Bperp = Bperp;
feat.power = power;
feat.rnorm = rnorm;
feat.cnorm = cnorm;
feat.distToTrusted = distToTrustedSeed;
feat.distToCRSupport = distToCRSupport;
feat.distToDirectSolved = distToDirectSolved;
feat.localCoh = localCoh;
feat.trustedMask = trustedMask;
feat.validFeatureMask = validFeatureMask;
feat.crSupportMask = crSupportMask;
feat.directSolvedMask = directSolvedMask;
end

% -------------------------------------------------------------------------
function [X, Y] = get_xy_from_dem(demData, nr, nc)

if isfield(demData, 'X') && isfield(demData, 'Y') && ...
        isequal(size(demData.X), [nr nc]) && isequal(size(demData.Y), [nr nc])
    X = double(demData.X);
    Y = double(demData.Y);
    return;
end

if isfield(demData, 'x') && isfield(demData, 'y')
    x = double(demData.x(:)');
    y = double(demData.y(:));
    if numel(x) == nc && numel(y) == nr
        [X, Y] = meshgrid(x, y);
        return;
    end
end

[X, Y] = meshgrid(1:nc, 1:nr);
X = double(X);
Y = double(Y);

end

% -------------------------------------------------------------------------
function Z = get_dem_z(demData, nr, nc)

cand = {'Z','dem','elev','elevation'};
for i = 1:numel(cand)
    f = cand{i};
    if isfield(demData, f) && isequal(size(demData.(f)), [nr nc])
        Z = double(demData.(f));
        return;
    end
end

Z = nan(nr, nc);

end

% -------------------------------------------------------------------------
function A = get_field_or_nan(S, fieldName, nr, nc)

if isfield(S, fieldName) && ~isempty(S.(fieldName)) && isequal(size(S.(fieldName)), [nr nc])
    A = double(S.(fieldName));
else
    A = nan(nr, nc);
end

end

% -------------------------------------------------------------------------
function M = local_mean_nan(A, win)
if nargin < 2, win = 7; end
W = ones(win, win);
valid = isfinite(A);
A0 = A;
A0(~valid) = 0;
num = conv2(A0, W, 'same');
den = conv2(double(valid), W, 'same');
M = num ./ max(den, 1);
M(den == 0) = nan;
end
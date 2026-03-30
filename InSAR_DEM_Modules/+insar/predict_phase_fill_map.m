function out = predict_phase_fill_map(feat, model, opts)
%PREDICT_PHASE_FILL_MAP Predict dense referenced phase and blend with data.

if nargin < 3
    opts = struct();
end

if ~isfield(opts, 'predCohThresh') || isempty(opts.predCohThresh)
    opts.predCohThresh = 0.15;
end
if ~isfield(opts, 'maxFillDistance') || isempty(opts.maxFillDistance)
    opts.maxFillDistance = 60;
end
if ~isfield(opts, 'smoothSigma') || isempty(opts.smoothSigma)
    opts.smoothSigma = 2;
end
if ~isfield(opts, 'preserveTrusted') || isempty(opts.preserveTrusted)
    opts.preserveTrusted = true;
end
if ~isfield(opts, 'maxPredictDistToCR') || isempty(opts.maxPredictDistToCR)
    opts.maxPredictDistToCR = 250;
end
if ~isfield(opts, 'maxPredictDistToDirectSolved') || isempty(opts.maxPredictDistToDirectSolved)
    opts.maxPredictDistToDirectSolved = 150;
end

[Xmat, ~] = make_feature_matrix_from_keepmask(feat, model.featureKeepMask);

% predictMask = feat.validFeatureMask(:) & all(isfinite(Xmat), 2) & ...
%               isfinite(feat.cor(:)) & feat.cor(:) >= opts.predCohThresh & ...
%               feat.distToTrusted(:) <= opts.maxFillDistance;
% 
% predictMask = feat.validFeatureMask(:) & ...
%               all(isfinite(Xmat), 2) & ...
%               isfinite(feat.cor(:)) & ...
%               feat.cor(:) >= opts.predCohThresh & ...
%               feat.distToTrusted(:) <= opts.maxFillDistance & ...
%               feat.distToCRSupport(:) <= opts.maxPredictDistToCR & ...
%               feat.distToDirectSolved(:) <= opts.maxPredictDistToDirectSolved;
useCRGate = any(feat.crSupportMask(:)) && isfinite(opts.maxPredictDistToCR);

if useCRGate
    crGate = feat.distToCRSupport(:) <= opts.maxPredictDistToCR;
else
    crGate = true(size(feat.phz(:)));
end

predictMask = feat.validFeatureMask(:) & ...
              all(isfinite(Xmat), 2) & ...
              isfinite(feat.cor(:)) & ...
              feat.cor(:) >= opts.predCohThresh & ...
              feat.distToTrusted(:) <= opts.maxFillDistance & ...
              crGate & ...
              feat.distToDirectSolved(:) <= opts.maxPredictDistToDirectSolved;

yhat = nan(size(feat.phz(:)));
if any(predictMask)
    yhat(predictMask) = predict(model.mdl, Xmat(predictMask, :));
end

predMap = reshape(yhat, size(feat.phz));

predMapSm = predMap;
if opts.smoothSigma > 0
    predMapSm = masked_gaussian_smooth(predMap, isfinite(predMap), opts.smoothSigma);
end

filled = predMapSm;

if opts.preserveTrusted
    trusted = feat.trustedMask;
    filled(trusted) = feat.phz(trusted);
end

% support = reshape(predictMask, size(feat.phz));
% confidence = exp(-feat.distToTrusted / max(opts.maxFillDistance,1));
% confidence(~support) = 0;
% confidence(feat.trustedMask) = 1;
support = reshape(predictMask, size(feat.phz));

c1 = exp(-feat.distToTrusted / max(opts.maxFillDistance,1));
c2 = exp(-feat.distToCRSupport / max(opts.maxPredictDistToCR,1));
c3 = exp(-feat.distToDirectSolved / max(opts.maxPredictDistToDirectSolved,1));
c4 = min(max(feat.cor, 0), 1);

confidence = c1 .* c2 .* c3 .* c4;
confidence(~support) = 0;
confidence(feat.trustedMask) = 1;

out = struct();
out.predictedResidual = predMapSm;
out.filledResidual = filled;
out.supportMask = support;
out.confidence = confidence;
out.trustedMask = feat.trustedMask;

end

% -------------------------------------------------------------------------
function [Xmat, names] = make_feature_matrix_from_keepmask(feat, keepMask)

% rawNames = { ...
%     'X', 'Y', 'Z', 'slope', 'aspect_sin', 'aspect_cos', 'curvature', ...
%     'slant', 'incidence', 'Bperp', 'power', 'cor', ...
%     'rnorm', 'cnorm', 'distToTrusted', 'localCoh'};
% 
% rawCols = { ...
%     feat.X(:), feat.Y(:), feat.Z(:), feat.slope(:), ...
%     feat.aspect_sin(:), feat.aspect_cos(:), feat.curvature(:), ...
%     feat.slant(:), feat.incidence(:), feat.Bperp(:), feat.power(:), ...
%     feat.cor(:), feat.rnorm(:), feat.cnorm(:), ...
%     feat.distToTrusted(:), feat.localCoh(:)};
rawNames = { ...
    'X', 'Y', 'Z', 'slope', 'aspect_sin', 'aspect_cos', 'curvature', ...
    'slant', 'incidence', 'Bperp', 'power', 'cor', ...
    'rnorm', 'cnorm', 'distToTrusted', 'distToCRSupport', ...
    'distToDirectSolved', 'localCoh'};

rawCols = { ...
    feat.X(:), feat.Y(:), feat.Z(:), feat.slope(:), ...
    feat.aspect_sin(:), feat.aspect_cos(:), feat.curvature(:), ...
    feat.slant(:), feat.incidence(:), feat.Bperp(:), feat.power(:), ...
    feat.cor(:), feat.rnorm(:), feat.cnorm(:), feat.distToTrusted(:), ...
    feat.distToCRSupport(:), feat.distToDirectSolved(:), feat.localCoh(:)};

names = rawNames(keepMask);
cols = rawCols(keepMask);

Xmat = [];
for k = 1:numel(cols)
    Xmat = [Xmat, cols{k}]; %#ok<AGROW>
end

end

function B = masked_gaussian_smooth(A, mask, sigma)

if sigma <= 0
    B = A;
    return;
end

sz = max(3, ceil(6*sigma));
if mod(sz,2) == 0
    sz = sz + 1;
end

g = fspecial('gaussian', [sz sz], sigma);

A0 = A;
A0(~mask) = 0;

num = imfilter(A0, g, 'replicate', 'same');
den = imfilter(double(mask), g, 'replicate', 'same');

B = num ./ max(den, eps);
B(den < 1e-6) = nan;

end
% function out = predict_phase_fill_map(feat, model, opts)
% %PREDICT_PHASE_FILL_MAP Predict dense referenced phase and blend with data.
% 
% if nargin < 3
%     opts = struct();
% end
% 
% if ~isfield(opts, 'predCohThresh') || isempty(opts.predCohThresh)
%     opts.predCohThresh = 0.15;
% end
% if ~isfield(opts, 'maxFillDistance') || isempty(opts.maxFillDistance)
%     opts.maxFillDistance = 60;
% end
% if ~isfield(opts, 'smoothSigma') || isempty(opts.smoothSigma)
%     opts.smoothSigma = 2;
% end
% if ~isfield(opts, 'preserveTrusted') || isempty(opts.preserveTrusted)
%     opts.preserveTrusted = true;
% end
% 
% [Xmat, ~] = make_feature_matrix(feat);
% 
% predictMask = feat.validFeatureMask(:) & all(isfinite(Xmat), 2) & ...
%               isfinite(feat.cor(:)) & feat.cor(:) >= opts.predCohThresh & ...
%               feat.distToTrusted(:) <= opts.maxFillDistance;
% 
% yhat = nan(size(feat.phz(:)));
% if any(predictMask)
%     yhat(predictMask) = predict(model.mdl, Xmat(predictMask, :));
% end
% 
% predMap = reshape(yhat, size(feat.phz));
% 
% % Optional smoothing of predictions only
% predMapSm = predMap;
% if opts.smoothSigma > 0
%     predMapSm = masked_gaussian_smooth(predMap, isfinite(predMap), opts.smoothSigma);
% end
% 
% filled = predMapSm;
% 
% if opts.preserveTrusted
%     trusted = feat.trustedMask;
%     filled(trusted) = feat.phz(trusted);
% end
% 
% % confidence/support
% support = reshape(predictMask, size(feat.phz));
% confidence = exp(-feat.distToTrusted / max(opts.maxFillDistance,1));
% confidence(~support) = 0;
% confidence(feat.trustedMask) = 1;
% 
% out = struct();
% out.predictedResidual = predMapSm;
% out.filledResidual = filled;
% out.supportMask = support;
% out.confidence = confidence;
% out.trustedMask = feat.trustedMask;
% 
% end
% 
% % -------------------------------------------------------------------------
% function [Xmat, names] = make_feature_matrix(feat)
% 
% names = { ...
%     'X', 'Y', 'Z', 'slope', 'aspect_sin', 'aspect_cos', 'curvature', ...
%     'slant', 'incidence', 'Bperp', 'power', 'cor', ...
%     'rnorm', 'cnorm', 'distToTrusted', 'localCoh'};
% 
% Xmat = [ ...
%     feat.X(:), feat.Y(:), feat.Z(:), feat.slope(:), ...
%     feat.aspect_sin(:), feat.aspect_cos(:), feat.curvature(:), ...
%     feat.slant(:), feat.incidence(:), feat.Bperp(:), feat.power(:), ...
%     feat.cor(:), feat.rnorm(:), feat.cnorm(:), ...
%     feat.distToTrusted(:), feat.localCoh(:)];
% end
% 
% % -------------------------------------------------------------------------
% function B = masked_gaussian_smooth(A, mask, sigma)
% 
% if sigma <= 0
%     B = A;
%     return;
% end
% 
% sz = max(3, ceil(6*sigma));
% if mod(sz,2) == 0
%     sz = sz + 1;
% end
% 
% g = fspecial('gaussian', [sz sz], sigma);
% 
% A0 = A;
% A0(~mask) = 0;
% 
% num = imfilter(A0, g, 'replicate', 'same');
% den = imfilter(double(mask), g, 'replicate', 'same');
% 
% B = num ./ max(den, eps);
% B(den < 1e-6) = nan;
% 
% end
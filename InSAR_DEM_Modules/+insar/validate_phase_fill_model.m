function val = validate_phase_fill_model(feat, fillOpts, predOpts, mode, varargin)
%VALIDATE_PHASE_FILL_MODEL Validate phase fill by holdout reconstruction.
%
% INPUTS
%   feat      - feature struct from build_phase_fill_features
%   fillOpts  - training options
%   predOpts  - prediction options
%   mode      - 'random' or 'block'
%
% OPTIONAL
%   For mode='random':
%       holdoutFrac  - fraction of trusted pixels to hold out (default 0.15)
%   For mode='block':
%       blockRadius  - patch radius in pixels (default 20)
%       numBlocks    - number of blocks (default 25)
%
% OUTPUT
%   val       - struct with metrics and masks

if nargin < 4
    error('Need feat, fillOpts, predOpts, and mode.');
end

trustedMask0 = logical(feat.trustedMask);
if nnz(trustedMask0) < 100
    error('Trusted mask too small for validation.');
end

switch lower(mode)
    case 'random'
        holdoutFrac = 0.15;
        if ~isempty(varargin), holdoutFrac = varargin{1}; end
        holdoutMask = insar.make_random_holdout_mask(trustedMask0, holdoutFrac);

    case 'block'
        blockRadius = 20;
        numBlocks = 25;
        if numel(varargin) >= 1, blockRadius = varargin{1}; end
        if numel(varargin) >= 2, numBlocks = varargin{2}; end
        holdoutMask = insar.make_block_holdout_mask(trustedMask0, blockRadius, numBlocks);

    otherwise
        error('Unknown mode: %s', mode);
end

trainMask = trustedMask0 & ~holdoutMask;

% % Build a training feature struct identical to feat, but with modified trustedMask
% featTrain = feat;
% featTrain.trustedMask = trainMask;
% 
% % Train
% model = insar.train_phase_fill_model(featTrain, fillOpts);
% 
% % Predict
% pred = insar.predict_phase_fill_map(featTrain, model, predOpts);

% Build a training feature struct identical to feat, but with modified trustedMask
featTrain = feat;
featTrain.trustedMask = trainMask;

% Recompute any mask-derived predictors so holdout pixels are not leaked back in
featTrain.distToTrusted = bwdist(featTrain.trustedMask);

% If directSolvedMask is intended to represent available support for validation,
% you can optionally intersect it with training support as well:
% featTrain.directSolvedMask = feat.directSolvedMask & trainMask;
% featTrain.distToDirectSolved = bwdist(featTrain.directSolvedMask & isfinite(featTrain.phz));

% Train
model = insar.train_phase_fill_model(featTrain, fillOpts);

% Predict using the training-support definition, not the original full trusted mask
pred = insar.predict_phase_fill_map(featTrain, model, predOpts);

% Observed target on held-out trusted pixels
yTrue = feat.phz(holdoutMask);
yPred = pred.predictedResidual(holdoutMask);
conf  = pred.confidence(holdoutMask);

% Raw metrics on all holdout pixels with finite prediction
valid = isfinite(yTrue) & isfinite(yPred);

metricsRaw = compute_metrics(yTrue(valid), yPred(valid));

% Confidence-threshold sweep
% thresholds = [0.3 0.5 0.7];

thresholds = [0.5 0.75 0.9];
metricsThresh = repmat(struct( ...
    'threshold', nan, ...
    'n', 0, ...
    'fracRetained', nan, ...
    'rmse', nan, ...
    'mae', nan, ...
    'bias', nan, ...
    'corr', nan), numel(thresholds), 1);

for ii = 1:numel(thresholds)
    t = thresholds(ii);
    keep = valid & (conf >= t);

    metricsThresh(ii).threshold = t;
    metricsThresh(ii).n = nnz(keep);
    metricsThresh(ii).fracRetained = nnz(keep) / max(nnz(valid), 1);

    if any(keep)
        mm = compute_metrics(yTrue(keep), yPred(keep));
        metricsThresh(ii).rmse = mm.rmse;
        metricsThresh(ii).mae  = mm.mae;
        metricsThresh(ii).bias = mm.bias;
        metricsThresh(ii).corr = mm.corr;
    end
end

val = struct();
val.mode = mode;
val.model = model;
val.pred = pred;
val.trainMask = trainMask;
val.holdoutMask = holdoutMask;
val.metricsRaw = metricsRaw;
val.metricsThresh = metricsThresh;
val.nTrusted = nnz(trustedMask0);
val.nTrain = nnz(trainMask);
val.nHoldout = nnz(holdoutMask);

end

% -------------------------------------------------------------------------
function m = compute_metrics(yTrue, yPred)

err = yPred - yTrue;

m = struct();
m.n = numel(err);
m.rmse = sqrt(mean(err.^2, 'omitnan'));
m.mae  = mean(abs(err), 'omitnan');
m.bias = mean(err, 'omitnan');

if numel(yTrue) > 1
    C = corrcoef(yTrue, yPred, 'Rows', 'complete');
    if numel(C) >= 4
        m.corr = C(1,2);
    else
        m.corr = nan;
    end
else
    m.corr = nan;
end

end
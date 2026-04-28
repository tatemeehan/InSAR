function out = inpaint_interferogram_ml_committee(insarPair, G, basinOwner, basinMeta, crResult, demData, opts)
%INPAINT_INTERFEROGRAM_ML
% Production phase inpainting using shared CR-anchored support.
%
% Supports:
%   opts.method = 'ml'       -> bagged trees / selected ML model
%   opts.method = 'gaussian' -> Gaussian phase interpolation from trusted mask
%
% Optional committee mode for ML:
%   opts.useCommittee = true/false
%
% Outputs:
%   phzFilled, phzWrappedFilled, phaseQuality, trusted/support masks, etc.

phz = double(insarPair.phzReferenced);
cor = double(insarPair.cor);

% -------------------------------------------------------------------------
% CR-connected support mask
% -------------------------------------------------------------------------
supportOpts = struct();
supportOpts.crRadiusPx = opts.crRadiusPx;
supportOpts.maxGraphHops = opts.maxGraphHops;
supportOpts.includeGapEdges = opts.includeGapEdges;

[crSupportMask, supportBasinIDs, seedBasinIDs] = insar.build_cr_connected_basin_support( ...
    basinOwner, crResult, basinMeta, supportOpts);

% -------------------------------------------------------------------------
% Direct solved mask
% -------------------------------------------------------------------------
directSolvedMask = true(size(phz));
if ~isempty(basinOwner) && isfield(basinMeta, 'solved') && ~isempty(basinMeta.solved)
    solved = basinMeta.solved(:);
    directSolvedMask = false(size(basinOwner));
    validOwner = basinOwner > 0 & double(basinOwner) <= numel(solved);
    directSolvedMask(validOwner) = solved(double(basinOwner(validOwner)));
end

% -------------------------------------------------------------------------
% Same-wrap mask relative to CR-supported phase
% -------------------------------------------------------------------------
sameWrapMask = true(size(phz));
if any(crSupportMask(:))
    z = exp(1i * phz(crSupportMask));
    phiCR = angle(mean(z, 'omitnan'));
    dphiCR = mod(phz - phiCR + pi, 2*pi) - pi;
    sameWrapMask = abs(dphiCR) <= opts.sameWrapThresh;
else
    phiCR = nan;
    dphiCR = nan(size(phz));
end

% -------------------------------------------------------------------------
% Shared support/trusted mask
% -------------------------------------------------------------------------
maskOpts = struct();
maskOpts.trainCohThresh = opts.trainCohThresh;
maskOpts.crSupportMask = crSupportMask;
maskOpts.directSolvedMask = directSolvedMask;
maskOpts.maxTrainDistToCR = opts.maxTrainDistToCR;
maskOpts.sameWrapMask = sameWrapMask;

support = insar.create_trusted_mask(phz, cor, maskOpts);
trustedMask = support.trustedMask;

% -------------------------------------------------------------------------
% Fill phase
% -------------------------------------------------------------------------
phaseFillSpread = [];

switch lower(opts.method)
    case 'gaussian'
        phzFilled = insar.phase_fill_gaussian_interp(phz, trustedMask, opts.phaseSigma);
        phzFilled(trustedMask) = phz(trustedMask);
        modelInfo = struct('method','gaussian');

    otherwise
        % ML feature build
        featOpts = maskOpts;
        featOpts.minValidFrac = opts.minValidFrac;

        feat = insar.build_phase_fill_features(demData, G, insarPair, featOpts);
        feat.trustedMask = trustedMask;  % enforce shared support

        % training options
        fillOpts = struct();
        fillOpts.trainCohThresh = opts.trainCohThresh;
        fillOpts.maxTrainSamples = opts.maxTrainSamples;
        fillOpts.modelType = opts.modelType;
        fillOpts.numTrees = opts.numTrees;
        fillOpts.minLeafSize = opts.minLeafSize;
        fillOpts.minValidFrac = opts.minValidFrac;
        fillOpts.learnRate = opts.learnRate;
        fillOpts.selectedFeatures = opts.selectedFeatures;

        % committee options
        if isfield(opts,'useCommittee'),         fillOpts.useCommittee = opts.useCommittee; end
        if isfield(opts,'nCommittee'),           fillOpts.nCommittee = opts.nCommittee; end
        if isfield(opts,'treesPerModel'),        fillOpts.treesPerModel = opts.treesPerModel; end
        if isfield(opts,'subsampleFrac'),        fillOpts.subsampleFrac = opts.subsampleFrac; end
        if isfield(opts,'baseSeed'),             fillOpts.baseSeed = opts.baseSeed; end
        if isfield(opts,'committeeAggregate'),   fillOpts.committeeAggregate = opts.committeeAggregate; end
        if isfield(opts,'storeCommitteeSpread'), fillOpts.storeCommitteeSpread = opts.storeCommitteeSpread; end

        % prediction options
        predOpts = struct();
        predOpts.predCohThresh = opts.predCohThresh;
        predOpts.maxFillDistance = opts.maxFillDistance;
        predOpts.maxPredictDistToCR = opts.maxPredictDistToCR;
        predOpts.maxPredictDistToDirectSolved = opts.maxPredictDistToDirectSolved;
        predOpts.smoothSigma = opts.smoothSigma;
        predOpts.preserveTrusted = opts.preserveTrusted;

        % -----------------------------
        % Committee or single model
        % -----------------------------
        useCommittee = isfield(fillOpts,'useCommittee') && fillOpts.useCommittee;

        if useCommittee
            committee = train_phase_fill_committee(feat, fillOpts);
            committeePred = predict_phase_fill_committee(feat, committee, predOpts);

            aggregateMode = 'median';
            if isfield(fillOpts,'committeeAggregate') && ~isempty(fillOpts.committeeAggregate)
                aggregateMode = lower(fillOpts.committeeAggregate);
            end

            switch aggregateMode
                case 'mean'
                    phzFilled = committeePred.filledResidualMean;
                otherwise
                    phzFilled = committeePred.filledResidualMedian;
            end

            storeSpread = false;
            if isfield(fillOpts,'storeCommitteeSpread') && fillOpts.storeCommitteeSpread
                storeSpread = true;
            end
            if storeSpread
                phaseFillSpread = committeePred.filledResidualMAD;
            end

            % preserve trusted pixels exactly
            phzFilled(trustedMask) = phz(trustedMask);

            modelInfo = struct();
            modelInfo.method = 'ml_committee';
            modelInfo.nCommittee = committee.opts.nCommittee;
            modelInfo.treesPerModel = committee.opts.treesPerModel;
            modelInfo.subsampleFrac = committee.opts.subsampleFrac;
            modelInfo.aggregate = aggregateMode;
            modelInfo.featureNames = committee.featureInfo.featureNames;

        else
            model = insar.train_phase_fill_model(feat, fillOpts);
            pred = insar.predict_phase_fill_map(feat, model, predOpts);
            phzFilled = pred.filledResidual;

            % preserve trusted pixels exactly
            phzFilled(trustedMask) = phz(trustedMask);

            modelInfo = rmfield(model, 'mdl');
            modelInfo.method = 'ml_single';
        end
end

% -------------------------------------------------------------------------
% Wrapped version of the filled phase
% -------------------------------------------------------------------------
phzWrappedFilled = angle(exp(1i * phzFilled));

% -------------------------------------------------------------------------
% Phase-quality map = synthetic coherence metric for the filled phase
% NOTE: using your updated signature:
%   phase_quality_gaussian(wrapped_phase, sigma, minValidFrac, validMask)
% -------------------------------------------------------------------------
validMaskForQuality = support.validMask & isfinite(phzFilled);
[phaseQuality, qualitySupportFrac] = insar.phase_quality_gaussian( ...
    phzWrappedFilled, opts.qualitySigma, opts.qualityMinValidFrac, validMaskForQuality);

% -------------------------------------------------------------------------
% Optional quality-thresholded masked product
% -------------------------------------------------------------------------
acceptMask = [];
phzFilledMasked = [];
if ~isempty(opts.qualityThreshold)
    acceptMask = phaseQuality >= opts.qualityThreshold;
    phzFilledMasked = phzFilled;
    phzFilledMasked(~acceptMask) = NaN;
end

% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
out = struct();
out.phzFilled = phzFilled;
out.phzWrappedFilled = phzWrappedFilled;
out.phaseQuality = phaseQuality;
out.qualitySupportFrac = qualitySupportFrac;

out.trustedMask = trustedMask;
out.crSupportMask = crSupportMask;
out.directSolvedMask = directSolvedMask;
out.sameWrapMask = sameWrapMask;
out.support = support;
out.supportBasinIDs = supportBasinIDs;
out.seedBasinIDs = seedBasinIDs;
out.phiCR = phiCR;
out.dphiCR = dphiCR;

out.acceptMask = acceptMask;
out.phzFilledMasked = phzFilledMasked;
out.modelInfo = modelInfo;
out.phaseFillSpread = phaseFillSpread;

end

%% Helper Functions
%--------------------------------------------------------------------------
function committee = train_phase_fill_committee(feat, opts)
%TRAIN_PHASE_FILL_COMMITTEE Train a committee of bagged ensembles.

if ~isfield(opts,'nCommittee') || isempty(opts.nCommittee)
    opts.nCommittee = 5;
end
if ~isfield(opts,'treesPerModel') || isempty(opts.treesPerModel)
    opts.treesPerModel = 75;
end
if ~isfield(opts,'subsampleFrac') || isempty(opts.subsampleFrac)
    opts.subsampleFrac = 0.7;
end
if ~isfield(opts,'baseSeed') || isempty(opts.baseSeed)
    opts.baseSeed = 1;
end
if ~isfield(opts,'committeeAggregate') || isempty(opts.committeeAggregate)
    opts.committeeAggregate = 'median';
end

committee = struct();
committee.models = cell(opts.nCommittee,1);
committee.opts = opts;

trainOptsBase = opts;
trainOptsBase.modelType = 'bag';
trainOptsBase.numTrees = opts.treesPerModel;

[Xall, yall, trainMask, featureInfo] = make_training_matrix_for_committee(feat, trainOptsBase);

nTrain = size(Xall,1);
if nTrain < 2
    error('train_phase_fill_committee:NotEnoughTrainingData', ...
        'Not enough finite trusted training rows for committee training.');
end

nSub = round(opts.subsampleFrac * nTrain);
nSub = max(10, nSub);
nSub = min(nTrain, nSub);

for k = 1:opts.nCommittee
    rng(opts.baseSeed + k - 1, 'twister');

    idx = randperm(nTrain, nSub);
    Xk = Xall(idx,:);
    yk = yall(idx);

    mdl = fitrensemble( ...
        Xk, yk, ...
        'Method', 'Bag', ...
        'NumLearningCycles', trainOptsBase.numTrees, ...
        'Learners', templateTree('MinLeafSize', trainOptsBase.minLeafSize));

    modelk = struct();
    modelk.mdl = mdl;
    modelk.featureKeepMask = featureInfo.featureKeepMask;
    modelk.featureNames = featureInfo.featureNames;
    modelk.opts = trainOptsBase;
    modelk.subsampleIdx = idx;

    committee.models{k} = modelk;
end

committee.featureInfo = featureInfo;
committee.trainMask = trainMask;
end

%--------------------------------------------------------------------------
function out = predict_phase_fill_committee(feat, committee, predOpts)
%PREDICT_PHASE_FILL_COMMITTEE Predict with multiple bagged ensembles and aggregate.

nCommittee = numel(committee.models);
predStack = [];

for k = 1:nCommittee
    predk = insar.predict_phase_fill_map(feat, committee.models{k}, predOpts);

    if k == 1
        predStack = nan([size(predk.filledResidual), nCommittee]);
    end

    predStack(:,:,k) = predk.filledResidual;
end

out = struct();
out.predStack = predStack;
out.filledResidualMedian = median(predStack, 3, 'omitnan');
out.filledResidualMean   = mean(predStack, 3, 'omitnan');
out.filledResidualMAD    = median(abs(predStack - out.filledResidualMedian), 3, 'omitnan');
out.filledResidualSTD    = std(predStack, 0, 3, 'omitnan');
end

%--------------------------------------------------------------------------
function [Xall, yall, trainMask, featureInfo] = make_training_matrix_for_committee(feat, opts)
%MAKE_TRAINING_MATRIX_FOR_COMMITTEE
% Mirror the feature-selection logic from train_phase_fill_model.

if ~isfield(opts,'selectedFeatures') || isempty(opts.selectedFeatures)
    opts.selectedFeatures = {'Z','slope','aspect_sin','aspect_cos','power','cor','incidence'};
end

allNames = {'Z','V','slope','aspect_sin','aspect_cos','curvature','power','cor','localCoh','incidence','Bperp','slant'};
allCols = { ...
    feat.Z(:), feat.V(:), feat.slope(:), feat.aspect_sin(:), feat.aspect_cos(:), feat.curvature(:), ...
    feat.power(:), feat.cor(:), feat.localCoh(:), feat.incidence(:), feat.Bperp(:), feat.slant(:)};

keep = ismember(allNames, opts.selectedFeatures);
rawNames = allNames(keep);
rawCols = allCols(keep);

Xraw = [];
for i = 1:numel(rawCols)
    Xraw = [Xraw, rawCols{i}]; %#ok<AGROW>
end

y = feat.phz(:);

% OG behavior: train from trustedMask unless a trainMask exists
if isfield(feat, 'trainMask') && ~isempty(feat.trainMask)
    trainMask = feat.trainMask(:) & feat.validFeatureMask(:) & isfinite(y);
else
    trainMask = feat.trustedMask(:) & feat.validFeatureMask(:) & isfinite(y);
end

% Match single-model fallback behavior
if nnz(trainMask) < 1000
    warning('Committee training mask too small; falling back to coherence-only training mask.');
    trainMask = isfinite(y) & feat.validFeatureMask(:) & isfinite(feat.cor(:)) & ...
                feat.cor(:) >= opts.trainCohThresh;
end

if ~any(trainMask)
    error('make_training_matrix_for_committee:EmptyTrainMask', ...
        'No usable training rows remain after fallback.');
end

validFrac = mean(isfinite(Xraw(trainMask,:)), 1);
featureKeepMask = validFrac >= opts.minValidFrac;

Xraw = Xraw(:, featureKeepMask);
featureNames = rawNames(featureKeepMask);

rowMask = trainMask & all(isfinite(Xraw), 2);

Xall = Xraw(rowMask,:);
yall = y(rowMask);

featureInfo = struct();
featureInfo.featureKeepMask = featureKeepMask;
featureInfo.featureNames = featureNames;
featureInfo.rowMask = rowMask;
end
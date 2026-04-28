function out = inpaint_interferogram_ml_og(insarPair, G, basinOwner, basinMeta, crResult, demData, opts)
%INPAINT_INTERFEROGRAM_ML
% Production phase inpainting using shared CR-anchored support.
%
% Supports:
%   opts.method = 'ml'       -> bagged trees / selected ML model
%   opts.method = 'gaussian' -> Gaussian phase interpolation from trusted mask
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

        % optionally prune features centrally
        fillOpts = struct();
        fillOpts.trainCohThresh = opts.trainCohThresh;
        fillOpts.maxTrainSamples = opts.maxTrainSamples;
        fillOpts.modelType = opts.modelType;
        fillOpts.numTrees = opts.numTrees;
        fillOpts.minLeafSize = opts.minLeafSize;
        fillOpts.minValidFrac = opts.minValidFrac;
        fillOpts.learnRate = opts.learnRate;
        fillOpts.selectedFeatures = opts.selectedFeatures;

        predOpts = struct();
        predOpts.predCohThresh = opts.predCohThresh;
        predOpts.maxFillDistance = opts.maxFillDistance;
        predOpts.maxPredictDistToCR = opts.maxPredictDistToCR;
        predOpts.maxPredictDistToDirectSolved = opts.maxPredictDistToDirectSolved;
        predOpts.smoothSigma = opts.smoothSigma;
        predOpts.preserveTrusted = opts.preserveTrusted;

        model = insar.train_phase_fill_model(feat, fillOpts);
        pred = insar.predict_phase_fill_map(feat, model, predOpts);
        phzFilled = pred.filledResidual;

        modelInfo = rmfield(model, 'mdl');
        modelInfo.method = 'ml';
    end

% -------------------------------------------------------------------------
% Wrapped version of the filled phase
% -------------------------------------------------------------------------
phzWrappedFilled = angle(exp(1i * phzFilled));

% -------------------------------------------------------------------------
% Phase-quality map = synthetic coherence metric for the filled phase
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

end
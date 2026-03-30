% fillOpts = struct();
% fillOpts.trainCohThresh = 0.65;
% fillOpts.maxTrainSamples = 150000;
% fillOpts.modelType = 'bag';      % or 'lsboost'
% fillOpts.numTrees = 100;
% fillOpts.minLeafSize = 10;
% 
% predOpts = struct();
% predOpts.predCohThresh = 0.10;
% predOpts.maxFillDistance = 75;
% predOpts.smoothSigma = 2;
% predOpts.preserveTrusted = true;
% 
% insarPair = struct();
% insarPair.phzReferenced = insarData(k).phzReferenced;
% insarPair.cor = insarData(k).coherence;
% 
% if isfield(insarData(k), 'power')
%     insarPair.power = insarData(k).power;
% elseif isfield(sarData(i), 'db') && numel(sarData(i).db) >= b1
%     insarPair.power = sarData(i).db{b1};
% end

fillOpts = struct();
fillOpts.trainCohThresh = 0.85;
fillOpts.maxTrainDistToCR = inf;
fillOpts.maxTrainSamples = 100000;
fillOpts.modelType = 'bag';
fillOpts.numTrees = 250;
fillOpts.minLeafSize = 5;
fillOpts.minValidFrac = 0.25;

% Support-aware additions
% fillOpts.maxTrainDistToCR = 150;

predOpts = struct();
predOpts.predCohThresh = 0.15;
predOpts.maxFillDistance = 50;
predOpts.maxPredictDistToCR = inf;
predOpts.maxPredictDistToDirectSolved = 75;
predOpts.smoothSigma = 3;
predOpts.preserveTrusted = true;

% % crSupportMask = false(size(insarData(k).phzReferenced));
% owner = insarData(k).basinOwner;
% anchorBasin = insarData(k).basinMeta.anchorBasin;
% 
% crSupportMask = (owner == anchorBasin);
% 
% if isfield(insarData(k), 'refMeta') && isstruct(insarData(k).refMeta)
%     if isfield(insarData(k).refMeta, 'referenceMask') && ...
%             isequal(size(insarData(k).refMeta.referenceMask), size(crSupportMask))
%         crSupportMask = logical(insarData(k).refMeta.referenceMask);
%     elseif isfield(insarData(k).refMeta, 'mask') && ...
%             isequal(size(insarData(k).refMeta.mask), size(crSupportMask))
%         crSupportMask = logical(insarData(k).refMeta.mask);
%     end
% end
% 
% % Fallback: use trusted CR result masks if available
% if ~any(crSupportMask(:)) && isfield(crResult, 'referenceMask') && ...
%         isequal(size(crResult.referenceMask), size(crSupportMask))
%     crSupportMask = logical(crResult.referenceMask);
% end

crMaskOpts = struct();
crMaskOpts.crRadiusPx = 40;          % tune
crMaskOpts.includeAnchorBasin = true;

basinOwner = [];
basinMeta = [];
if isfield(insarData(k),'basinOwner')
    basinOwner = insarData(k).basinOwner;
end
if isfield(insarData(k),'basinMeta')
    basinMeta = insarData(k).basinMeta;
end

insarPair = struct();
insarPair.phzReferenced = insarData(k).phzReferenced;

% crSupportMask = insar.build_cr_support_mask_from_index( ...
    % insarPair, crResult, basinOwner, basinMeta, crMaskOpts);
    supportOpts = struct();
supportOpts.crRadiusPx = 100;
supportOpts.maxGraphHops = 50;
supportOpts.includeGapEdges = false;

[crSupportMask, supportBasinIDs, seedBasinIDs] = ...
    insar.build_cr_connected_basin_support( ...
        insarData(k).basinOwner, crResult, insarData(k).basinMeta, supportOpts);

% fillOpts.crSupportMask = crSupportMask;

fillOpts.crSupportMask = crSupportMask;

directSolvedMask = true(size(insarData(k).phzReferenced));
if isfield(insarData(k), 'basinOwner') && isfield(insarData(k), 'basinMeta') && ...
        isfield(insarData(k).basinMeta, 'solved')

    owner = insarData(k).basinOwner;
    solved = insarData(k).basinMeta.solved(:);

    directSolvedMask = false(size(owner));
    validOwner = owner > 0;
    directSolvedMask(validOwner) = solved(double(owner(validOwner)));
end

insarPair = struct();
insarPair.phzReferenced = insarData(k).phzReferenced;
insarPair.cor = insarData(k).coherence;

if isfield(insarData(k), 'power')
    insarPair.power = insarData(k).power;
elseif isfield(sarData(i), 'db') && numel(sarData(i).db) >= b1
    insarPair.power = sarData(i).db{b1};
end

fillOpts.crSupportMask = crSupportMask;
fillOpts.directSolvedMask = directSolvedMask;

feat = insar.build_phase_fill_features(demData, geomData(k), insarPair, fillOpts);
phaseFillModel = insar.train_phase_fill_model(feat, fillOpts);
phaseFillOut = insar.predict_phase_fill_map(feat, phaseFillModel, predOpts);

% feat = insar.build_phase_fill_features(demData, geomData, insarPair, fillOpts);
% fprintf('trustedMask count: %d\n', nnz(feat.trustedMask));
% fprintf('validFeatureMask count: %d\n', nnz(feat.validFeatureMask));
% fprintf('finite phz count: %d\n', nnz(isfinite(feat.phz)));
% phaseFillModel = insar.train_phase_fill_model(feat, fillOpts);
% phaseFillOut = insar.predict_phase_fill_map(feat, phaseFillModel, predOpts);

% insarData(k).phaseFillModelInfo = rmfield(phaseFillModel, 'mdl');
% insarData(k).phzReferencedFilled = phaseFillOut.filledResidual;
% insarData(k).phzReferencedPredicted = phaseFillOut.predictedResidual;
% insarData(k).phaseFillSupportMask = phaseFillOut.supportMask;
% insarData(k).phaseFillConfidence = phaseFillOut.confidence;
% insarData(k).phaseFillTrustedMask = phaseFillOut.trustedMask;

obs = feat.phz;
pred = phaseFillOut.predictedResidual;

evalMask = feat.trustedMask & isfinite(pred);

rmse = sqrt(mean((pred(evalMask) - obs(evalMask)).^2, 'omitnan'));
mae  = mean(abs(pred(evalMask) - obs(evalMask)), 'omitnan');

fprintf('Phase fill RMSE: %.3f rad\n', rmse);
fprintf('Phase fill MAE : %.3f rad\n', mae);
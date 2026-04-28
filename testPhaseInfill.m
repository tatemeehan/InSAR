%% Load Data And Such
matDir = 'E:\MCS\MCS022326\GLSAR\Export\MCS__lband-mcs_2026_TMA_2ms_200MHz_20260223T224040_TMA_TX2_RX2_V_V_V__lband-mcs_2026_TMA_2ms_200MHz_20260226T184834_TMA_TX2_RX2_V_V_V__standard__20260331T124039\mats';
insarData = io.load_export_struct(fullfile(matDir,'insarData_standard.mat'), 'insarSave');
geomData = io.load_export_struct(fullfile(matDir,'geomData_standard.mat'), 'geomSave');
sarData = io.load_export_struct(fullfile(matDir,'sarData_standard.mat'), 'sarSave');
demData = io.load_export_struct(fullfile(matDir,'demData_standard.mat'), 'demSave');
insarClosureData = io.load_export_struct(fullfile(matDir,'insarClosureData_standard.mat'), 'closureSave');

k = 1; p = 1; i = 1;
b1 = 1;

% Processing Parameters
params.c          = 0.3;      % Wave speed sonstant (m/ns)
params.f          = 1.3;      % Frequency (GHz)
params.lambda     = params.c / params.f;    % Radar wavelength (m)
params.quality    = 0.75;     % Coherence threshold for unwrapping
params.dx         = 5;        % Trajectory resampling resolution (m)
params.filterSize = 5;        % Multilook filter window size (pixels)

crDir = 'E:\MCS\MCS_CR_SM_25';
crFn = 'MCS_CRs_2026_Corrected.csv';
crPath = fullfile(crDir,crFn);
crName = {'CR1','CR2'};

params.k = 5; % Nearest Neighborscr
params.weightType = 'gaussian'; % Use Gaussian Weighting Kernel
params.sigma = 1; % meters
% Check for CR existance and if no CR make empty CR table
crResult = utils.CRpower(crPath, crName, demData, sarData, geomData, params.lambda, params.k, params.weightType, params.sigma);
%% Test ML
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
fillOpts.trainCohThresh = 0.75;
fillOpts.maxTrainDistToCR = 100;
fillOpts.maxTrainSamples = 100000;
fillOpts.modelType = 'bag';
fillOpts.numTrees = 250;
fillOpts.minLeafSize = 15;
fillOpts.minValidFrac = 0.25;

% fillOpts.modelType = 'lsboost';
% fillOpts.numTrees = 300;
% fillOpts.minLeafSize = 5;
% fillOpts.learnRate = 0.05;

% Support-aware additions
% fillOpts.maxTrainDistToCR = 150;

predOpts = struct();
predOpts.predCohThresh = 0.25;
predOpts.maxFillDistance = 500;
predOpts.maxPredictDistToCR = inf;
predOpts.maxPredictDistToDirectSolved = 150;
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

% CR Phase
z = exp(1i * insarPair.phzReferenced(crSupportMask));
phiCR = angle(mean(z, 'omitnan'));
% Wrapped Phase Difference to Reference
dphiCR = mod(insarPair.phzReferenced - phiCR + pi, 2*pi) - pi;
sameWrapMask = abs(dphiCR) <= 1;   % radians, tune
idxG = 1;
G = geomData.sarGeometry{idxG};
feat = insar.build_phase_fill_features(demData, G, insarPair, fillOpts);
[dx, dy] = gradient(feat.phz);
jumpMask = hypot(dx, dy) > .625;   % tune
jumpMask = bwdist(jumpMask) <= 10; % 3-pixel buffer

% feat.trustedMask = feat.trustedMask & ~jumpMask;
tic
feat.trustedMask = feat.trustedMask & sameWrapMask;
phaseFillModel = insar.train_phase_fill_model(feat, fillOpts);
phaseFillOut = insar.predict_phase_fill_map(feat, phaseFillModel, predOpts);
toc

z = exp(1i * phaseFillOut.filledResidual);
tmprewrapped = angle(z);
rewrapped = mod(tmprewrapped+ pi, 2*pi) - pi;
Q = insar.phase_quality_gaussian(rewrapped,11);

% acceptMask = phaseFillOut.supportMask & (phaseFillOut.confidence >= 0.25);
acceptMask = Q >= 0.7;


filledResidualMasked = phaseFillOut.filledResidual;
filledResidualMasked(~acceptMask) = NaN;
% filledResidualMasked(~acceptMask & ~phaseFillOut.trustedMask) = NaN;
% filledResidualMasked(phaseFillOut.trustedMask) = feat.phz(phaseFillOut.trustedMask);

%% Inpaint
phi = double(insarPair.phzReferenced);

phiCR = median(phi(crSupportMask), 'omitnan');
dphiCR = mod(phi - phiCR + pi, 2*pi) - pi;
sameWrapMask = abs(dphiCR) <= 1.0;

interpOpts = struct();
interpOpts.trainCohThresh = 0.5;
interpOpts.crSupportMask = crSupportMask;
interpOpts.directSolvedMask = directSolvedMask;
interpOpts.maxTrainDistToCR = inf;
interpOpts.sameWrapMask = sameWrapMask;

interpOpts.phaseSigma = 5;
interpOpts.cohSigma = 5;
interpOpts.fillCoherence = true;

interpOpts.qualitySigma = 5;
interpOpts.qualityMinValidFrac = 0.35;

out = insar.inpaint_interferogram_gaussian(insarPair, interpOpts);

phzFilled = out.phzFilled;
cohFilled = out.cohFilled;
phaseQuality = out.phaseQuality;
trustedMask = out.trustedMask;
%% Validation
valRand = insar.validate_phase_fill_model(feat, fillOpts, predOpts, 'random', 0.15);
valBlock = insar.validate_phase_fill_model(feat, fillOpts, predOpts, 'block', 50, 15);

% Report Validaiton
fprintf('\n--- Random holdout ---\n');
fprintf('Trusted: %d, Train: %d, Holdout: %d\n', ...
    valRand.nTrusted, valRand.nTrain, valRand.nHoldout);
fprintf('Raw RMSE: %.3f rad\n', valRand.metricsRaw.rmse);
fprintf('Raw MAE : %.3f rad\n', valRand.metricsRaw.mae);
fprintf('Raw Bias: %.3f rad\n', valRand.metricsRaw.bias);
fprintf('Raw Corr: %.3f\n', valRand.metricsRaw.corr);

for ii = 1:numel(valRand.metricsThresh)
    m = valRand.metricsThresh(ii);
    fprintf('Conf >= %.2f : retained=%.2f, RMSE=%.3f, MAE=%.3f, Bias=%.3f, Corr=%.3f\n', ...
        m.threshold, m.fracRetained, m.rmse, m.mae, m.bias, m.corr);
end

fprintf('\n--- Block holdout ---\n');
fprintf('Trusted: %d, Train: %d, Holdout: %d\n', ...
    valBlock.nTrusted, valBlock.nTrain, valBlock.nHoldout);
fprintf('Raw RMSE: %.3f rad\n', valBlock.metricsRaw.rmse);
fprintf('Raw MAE : %.3f rad\n', valBlock.metricsRaw.mae);
fprintf('Raw Bias: %.3f rad\n', valBlock.metricsRaw.bias);
fprintf('Raw Corr: %.3f\n', valBlock.metricsRaw.corr);

for ii = 1:numel(valBlock.metricsThresh)
    m = valBlock.metricsThresh(ii);
    fprintf('Conf >= %.2f : retained=%.2f, RMSE=%.3f, MAE=%.3f, Bias=%.3f, Corr=%.3f\n', ...
        m.threshold, m.fracRetained, m.rmse, m.mae, m.bias, m.corr);
end

% Plot Validation
figure();
imagesc(valBlock.holdoutMask);
axis image;
title('Block Holdout Mask');
colorbar;

figure();
imagesc(valBlock.pred.predictedResidual);
axis image;
title('Predicted Residual (Block Validation)');
colorbar;

figure();
imagesc(valBlock.pred.confidence);
axis image;
title('Confidence');
colorbar;

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

%% Benchmark Linear and Interpolation
baseRand = insar.validate_phase_fill_baselines( ...
    feat, valRand.holdoutMask, {'Z','slope','aspect_sin','aspect_cos','power','cor'}, 8);

baseBlock = insar.validate_phase_fill_baselines( ...
    feat, valBlock.holdoutMask, {'Z','slope','aspect_sin','aspect_cos','power','cor'}, 8);

fprintf('\n--- Random holdout baselines ---\n');
fprintf('Linear   RMSE: %.3f  MAE: %.3f  Bias: %.3f  Corr: %.3f\n', ...
    baseRand.linear.metrics.rmse, baseRand.linear.metrics.mae, ...
    baseRand.linear.metrics.bias, baseRand.linear.metrics.corr);
fprintf('Gaussian RMSE: %.3f  MAE: %.3f  Bias: %.3f  Corr: %.3f\n', ...
    baseRand.gaussian.metrics.rmse, baseRand.gaussian.metrics.mae, ...
    baseRand.gaussian.metrics.bias, baseRand.gaussian.metrics.corr);

fprintf('\n--- Block holdout baselines ---\n');
fprintf('Linear   RMSE: %.3f  MAE: %.3f  Bias: %.3f  Corr: %.3f\n', ...
    baseBlock.linear.metrics.rmse, baseBlock.linear.metrics.mae, ...
    baseBlock.linear.metrics.bias, baseBlock.linear.metrics.corr);
fprintf('Gaussian RMSE: %.3f  MAE: %.3f  Bias: %.3f  Corr: %.3f\n', ...
    baseBlock.gaussian.metrics.rmse, baseBlock.gaussian.metrics.mae, ...
    baseBlock.gaussian.metrics.bias, baseBlock.gaussian.metrics.corr);
%% Try a Sigma
    baseRand = insar.validate_phase_fill_baselines( ...
        feat, valRand.holdoutMask, {'Z','slope','aspect_sin','aspect_cos','power','cor'}, 5);
    z = exp(1i * baseRand.gaussian.predMap);
tmprewrapped = angle(z);
rewrapped = mod(tmprewrapped+ pi, 2*pi) - pi;
Q = insar.phase_quality_gaussian(rewrapped,51);
acceptMask = Q >= 0.9;
phz = baseRand.gaussian.predMap;
phz(~acceptMask) = NaN;

% Interpolate Phase and Coherence
phi = double(insarPair.phzReferenced);
cor = double(insarPair.cor);

% basic valid data mask
validMask = isfinite(phi) & isfinite(cor);

% coherence threshold
cohMask = cor >= 0.75;

% CR-referenced wrap-consistency
phiCR = median(phi(crSupportMask), 'omitnan');
dphiCR = mod(phi - phiCR + pi, 2*pi) - pi;
sameWrapMask = abs(dphiCR) <= 1.0;   % radians, tune

% trusted training mask
trustedMask = validMask & cohMask & sameWrapMask;
philledPhi = insar.phase_fill_gaussian_interp(phi, feat.trustedMask, 5);
philledCoh = insar.gaussian_interp_coherence(cor, feat.trustedMask, 5);
%% Gaussian Sigma Sweep
sigmas = [2 3 4 5 6];

randResults  = struct([]);
blockResults = struct([]);

for ii = 1:numel(sigmas)
    sigma = sigmas(ii);

    baseRand = insar.validate_phase_fill_baselines( ...
        feat, valRand.holdoutMask, {'Z','slope','aspect_sin','aspect_cos','power','cor','incidence'}, sigma);

    baseBlock = insar.validate_phase_fill_baselines( ...
        feat, valBlock.holdoutMask, {'Z','slope','aspect_sin','aspect_cos','power','cor'}, sigma);

    randResults(ii).sigma = sigma;
    randResults(ii).rmse  = baseRand.gaussian.metrics.rmse;
    randResults(ii).mae   = baseRand.gaussian.metrics.mae;
    randResults(ii).bias  = baseRand.gaussian.metrics.bias;
    randResults(ii).corr  = baseRand.gaussian.metrics.corr;

    blockResults(ii).sigma = sigma;
    blockResults(ii).rmse  = baseBlock.gaussian.metrics.rmse;
    blockResults(ii).mae   = baseBlock.gaussian.metrics.mae;
    blockResults(ii).bias  = baseBlock.gaussian.metrics.bias;
    blockResults(ii).corr  = baseBlock.gaussian.metrics.corr;

    fprintf('sigma=%2d | Rand RMSE=%.3f MAE=%.3f Corr=%.3f | Block RMSE=%.3f MAE=%.3f Corr=%.3f\n', ...
        sigma, ...
        randResults(ii).rmse, randResults(ii).mae, randResults(ii).corr, ...
        blockResults(ii).rmse, blockResults(ii).mae, blockResults(ii).corr);
end

T = table( ...
    [blockResults.sigma]', ...
    [blockResults.rmse]', ...
    [blockResults.mae]', ...
    [blockResults.bias]', ...
    [blockResults.corr]', ...
    'VariableNames', {'sigma','RMSE','MAE','Bias','Corr'});

disp(T)

figure;
plot([blockResults.sigma], [blockResults.rmse], '-o');
xlabel('\sigma');
ylabel('Block holdout RMSE (rad)');
grid on;
title('Gaussian interpolation sigma sweep');
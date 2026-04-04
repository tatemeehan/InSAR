function [insarData, unwrapOpts] = process_interferometric_phase(sarData, geomData, demData, crResult, params, unwrapOpts)
%PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% Reference hierarchy is now sensor-perspective and CR-driven:
%   0 valid CRs -> sceneMean
%   1 valid CR  -> crConstant
%   2 valid CRs -> crAzRamp, crRangeRamp, or crAvgRamp (auto-selected or forced)
%   3+ valid CRs -> crSensorPlane
%
% Optional PCA phase-screen removal remains available as a pre-reference step.

if nargin < 5, unwrapOpts = struct(); end

% -------------------------
% Defaults
% -------------------------
if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
if ~isfield(unwrapOpts, 'useRobustWeighting'), unwrapOpts.useRobustWeighting = false; end
if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
if ~isfield(unwrapOpts, 'postprocessRegions'), unwrapOpts.postprocessRegions = false; end
if ~isfield(unwrapOpts, 'stitchMethod'), unwrapOpts.stitchMethod = 'mst'; end
if ~isfield(unwrapOpts, 'tinyRegionSize'), unwrapOpts.tinyRegionSize = 30; end
if ~isfield(unwrapOpts, 'stitchAlpha'), unwrapOpts.stitchAlpha = 0.5; end
if ~isfield(unwrapOpts, 'crSlcPolicy'), unwrapOpts.crSlcPolicy = 'burst'; end

% Multiseed Bucket Unwrap Options
if ~isfield(unwrapOpts,'bucketBins'), unwrapOpts.bucketBins = 256; end
if ~isfield(unwrapOpts,'seedQualityThresh'), unwrapOpts.seedQualityThresh = 0.75; end
if ~isfield(unwrapOpts,'validQualityThresh'), unwrapOpts.validQualityThresh = unwrapOpts.qualityThresh; end
if ~isfield(unwrapOpts,'minSeedRegionSize'), unwrapOpts.minSeedRegionSize = 1; end
% Phase Region Alignment Options
if ~isfield(unwrapOpts,'reconcileBasins'), unwrapOpts.reconcileBasins = true; end
if ~isfield(unwrapOpts,'reconcileGapBasins'), unwrapOpts.reconcileGapBasins = true; end
if ~isfield(unwrapOpts,'basinMinEdgeSupport'), unwrapOpts.basinMinEdgeSupport = 5; end
if ~isfield(unwrapOpts,'basinGapMaxNeighbors'),unwrapOpts.basinGapMaxNeighbors = 2; end
if ~isfield(unwrapOpts,'basinGapMaxPixels'), unwrapOpts.basinGapMaxPixels = 25; end
if ~isfield(unwrapOpts,'basinUseReferenceScreen'), unwrapOpts.basinUseReferenceScreen = true; end
if ~isfield(unwrapOpts,'basinScreenCohThresh'), unwrapOpts.basinScreenCohThresh = max(0.4, unwrapOpts.qualityThresh); end
if ~isfield(unwrapOpts,'basinAnchorPolicy'), unwrapOpts.basinAnchorPolicy = 'crFirst'; end
if ~isfield(unwrapOpts,'storeUnwrapMeta'), unwrapOpts.storeUnwrapMeta = true; end

% -------------------------
% ML / Gaussian phase inpainting defaults
% -------------------------
if ~isfield(unwrapOpts, 'phaseFill')
    unwrapOpts.phaseFill = struct();
end

if ~isfield(unwrapOpts.phaseFill, 'enable'), unwrapOpts.phaseFill.enable = false; end
if ~isfield(unwrapOpts.phaseFill, 'method'), unwrapOpts.phaseFill.method = 'ml'; end   % 'ml' or 'gaussian'

% --- shared support defaults
if ~isfield(unwrapOpts.phaseFill, 'trainCohThresh'), unwrapOpts.phaseFill.trainCohThresh = 0.75; end
if ~isfield(unwrapOpts.phaseFill, 'maxTrainDistToCR'), unwrapOpts.phaseFill.maxTrainDistToCR = 250; end
if ~isfield(unwrapOpts.phaseFill, 'sameWrapThresh'), unwrapOpts.phaseFill.sameWrapThresh = 1.0; end

% --- CR support graph defaults
if ~isfield(unwrapOpts.phaseFill, 'crRadiusPx'), unwrapOpts.phaseFill.crRadiusPx = 100; end
if ~isfield(unwrapOpts.phaseFill, 'maxGraphHops'), unwrapOpts.phaseFill.maxGraphHops = 50; end
if ~isfield(unwrapOpts.phaseFill, 'includeGapEdges'), unwrapOpts.phaseFill.includeGapEdges = false; end

% --- ML defaults
if ~isfield(unwrapOpts.phaseFill, 'maxTrainSamples'), unwrapOpts.phaseFill.maxTrainSamples = 100000; end
if ~isfield(unwrapOpts.phaseFill, 'modelType'), unwrapOpts.phaseFill.modelType = 'bag'; end
if ~isfield(unwrapOpts.phaseFill, 'numTrees'), unwrapOpts.phaseFill.numTrees = 250; end
if ~isfield(unwrapOpts.phaseFill, 'minLeafSize'), unwrapOpts.phaseFill.minLeafSize = 15; end
if ~isfield(unwrapOpts.phaseFill, 'minValidFrac'), unwrapOpts.phaseFill.minValidFrac = 0.25; end
if ~isfield(unwrapOpts.phaseFill, 'learnRate'), unwrapOpts.phaseFill.learnRate = 0.05; end

% best working feature set, incidence optional
if ~isfield(unwrapOpts.phaseFill, 'selectedFeatures')
    unwrapOpts.phaseFill.selectedFeatures = {'Z','V','slope','aspect_sin','aspect_cos','power','cor','incidence'};
end

% --- prediction defaults
if ~isfield(unwrapOpts.phaseFill, 'predCohThresh'), unwrapOpts.phaseFill.predCohThresh = 0.25; end
if ~isfield(unwrapOpts.phaseFill, 'maxFillDistance'), unwrapOpts.phaseFill.maxFillDistance = 500; end
if ~isfield(unwrapOpts.phaseFill, 'maxPredictDistToCR'), unwrapOpts.phaseFill.maxPredictDistToCR = inf; end
if ~isfield(unwrapOpts.phaseFill, 'maxPredictDistToDirectSolved'), unwrapOpts.phaseFill.maxPredictDistToDirectSolved = 150; end
if ~isfield(unwrapOpts.phaseFill, 'smoothSigma'), unwrapOpts.phaseFill.smoothSigma = 3; end
if ~isfield(unwrapOpts.phaseFill, 'preserveTrusted'), unwrapOpts.phaseFill.preserveTrusted = true; end

% --- quality defaults
if ~isfield(unwrapOpts.phaseFill, 'qualitySigma'), unwrapOpts.phaseFill.qualitySigma = 5; end
if ~isfield(unwrapOpts.phaseFill, 'qualityMinValidFrac'), unwrapOpts.phaseFill.qualityMinValidFrac = 0.35; end
if ~isfield(unwrapOpts.phaseFill, 'qualityThreshold'), unwrapOpts.phaseFill.qualityThreshold = []; end   % e.g. 0.7 or []
if ~isfield(unwrapOpts.phaseFill, 'storeModelInfo'), unwrapOpts.phaseFill.storeModelInfo = false; end

% Pair construction
if ~isfield(unwrapOpts,'pairingMode'), unwrapOpts.pairingMode = 'inter'; end
pairOpts = struct;
switch lower(unwrapOpts.pairingMode)
    case 'inter'
        pairOpts.includeInter = true;
        pairOpts.includeIntra = false;
        pairOpts.burstMode = 'sameIndex';
    case 'intra'
        pairOpts.includeInter = false;
        pairOpts.includeIntra = true;
        pairOpts.burstMode = 'allCombos';
    otherwise
        pairOpts.includeInter = true;
        pairOpts.includeIntra = true;
        pairOpts.burstMode = 'allCombos';
end
pairOpts.filterFn = @(i1,b1,i2,b2) true;

% Optional pre-reference screen removal (PCA, etc.)
if ~isfield(unwrapOpts, 'preReferenceScreenRemoval')
    unwrapOpts.preReferenceScreenRemoval = struct();
end
if ~isfield(unwrapOpts.preReferenceScreenRemoval, 'enable')
    unwrapOpts.preReferenceScreenRemoval.enable = false;
end
if ~isfield(unwrapOpts.preReferenceScreenRemoval, 'method')
    unwrapOpts.preReferenceScreenRemoval.method = 'pcaAz';
end
if ~isfield(unwrapOpts.preReferenceScreenRemoval, 'cohThresh')
    unwrapOpts.preReferenceScreenRemoval.cohThresh = 0.3;
end

% Reference options
if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'auto'; end
if ~isfield(unwrapOpts, 'referenceFallback'), unwrapOpts.referenceFallback = 'sceneMean'; end
if ~isfield(unwrapOpts, 'referenceImgForCR'), unwrapOpts.referenceImgForCR = 'both'; end
if ~isfield(unwrapOpts, 'referenceMinCoherence'), unwrapOpts.referenceMinCoherence = 0.3; end
if ~isfield(unwrapOpts, 'referenceUseRobustFit'), unwrapOpts.referenceUseRobustFit = true; end
if ~isfield(unwrapOpts, 'referenceAutoTwoCRMode'), unwrapOpts.referenceAutoTwoCRMode = 'score'; end
if ~isfield(unwrapOpts, 'referencePreferSensorCoords'), unwrapOpts.referencePreferSensorCoords = true; end

% Penetration Options
if ~isfield(unwrapOpts, 'targetSensitivity'), unwrapOpts.targetSensitivity = 5; end
if ~isfield(unwrapOpts, 'phaseNoiseStd') || isempty(unwrapOpts.phaseNoiseStd),unwrapOpts.phaseNoiseStd = 0.1; end
if ~isfield(unwrapOpts, 'estimatePhaseNoiseFromCoherence'), unwrapOpts.estimatePhaseNoiseFromCoherence = true; end
if ~isfield(unwrapOpts, 'nLooksEffective'),unwrapOpts.nLooksEffective = 1; end
if ~isfield(unwrapOpts.phaseFill, 'useFilledForPenetration')
    unwrapOpts.phaseFill.useFilledForPenetration = true;
end
if ~isfield(unwrapOpts.phaseFill, 'useMaskedFilledForPenetration')
    unwrapOpts.phaseFill.useMaskedFilledForPenetration = false;
end


pairs = utils.build_pairs(sarData, pairOpts);
if isempty(pairs)
    warning('No interferometric pairs formed.');
    insarData = struct([]);
    return;
end

fprintf('Interferogram Queue (%d total):\n', numel(pairs));
for p = 1:numel(pairs)
    fprintf('  (%d,b%02d) × (%d,b%02d)\n', pairs(p).i1, pairs(p).b1, pairs(p).i2, pairs(p).b2);
end

% -------------------------
% Main loop
% -------------------------
insarData = struct([]);
k = 0;

for p = 1:numel(pairs)
    tic

    % clear pair-local vars
    phzCorrected = [];
    cycleMap = [];
    correctionMap = [];
    regionLabels = [];
    unwrapMeta = struct();
    basinMeta = struct();

    i  = pairs(p).i1;
    j  = pairs(p).i2;
    b1 = pairs(p).b1;
    b2 = pairs(p).b2;

    % -------------------------
    % Pair geometry
    % -------------------------
    idxG = utils.lookup_pair_geom(geomData, i, b1, j, b2);
    if isnan(idxG)
        warning('No geometry for (%d,b%02d)×(%d,b%02d). Skipping.', i, b1, j, b2);
        continue;
    end
    G = geomData.sarGeometry{idxG};

    slc1 = sarData(i).slc{b1};
    slc2 = sarData(j).slc{b2};

    if ~isequal(size(slc1), size(slc2))
        warning('Skipping (%d,b%02d)×(%d,b%02d): SLC size mismatch [%s] vs [%s].', ...
            i, b1, j, b2, mat2str(size(slc1)), mat2str(size(slc2)));
        continue;
    end

    % -------------------------
    % Mask
    % -------------------------
    mask = true(size(slc1));
    if isfield(G,'lookMask') && ~isempty(G.lookMask)
        mask = mask & logical(G.lookMask);
    end
    slc1(~mask) = NaN;
    slc2(~mask) = NaN;

    % -------------------------
    % Interferogram
    % -------------------------
    [phzWrapped, cor, ccor, cphz, meta] = insar.compute_interferogram( ...
        slc1, slc2, unwrapOpts.qualityThresh, unwrapOpts.filterSize, ...
        unwrapOpts.useRobustWeighting, unwrapOpts.sigma, unwrapOpts.alpha);

    phzWrapped(~mask) = NaN;
    cor(~mask)        = NaN;
    ccor(~mask)       = NaN;
    cphz(~mask)       = NaN;
    unwrapOpts.nLooksEffective = meta.nLooksEffective;

    % -------------------------
    % Optional pre-reference screen removal
    % -------------------------
    preRefModel = [];
    if unwrapOpts.preReferenceScreenRemoval.enable
        switch lower(unwrapOpts.preReferenceScreenRemoval.method)
            case 'pcaaz'
                v = diff(sarData(i).traj{b1}(:,1:2),1,1);
                vm = mean(v,1,'omitnan');
                headingDeg = rad2deg(atan2(vm(2), vm(1)));

                rampOpts = struct();
                rampOpts.cohThresh = unwrapOpts.preReferenceScreenRemoval.cohThresh;

                [cphz, preRefModel] = insar.remove_phzramp_azpca_rot(cphz, cor, headingDeg, rampOpts);
                phzWrapped = angle(cphz);

            otherwise
                error('Unknown preReferenceScreenRemoval.method: %s', ...
                    unwrapOpts.preReferenceScreenRemoval.method);
        end
    end

    % -------------------------
    % Unwrap
    % -------------------------
    switch lower(unwrapOpts.method)
        case 'multiseed'
            % phzUnwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor);
            [phzUnwrapped, unwrapMeta] = insar.region_grow_unwrap_bucket( ...
                phzWrapped, cor, mask, unwrapOpts);

            if unwrapOpts.postprocessRegions
                [phzCorrected, cycleMap, correctionMap, regionLabels] = ...
                    insar.correct_phase_cycles_regionwise(phzUnwrapped, cor, unwrapOpts.tinyRegionSize);
            end

        case 'goldstein'
            phzUnwrapped = insar.unwrap_goldstein(phzWrapped);

        case 'none'
            phzUnwrapped = phzWrapped;

        otherwise
            error('Unknown unwrap method: %s', unwrapOpts.method);
    end

    % -------------------------
    % Build sensor-perspective predictors
    % -------------------------
    pred = insar.build_reference_predictors(G, mask);

    % Align Unwrapped Basins with Phase Screen Predictors
    if strcmpi(unwrapOpts.method, 'multiseed') && unwrapOpts.reconcileBasins && ...
            isfield(unwrapMeta, 'owner') && ~isempty(unwrapMeta.owner)
        [phzUnwrapped, basinMeta] = insar.align_unwrap_basins_with_reference( ...
            phzWrapped, phzUnwrapped, cor, unwrapMeta.owner, ...
            crResult, pred, i, b1, j, b2, unwrapOpts);
    end

    % -------------------------
    % CR-based spatial reference screen
    % -------------------------
    [refScreen, refScreenWrapped, refMeta] = insar.build_phase_reference_screen_sensor( ...
        phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, ...
        i, b1, j, b2, unwrapOpts);

    phzReferenced = phzUnwrapped - refScreen;
    phzWrappedReferenced = angle(exp(1i * (phzWrapped - refScreenWrapped)));

    % -------------------------
    % ML/Gaussian interferogram inpainting
    % -------------------------
    phaseFill = struct();
    if unwrapOpts.phaseFill.enable
        insarPair = struct();
        insarPair.phzReferenced = phzReferenced;
        insarPair.cor = cor;

        if isfield(sarData(i), 'db') && numel(sarData(i).db) >= b1
            insarPair.power = sarData(i).db{b1};
        end

        basinOwnerLocal = [];
        if isfield(unwrapMeta, 'owner') && ~isempty(unwrapMeta.owner)
            basinOwnerLocal = unwrapMeta.owner;
        end

        phaseFill = insar.inpaint_interferogram_ml( ...
            insarPair, G, basinOwnerLocal, basinMeta, crResult, demData, unwrapOpts.phaseFill);
    end

    % -------------------------
    % Store
    % -------------------------
    k = k + 1;
    insarData(k).pair = [i j];
    insarData(k).burst = [b1 b2];

    insarData(k).phzWrapped = phzWrapped;
    insarData(k).phzWrappedReferenced = phzWrappedReferenced;
    insarData(k).phzUnwrapped = phzUnwrapped;
    insarData(k).phzReferenced = phzReferenced;

    insarData(k).coherence = cor;
    insarData(k).complexCoherence = ccor;
    insarData(k).complexPhase = cphz;

    if unwrapOpts.postprocessRegions && ~isempty(phzCorrected)
        insarData(k).phzCorrected = phzCorrected;
        insarData(k).phzCycleMap = cycleMap;
        insarData(k).phzCorrectionMap = correctionMap;
        insarData(k).regionLabels = regionLabels;
    end

    insarData(k).referencePredictors = pred;
    insarData(k).refScreen = refScreen;
    insarData(k).refScreenWrapped = refScreenWrapped;
    insarData(k).refMeta = refMeta;
    insarData(k).refMethodRequested = unwrapOpts.referenceMethod;
    insarData(k).refMethodUsed = refMeta.methodUsed;
    insarData(k).refValue = refMeta.scalarEquivalent;
    insarData(k).refValueWrapped = refMeta.scalarEquivalentWrapped;
    insarData(k).unwrapMethod = unwrapOpts.method;
    insarData(k).preReferenceScreenRemoval = preRefModel;

    if unwrapOpts.storeUnwrapMeta
        insarData(k).unwrapMeta = unwrapMeta;
        insarData(k).basinMeta = basinMeta;

        if isfield(unwrapMeta, 'owner')
            insarData(k).basinOwner = unwrapMeta.owner;
        end
        if isfield(unwrapMeta, 'phzLocal')
            insarData(k).phzUnwrappedLocal = unwrapMeta.phzLocal;
        end
        if isfield(basinMeta, 'shift')
            insarData(k).basinShift = basinMeta.shift;
        end
        if isfield(basinMeta, 'anchorBasin')
            insarData(k).basinAnchor = basinMeta.anchorBasin;
        end
        if isfield(basinMeta, 'solved')
            insarData(k).basinSolved = basinMeta.solved;
        end
    end

    if unwrapOpts.phaseFill.enable
        insarData(k).phzReferencedFilled = phaseFill.phzFilled;
        insarData(k).phzWrappedReferencedFilled = phaseFill.phzWrappedFilled;
        insarData(k).phaseQuality = phaseFill.phaseQuality;  % synthetic coherence metric
        insarData(k).phaseQualitySupportFrac = phaseFill.qualitySupportFrac;
        insarData(k).phaseFillTrustedMask = phaseFill.trustedMask;
        insarData(k).phaseFillCRSupportMask = phaseFill.crSupportMask;
        insarData(k).phaseFillSameWrapMask = phaseFill.sameWrapMask;
        insarData(k).phaseFillSupport = phaseFill.support;

        if ~isempty(phaseFill.acceptMask)
            insarData(k).phzReferencedFilledMasked = phaseFill.phzFilledMasked;
            insarData(k).phaseFillAcceptMask = phaseFill.acceptMask;
        end

        if unwrapOpts.phaseFill.storeModelInfo && isfield(phaseFill, 'modelInfo')
            insarData(k).phaseFillModelInfo = phaseFill.modelInfo;
        end
    end
    % -------------------------
    % Phase source for penetration
    % -------------------------
    phzForPenetration = phzReferenced;
    phzForPenetrationSource = 'referenced';

    if unwrapOpts.phaseFill.enable && isfield(unwrapOpts.phaseFill, 'useFilledForPenetration') && ...
            unwrapOpts.phaseFill.useFilledForPenetration

        if isfield(unwrapOpts.phaseFill, 'useMaskedFilledForPenetration') && ...
                unwrapOpts.phaseFill.useMaskedFilledForPenetration && ...
                isfield(phaseFill, 'phzFilledMasked') && ~isempty(phaseFill.phzFilledMasked)

            phzForPenetration = phaseFill.phzFilledMasked;
            phzForPenetrationSource = 'filledMasked';

        elseif isfield(phaseFill, 'phzFilled') && ~isempty(phaseFill.phzFilled)

            phzForPenetration = phaseFill.phzFilled;
            phzForPenetrationSource = 'filled';
        end
    end

    % -------------------------
    % Penetration
    % -------------------------
    if isfield(G,'Bperp') && ~isempty(G.Bperp)

        % --- Pair-average slant range ---
        if isfield(G,'slant') && ~isempty(G.slant) && isfield(G,'slant2') && ~isempty(G.slant2)
            rslant = 0.5 * (G.slant + G.slant2);
        elseif isfield(G,'slant') && ~isempty(G.slant)
            rslant = G.slant;
        elseif isfield(G,'slant2') && ~isempty(G.slant2)
            rslant = G.slant2;
        else
            rslant = [];
        end

        % --- Pair-average incidence angle ---
        if isfield(G,'incidence') && ~isempty(G.incidence) && isfield(G,'incidence2') && ~isempty(G.incidence2)
            inc = 0.5 * (G.incidence + G.incidence2);
        elseif isfield(G,'incidence') && ~isempty(G.incidence)
            inc = G.incidence;
        elseif isfield(G,'incidence2') && ~isempty(G.incidence2)
            inc = G.incidence2;
        else
            inc = [];
        end

        Bperp = G.Bperp;

        if isempty(rslant) || isempty(inc)
            warning('Missing slant/incidence fields for geom idx %d. Skipping penetration.', idxG);
            insarData(k).penetration = single(NaN(size(phzForPenetration), 'single'));
            insarData(k).penetrationValid = false;
            insarData(k).penetrationValidMask = false(size(phzForPenetration));
            insarData(k).penetrationSensitivity = single(NaN(size(phzForPenetration), 'single'));      % m/rad
            insarData(k).penetrationDetectableHeight = single(NaN(size(phzForPenetration), 'single')); % m
            continue
        end

    % Degrees -> radians guard
    if any(inc(:) > pi, 'all')
        inc = deg2rad(inc);
    end

    target_sensitivity = unwrapOpts.targetSensitivity; % m/rad

    % Phase noise for minimum detectable height map
    if unwrapOpts.estimatePhaseNoiseFromCoherence
        phase_noise_std = insar.coherence_to_phase_std(cor, unwrapOpts.nLooksEffective);
    else
        phase_noise_std = unwrapOpts.phaseNoiseStd .* ones(size(cor), 'like', cor);
    end

    sin_theta = sin(inc);
    sin_theta(~isfinite(sin_theta) | abs(sin_theta) < 1e-6) = eps;

    absBperp = abs(Bperp);
    absBperp(~isfinite(absBperp) | absBperp < 1e-6) = NaN;

    % --- Sensitivity map: meters per radian ---
    sensitivity_map = (params.lambda .* rslant .* sin_theta) ./ ...
                      (4 .* pi .* absBperp);

    % --- Baseline threshold required to achieve target sensitivity ---
    Bperp_thresh = (params.lambda .* rslant) ./ ...
                   (4 .* pi .* target_sensitivity .* sin_theta);

    % --- Pixelwise validity mask ---
    validMask = isfinite(phzForPenetration) & isfinite(sensitivity_map) & ...
        isfinite(Bperp_thresh) & isfinite(abs(Bperp)) & ...
        (abs(Bperp) >= Bperp_thresh);

    % --- Minimum detectable height / penetration change map ---
    detectable_height_map = sensitivity_map .* phase_noise_std;

    % --- Summary stats ---
    mean_thresh = median(Bperp_thresh(:), 'omitnan');
    mean_Bperp  = median(abs(Bperp(:)), 'omitnan');
    frac_valid  = mean(validMask(:), 'omitnan');

    % --- Penetration estimate ---
    penetration = nan(size(phzForPenetration), 'single');

    pen_valid = insar.phase_to_penetration( ...
        phzForPenetration(validMask), params.lambda, Bperp(validMask), ...
        rslant(validMask), inc(validMask));

    penetration(validMask) = single(pen_valid);

    % --- Store compactly ---
    insarData(k).penetration = penetration;
    insarData(k).penetrationValid = any(validMask(:));
    insarData(k).penetrationValidMask = validMask;
    insarData(k).penetrationPhaseSource = phzForPenetrationSource;

    insarData(k).penetrationSensitivity = single(sensitivity_map);            % m/rad
    insarData(k).penetrationDetectableHeight = single(detectable_height_map); % m
    insarData(k).meanPenetrationSensitivity = median(sensitivity_map(:), 'omitnan');
    insarData(k).meanDetectableHeight = median(detectable_height_map(:), 'omitnan');

    insarData(k).meanBperp = mean_Bperp;
    insarData(k).minRequiredBperp = mean_thresh;
    insarData(k).sensitivityTarget = target_sensitivity;
    insarData(k).phaseNoiseStd = single(phase_noise_std);
    insarData(k).validSensitivityFraction = frac_valid;
    insarData(k).geomIdx = idxG;

    if ~any(validMask(:))
        warning(['No pixels meet target sensitivity %.2f m/rad ' ...
                 '(pair %d-%d bursts [%d %d], geom %d)'], ...
                 target_sensitivity, i, j, b1, b2, idxG);
    end

else
    warning('Missing geometry fields for geom idx %d. Skipping penetration.', idxG);
    insarData(k).penetration = single(NaN(size(phzForPenetration), 'single'));
    insarData(k).penetrationValid = false;
    insarData(k).penetrationValidMask = false(size(phzForPenetration));
    insarData(k).penetrationSensitivity = single(NaN(size(phzForPenetration), 'single'));
    insarData(k).penetrationDetectableHeight = single(NaN(size(phzForPenetration), 'single'));
end

    fprintf('Processed (dir %d,b%02d) × (dir %d,b%02d). Elapsed: %.2f s\n', ...
        i, b1, j, b2, toc);
end
end
% function [insarData, unwrapOpts] = process_interferometric_phase(sarData, geomData, crResult, params, unwrapOpts)
% %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% %   Includes multiseed unwrapping with optional post-processing upgrades:
% %       - Region bridging (cross-region stitching)
% %       - Phase gradient blending
% %       - MST-based merging
% %       - Fallback for tiny regions
% 
% if nargin < 5, unwrapOpts = struct(); end
% if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'none'; end
% if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% if ~isfield(unwrapOpts, 'useRobustWeighting'), unwrapOpts.useRobustWeighting = false; end
% if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize ./ 2; end
% if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
% if ~isfield(unwrapOpts, 'rampRemoval'), unwrapOpts.rampRemoval.enable = false; end
% if ~isfield(unwrapOpts, 'postprocessRegions'), unwrapOpts.postprocessRegions = false; end
% if ~isfield(unwrapOpts, 'stitchMethod'), unwrapOpts.stitchMethod = 'mst'; end
% if ~isfield(unwrapOpts, 'tinyRegionSize'), unwrapOpts.tinyRegionSize = 30; end
% if ~isfield(unwrapOpts, 'stitchAlpha'), unwrapOpts.stitchAlpha = 0.5; end
% if ~isfield(unwrapOpts,'crSlcPolicy'), unwrapOpts.crSlcPolicy = 'burst'; end
% % Pair construction (flexible)
% if ~isfield(unwrapOpts,'pairingMode'), unwrapOpts.pairingMode = 'inter'; end
% pairOpts = struct;
% switch lower(unwrapOpts.pairingMode)
%     case 'inter'
%         pairOpts.includeInter = true;  pairOpts.includeIntra = false;
%         pairOpts.burstMode = 'sameIndex';
%     case 'intra'
%         pairOpts.includeInter = false; pairOpts.includeIntra = true;
%         pairOpts.burstMode = 'allCombos';  % or 'sameIndex' neighbor-only
%     otherwise % 'auto' -> do both
%         pairOpts.includeInter = true;
%         pairOpts.includeIntra = true;
%         pairOpts.burstMode = 'allCombos';  % common choice; change if you want denser pairing
% end
% 
% % Optional: attach a filter predicate (e.g., geometry availability)
% pairOpts.filterFn = @(i1,b1,i2,b2) true;
% 
% pairs = utils.build_pairs(sarData, pairOpts);
% if isempty(pairs)
%     warning('No interferometric pairs formed.'); insarData = struct([]); return;
% end
% 
% % Print InSAR Queue
% fprintf('Interferogram Queue (%d total):\n', numel(pairs));
% for p = 1:numel(pairs)
%     fprintf('  (%d,b%02d) × (%d,b%02d)\n', pairs(p).i1, pairs(p).b1, pairs(p).i2, pairs(p).b2);
% end
% 
% 
% % === Main processing loop over pairs ===
% insarData = struct([]);
% k = 0;
% 
% for p = 1:numel(pairs)
%     i = pairs(p).i1; j = pairs(p).i2;
%     b1 = pairs(p).b1; b2 = pairs(p).b2;
%     tic
% 
%     % Apply lookMask using geomData
%     idxG = utils.lookup_pair_geom(geomData, i, b1, j, b2);
%     % idxG = utils.lookup_or_build_pair_geom(geomData, i,b1, j,b2, @buildGeomPair);
%     if isnan(idxG)
%         warning('No geometry for (%d,b%02d)×(%d,b%02d). Skipping.', i,b1,j,b2);
%         continue;
%     end
%     G = geomData.sarGeometry{idxG};
% 
%     slc1 = sarData(i).slc{b1};
%     slc2 = sarData(j).slc{b2};
% 
%     mask = true(size(slc1));
%     if isfield(G,'lookMask') && ~isempty(G.lookMask)
%         mask = mask & logical(G.lookMask);
%     end
%     slc1(~mask) = NaN;  slc2(~mask) = NaN;
% 
%     % mask = true(size(slc1));
%     % if isfield(geomData,"slcGeomIndex") && isfield(geomData,"sarGeometry")
%     %     idx1 = geomData.slcGeomIndex(i, b1).idx;
%     %     idx2 = geomData.slcGeomIndex(j, b2).idx;
%     %     if ~isnan(idx1), mask = mask & geomData.sarGeometry{idx1}.lookMask; end
%     %     if ~isnan(idx2), mask = mask & geomData.sarGeometry{idx2}.lookMask; end
%     % end
%     % slc1(~mask) = NaN; slc2(~mask) = NaN;
%     % if isfield(geomData, "slcGeomIndex") && isfield(geomData, "sarGeometry")
%         % idx1 = geomData.slcGeomIndex(i, b1).idx;
%         % idx2 = geomData.slcGeomIndex(j, b2).idx;
%     %     if ~isnan(idx1) && ~isnan(idx2) && idx1 == idx2
%     %         mask = geomData.sarGeometry{idx1}.lookMask;
%     %         slc1(~mask) = NaN;
%     %         slc2(~mask) = NaN;
%     %     end
%     % end
% 
%     % After loading slc1/slc2: Check for Compatibility
%     if ~isequal(size(slc1), size(slc2))
%         warning('Skipping (%d,b%02d)×(%d,b%02d): SLC size mismatch [%s] vs [%s].', ...
%             i,b1,j,b2, mat2str(size(slc1)), mat2str(size(slc2)));
%         continue;
%     end
% 
%     % Compute Multi-Looked Interferograms
%     [phzWrapped, cor, ccor, cphz] = insar.compute_interferogram(slc1, slc2, ...
%         unwrapOpts.qualityThresh, unwrapOpts.filterSize, unwrapOpts.useRobustWeighting, ...
%         unwrapOpts.sigma, unwrapOpts.alpha);
%     % Apply Look Mask
%     phzWrapped(~mask) = NaN; cor(~mask) = NaN; ccor(~mask) = NaN; cphz(~mask) = NaN;
% 
%     % Remove Phase Ramp
%     if isfield(unwrapOpts,'rampRemoval') && unwrapOpts.rampRemoval.enable
%         v   = diff(sarData(i).traj{b1}(:,1:2),1,1);          % horizontal velocity samples
%         vm  = mean(v,1,'omitnan');                  % mean along-track vector
%         headingDeg = rad2deg(atan2(vm(2), vm(1))); % 0°=East, 90°=Nort
%         rampOpts.cohThresh = 0.3;
%         [cphz, ~] = insar.remove_phzramp_azpca_rot(cphz, cor, headingDeg,rampOpts);
%         phzWrapped = angle(cphz);
%         % RR = unwrapOpts.rampRemoval;   % options (sigmaPred, sigmaPhi, cohThresh, robust, use)
%         % % Build predictor struct from geometry you have for THIS pair/burst:
%         % geomPred = struct();
%         % geomPred.dR = geomData.sarGeometry{idx1}.slant - geomData.sarGeometry{idx2}.slant2;
%         % geomPred.dInc = geomData.sarGeometry{idx1}.incidence - geomData.sarGeometry{idx2}.incidence2;
%         % geomPred.Bperp = geomData.sarGeometry{idx1}.Bperp;
%         % geomPred.lookMask = mask;%geomData.sarGeometry{idx1}.lookMask;
% 
%         % % Provide XY to catch residual planar ramps (optional, function will add if missing):
%         % % [X,Y] = meshgrid(1:size(Ifilt,2), 1:size(Ifilt,1));
%         % geomPred.X = demData.X; geomPred.Y = demData.Y; geomPred.Z = demData.dem;
%         % geomPred.closestIndex_master = geomData.sarGeometry{idx1}.closestIndex_master;
%         % geomPred.closestIndex_slave = geomData.sarGeometry{idx2}.closestIndex_slave;
%         % geomPred.height = demData.dem; 
%         % geomPred.slant = geomData.sarGeometry{idx1}.slant;
%         % geomPred.slant2 = geomData.sarGeometry{idx2}.slant2;
%         % geomPred.incidence = geomData.sarGeometry{idx1}.incidence;
%         %  geomPred.incidence2 = geomData.sarGeometry{idx2}.incidence2;
%          % G.lambda = params.lambda;
%         % out = insar.remove_phase_waves_geomdiff(phzWrapped, cor, demData.X, demData.Y, demData.dem, G, struct('cohThresh',0.6,'usePlane',false))
%         % trajM.pos = sarData(i).traj{b1};
%         % trajS.pos = sarData(j).traj{b2};
%         % [I_corr, model] = insar.remove_phzramp_geometry_gammaColumn(ccor, cor, geomPred, trajM, trajS, params.lambda);
%         % [I_corr, model] = insar.remove_phzramp_deltaB_column(ccor, cor, geomPred, trajM, params.lambda);
%         % [I_corr, model] = insar.remove_phzramp_geometry(ccor, cor, geomPred, trajM, trajS, params.lambda);
%         % [Icorr, mdl] = insar.remove_phzramp_rangeRidge(ccor, cor, geomPred.slant);
%         % [Icorr, model] = insar.remove_phzramp_sinusoid(ccor, cor, 10);
%         % [Icorr, model] = insar.remove_phzramp_azRamp(ccor, cor);
%         % [Icorr, model] = insar.remove_phzramp_geomScaleCol(ccor, cor, geomPred.dR);
%         % [Icorr, model] = insar.remove_phzramp_geom2col(ccor, cor, geomPred.dR, deg2rad(geomPred.dInc));
%         % [Icorr, model] = insar.remove_phzramp_azpca(ccor, cor);
% 
% %         [IfiltCorr, rampModel] = insar.remove_phzramp_multivar(ccor, cor, geomPred, RR);
% %         % Use corrected wrapped phase for unwrapping:
% %         ccor = IfiltCorr;
% %         cor = abs(IfiltCorr);
% %         phzWrapped = angle(IfiltCorr);
% % 
% %         [IfiltCorr, model] = insar.remove_phzramp_poisson(ccor, cor);
% %         ccor = IfiltCorr; cor = abs(IfiltCorr);
% %                 [IfiltCorr, model] = insar.remove_phzramp_poisson(ccor, cor);
% %                         ccor = IfiltCorr; cor = abs(IfiltCorr);
% %         [IfiltCorr, model] = insar.remove_phzramp_poisson(ccor, cor);
% %               ccor = IfiltCorr; cor = abs(IfiltCorr);
% %                 [IfiltCorr, model] = insar.remove_phzramp_poisson(ccor, cor);
% %                         ccor = IfiltCorr; cor = abs(IfiltCorr);
% %         [IfiltCorr, model] = insar.remove_phzramp_poisson(ccor, cor);
% % 
% %         [I_out, phi_total, logtab] = insar.remove_phzramp_multiscale(ccor, cor);
% % 
% %   ccor = I_out;
% %         cor = abs(ccor);
% %         phzWrapped = angle(ccor);
% % 
% %         [I_corr, prof, validR] = insar.remove_phzramp_rangeprofile(ccor, cor);%, orient, sigmaR, cohThresh)
% %         ccor = I_corr;
% %                 cor = abs(ccor);
% %         phzWrapped = angle(ccor);
% % 
% %         [I_corr, model] = insar.remove_phzramp_poisson_aniso(ccor, cor);
% %                 ccor = I_corr;
% %                 cor = abs(ccor);
% %         phzWrapped = angle(ccor);
% %         % trajM.pos, trajS.pos are Ns×3 ENU (lever arms already applied)
% % v   = diff(sarData.traj{1}(:,1:2),1,1);          % horizontal velocity samples
% % vm  = mean(v,1,'omitnan');                  % mean along-track vector
% % headingDeg = rad2deg(atan2(vm(2), vm(1))); % 0°=East, 90°=Nort
% % 
% %         [I_corr, phi_pred_back, thetaDeg, mdl] = ...
% %     insar.remove_phzramp_poisson_rot(ccor, cor, headingDeg);
% %         tmp = angle(I_corr);
%     end
% 
% 
% 
% 
%     % Compute unwrapped phase
%     switch lower(unwrapOpts.method)
%         case 'multiseed'
%             quality = cor;
%             % phzUnwrapped = zeros(size(cor)); % For Debugging
%             phzUnwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, quality);
%             % unwrapped = phzWrapped;
%             if unwrapOpts.postprocessRegions
%                 [phzCorrected, cycleMap, correctionMap, regionLabels] = insar.correct_phase_cycles_regionwise(phzUnwrapped, cor, unwrapOpts.tinyRegionSize);
%                 % [phzCorrected, cycleMap, correctionMap] = insar.correct_phase_cycles_regionwise(phzUnwrapped, cor, unwrapOpts.tinyRegionSize);
%                 % [phzCorrected, cycleMap] = insar.fix_integer_phase_offsets(phzUnwrapped);
%                 % unwrapped = insar.correct_nan_isolated_regions(unwrapped, unwrapOpts.filterSize);
%                 % unwrapped = insar.postprocess_unwrap_regions(unwrapped, quality,...
%                 %     unwrapOpts);
%             end
% 
%         case 'goldstein'
%             phzUnwrapped = insar.unwrap_goldstein(phzWrapped);
%         case 'none'
%             phzUnwrapped = phzWrapped;
%         otherwise
%             error('Unknown unwrap method: %s', unwrapOpts.method);
%     end
% 
%     % Reference phase
%     imgForCR = 'both'; % 'first' | 'second' | 'both'
%     switch lower(unwrapOpts.referenceMethod)
%         case 'crweightedmean'
%             refNum  = 0; refDen  = 0;     % unwrapped weighted mean
%             refNumW = 0; refDenW = 0;     % wrapped circular mean
% 
%             for c = 1:numel(crResult)
%                 % wEntry = utils.pickCRWeights(crResult(c), i, b, unwrapOpts.crSlcPolicy);
%                 switch imgForCR
%                     case 'first'
%                         wEntry = utils.pickCRWeights(crResult(c), i, b1, unwrapOpts.crSlcPolicy);
%                     case 'second'
%                         wEntry = utils.pickCRWeights(crResult(c), j, b2, unwrapOpts.crSlcPolicy);
%                     case 'both'
%                         % combine both entries (simple average here)
%                         e1 = utils.pickCRWeights(crResult(c), i, b1, unwrapOpts.crSlcPolicy);
%                         e2 = utils.pickCRWeights(crResult(c), j, b2, unwrapOpts.crSlcPolicy);
%                         wEntry = utils.combineWeightEntries(e1, e2); % see note below
%                     otherwise
%                         wEntry = utils.pickCRWeights(crResult(c), i, b1, unwrapOpts.crSlcPolicy);
%                 end
%                 if isempty(wEntry), continue; end
% 
%                 px = wEntry.indices(:);
%                 wt = wEntry.value(:);
% 
%                 % Unwrapped: weighted arithmetic mean
%                 vU = isfinite(phzUnwrapped(px));
%                 if any(vU)
%                     refNum = refNum + sum(wt(vU) .* phzUnwrapped(px(vU)));
%                     refDen = refDen + sum(wt(vU));
%                 end
% 
%                 % Wrapped: weighted circular mean
%                 vW = isfinite(phzWrapped(px));
%                 if any(vW)
%                     refNumW = refNumW + sum(wt(vW) .* exp(1i * phzWrapped(px(vW))));
%                     refDenW = refDenW + sum(wt(vW));
%                 end
%             end
% 
%             if refDen > 0
%                 referencePhase = refNum / refDen;
%             else
%                 referencePhase = mean(phzUnwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
%             end
% 
%             if refDenW > 0
%                 referencePhaseWrapped = angle(refNumW / refDenW);
%             else
%                 mask = cor >= unwrapOpts.qualityThresh;
%                 referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(mask)), 'omitnan'));
%             end
% 
%             % case 'crweightedmean'
%             %     refPhasesUnwrapped = []; refPhasesWrapped = [];
%             %     weights = [];
%             %
%             %     for c = 1:numel(crResult)
%             %         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
%             %         wEntry = crResult(c).Weights{i,j};
%             %         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
%             %         px = wEntry.indices;
%             %         wt = wEntry.value;
%             %
%             %         % Skip if unwrapped is invalid
%             %         if any(isnan(phzUnwrapped(px))), continue; end
%             %
%             %         % Weighted mean phase (unwrapped)
%             %         refPhasesUnwrapped(end+1) = sum(wt .* phzUnwrapped(px));
%             %
%             %         % Weighted mean phase (wrapped)
%             %         phzWrap = phzWrapped(px);
%             %         refPhasesWrapped(end+1) = angle(sum(wt .* exp(1i * phzWrap)));
%             %
%             %         weights(end+1) = sum(wt);
%             %     end
%             %
%             %     if ~isempty(refPhasesUnwrapped)
%             %         referencePhase        = sum(refPhasesUnwrapped .* weights) / sum(weights);
%             %         referencePhaseWrapped = angle(sum(exp(1i * refPhasesWrapped) .* weights));
%             %     else
%             %         referencePhase        = mean(phzUnwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
%             %         referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(cor >= unwrapOpts.qualityThresh)), 'omitnan'));
%             %     end
%         case 'crmean'
%             refPx = [];
%             for c = 1:numel(crResult)
%                 switch imgForCR
%                     case 'first'
%                         wEntry = utils.pickCRWeights(crResult(c), i, b1, unwrapOpts.crSlcPolicy);
%                     case 'second'
%                         wEntry = utils.pickCRWeights(crResult(c), j, b2, unwrapOpts.crSlcPolicy);
%                     case 'both'
%                         % combine both entries (simple average here)
%                         e1 = utils.pickCRWeights(crResult(c), i, b1, unwrapOpts.crSlcPolicy);
%                         e2 = utils.pickCRWeights(crResult(c), j, b2, unwrapOpts.crSlcPolicy);
%                         wEntry = utils.combineWeightEntries(e1, e2); % see note below
%                     otherwise
%                         wEntry = utils.pickCRWeights(crResult(c), i, b1, unwrapOpts.crSlcPolicy);
%                 end
%                 % wEntry = utils.pickCRWeights(crResult(c), i, b, unwrapOpts.crSlcPolicy);
%                 if isempty(wEntry) || ~isfield(wEntry,'indices'), continue; end
%                 refPx = [refPx; wEntry.indices(:)];
%             end
%             refPx = unique(refPx);
% 
%             referencePhase        = mean(phzUnwrapped(refPx), 'omitnan');
%             referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(refPx)), 'omitnan'));
% 
%             % case 'crmean'
%             %     refPx = [];
%             %     for c = 1:numel(crResult)
%             %         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
%             %         wEntry = crResult(c).Weights{i,j};
%             %         if ~isfield(wEntry, 'indices'), continue; end
%             %         refPx = [refPx; wEntry.indices(:)];
%             %     end
%             %     refPx = unique(refPx);
%             %
%             %     referencePhase        = mean(phzUnwrapped(refPx), 'omitnan');
%             %     referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(refPx)), 'omitnan'));
% 
%         case 'scenemean'
%             mask = cor >= unwrapOpts.qualityThresh;
%             referencePhase        = mean(phzUnwrapped(mask), 'omitnan');
%             referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(mask)), 'omitnan'));
% 
%         case 'none'
%             referencePhase        = 0;
%             referencePhaseWrapped = 0;
% 
%         otherwise
%             error('Unknown reference method: %s', unwrapOpts.referenceMethod);
%     end
% 
%     % Reference Wrapped Phase
%     % referencePhaseWrapped = angle(exp(1i * referencePhase));
%     phzReferenced = phzUnwrapped - referencePhase;
%     phzWrappedReferenced = angle(exp(1i * (phzWrapped - referencePhaseWrapped)));
%     % Store Results
%     k = k + 1;
%     insarData(k).pair = [i j];
%     insarData(k).burst = [b1 b2];
%     insarData(k).phzWrapped = phzWrapped;
%     insarData(k).phzWrappedReferenced = phzWrappedReferenced;
%     insarData(k).refValueWrapped = referencePhaseWrapped;
%     insarData(k).coherence = cor;
%     insarData(k).complexCoherence = ccor;
%     insarData(k).complexPhase = cphz;
%     insarData(k).phzUnwrapped = phzUnwrapped;
%     if unwrapOpts.postprocessRegions && exist("phzCorrected","var")
%         insarData(k).phzCorrected = phzCorrected;
%         insarData(k).phzCycleMap = cycleMap;
%         insarData(k).phzCorrectionMap = correctionMap;
%         insarData(k).regionLabels = regionLabels;
%         insarData(k).phzReferenced = phzReferenced;
%     else
%         insarData(k).phzReferenced = phzReferenced;
%     end
%     insarData(k).refCR = unwrapOpts.referenceMethod;
%     insarData(k).refValue = referencePhase;
%     insarData(k).unwrapMethod = unwrapOpts.method;
%     insarData(k).refMethod = unwrapOpts.referenceMethod;
% 
%     % -- Estimate signal penetration using phase and geometry (pair-specific) --
%     idxG = utils.lookup_pair_geom(geomData, i, b1, j, b2);   % exact pair (i,b1)-(j,b2)
% 
%     if ~isnan(idxG)
%         geom = geomData.sarGeometry{idxG};
% 
%         % (Optional but recommended) ensure we used the same mask for this IFG
%         % If you apply the mask earlier, you can skip this.
%         if isfield(geom,'lookMask') && ~isempty(geom.lookMask)
%             if ~isequal(size(geom.lookMask), size(phzUnwrapped))
%                 % guard for size mismatch if needed (e.g., crop/reshape)
%                 % geom.lookMask = imresize(geom.lookMask, size(phzUnwrapped), 'nearest');
%             end
%             % Do not re-mask phzUnwrapped here; assume you masked slc1/slc2 earlier.
%         end
% 
%         % Use the geometry for THIS pair
%         if isfield(geom,'slant') && isfield(geom,'incidence') && isfield(geom,'Bperp')
%             rslant = geom.slant;          % same size as the DEM grid
%             inc    = geom.incidence;      % radians (convert if stored in deg)
%             Bperp  = geom.Bperp;          % per-pixel or scalar
% 
%             % Degrees → radians guard
%             if any(inc(:) > pi, 'all'), inc = deg2rad(inc); end
% 
%             % Sensitivity check
%             target_sensitivity = 5; % m per radian
%             sin_theta = sin(inc);
%             sin_theta(~isfinite(sin_theta) | abs(sin_theta) < 1e-6) = eps;
% 
%             Bperp_thresh = (params.lambda .* rslant) ./ (4 .* pi .* target_sensitivity .* sin_theta);
% 
%             mean_thresh = median(Bperp_thresh(:), 'omitnan');
%             mean_Bperp  = median(abs(Bperp(:)),     'omitnan');
% 
%             if mean_Bperp < mean_thresh
%                 warning('Mean Bperp %.2f m < min required %.2f m (pair %d-%d bursts [%d %d], geom %d)', ...
%                     mean_Bperp, mean_thresh, i, j, b1, b2, idxG);
%                 insarData(k).penetration      = NaN(size(phzUnwrapped));
%                 insarData(k).penetrationValid = false;
%             else
%                 % Use referenced, unwrapped phase for penetration
%                 phi = phzUnwrapped - referencePhase;
% 
%                 % Your existing penetration conversion
%                 penetration = insar.phase_to_penetration(phi, params.lambda, Bperp, rslant, inc);
% 
%                 insarData(k).penetration      = penetration;
%                 insarData(k).penetrationValid = true;
%             end
% 
%             % Bookkeeping
%             insarData(k).meanBperp         = mean_Bperp;
%             insarData(k).minRequiredBperp  = mean_thresh;
%             insarData(k).sensitivityTarget = target_sensitivity;
%             insarData(k).geomIdx           = idxG;    % <-- remember which geometry was used
%         else
%             warning('Missing geometry fields for geom idx %d. Skipping penetration.', idxG);
%         end
%     else
%         warning('No geometry for pair (%d,b%02d)×(%d,b%02d). Skipping penetration.', i,b1,j,b2);
%     end
%     % -- Estimate signal penetration using phase and geometry --
%     % if isfield(geomData,"slcGeomIndex") && isfield(geomData,"sarGeometry")
%     %     idx1 = geomData.slcGeomIndex(i, b1).idx;
%     %     % For Bperp/inc/rslant you likely want the geometry used for the interferogram.
%     %     % If your pipeline stores pair-geometry only when idx1==idx2, guard it:
%     %     if ~isnan(idx1)
%     %         geom = geomData.sarGeometry{idx1};
%     %         if isfield(geom,'slant') && isfield(geom,'incidence') && isfield(geom,'Bperp')
%     %             rslant = geom.slant;
%     %             inc    = geom.incidence;
%     %             Bperp  = geom.Bperp;
%     % 
%     %             target_sensitivity = 2; % m/rad
%     %             if any(inc(:) > pi), inc = deg2rad(inc); end
%     %             sin_theta = sin(inc);
%     %             sin_theta(~isfinite(sin_theta) | abs(sin_theta) < 1e-6) = eps;
%     % 
%     %             Bperp_thresh = (params.lambda.*rslant) ./ (4.*pi.*target_sensitivity.*sin_theta);
%     % 
%     %             mean_thresh = median(Bperp_thresh(:), 'omitnan');
%     %             mean_Bperp  = median(abs(Bperp(:)), 'omitnan');
%     % 
%     %             if mean_Bperp < mean_thresh
%     %                 warning('Mean Bperp %.2f m < min required %.2f m (pair %d-%d bursts [%d %d])', ...
%     %                     mean_Bperp, mean_thresh, i, j, b1, b2);
%     %                 insarData(k).penetration = NaN(size(phzUnwrapped));
%     %                 insarData(k).penetrationValid = false;
%     %             else
%     %                 penetration = insar.phase_to_penetration(phzUnwrapped - referencePhase, params.lambda, Bperp, rslant, inc);
%     %                 insarData(k).penetration = penetration;
%     %                 insarData(k).penetrationValid = true;
%     %             end
%     %             insarData(k).meanBperp = mean_Bperp;
%     %             insarData(k).minRequiredBperp = mean_thresh;
%     %             insarData(k).sensitivityTarget = target_sensitivity;
%     %         else
%     %             warning('Missing geometry fields for (dir %d, burst %d). Skipping penetration.', i, b1);
%     %         end
%     %     else
%     %         warning('No geometry index for (dir %d, burst %d). Skipping penetration.', i, b1);
%     %     end
%     % end
% 
%     fprintf('Processed (dir %d,b%02d) × (dir %d,b%02d). Elapsed: %.2f s\n', ...
%         i, b1, j, b2, toc);
% 
%     % if isfield(geomData, "slcGeomIndex") && isfield(geomData, "sarGeometry")
%     %     idx = geomData.slcGeomIndex(i, b).idx;
%     %     geom = geomData.sarGeometry{idx};
%     % 
%     %     rslant = geom.slant;
%     %     inc = geom.incidence;
%     %     Bperp = geom.Bperp;
%     % 
%     %     target_sensitivity = 2; % m per radian
%     % 
%     %     % Convert incidence to radians if needed
%     %     if any(inc(:) > pi), inc = deg2rad(inc); end
%     %     % Guard against zero or invalid sin(θ)
%     %     sin_theta = sin(inc);
%     %     sin_theta(~isfinite(sin_theta) | abs(sin_theta) < 1e-6) = eps;
%     % 
%     %     % Avoid division by zero in sensitivity
%     %     if target_sensitivity <= 0
%     %         error('Invalid target_sensitivity: must be > 0');
%     %     end
%     % 
%     %     % Compute threshold Bperp per pixel
%     %     Bperp_thresh = (params.lambda.*rslant) ./ (4.*pi.* target_sensitivity.*sin_theta);
%     % 
%     %     % Use average Bperp across scene
%     %     mean_thresh = median(Bperp_thresh(:), 'omitnan');
%     %     mean_Bperp  = median(abs(Bperp(:)), 'omitnan');
%     % 
%     %     if mean_Bperp < mean_thresh
%     %         warning('Mean Bperp %.2f m < min required %.2f m (pair %d-%d burst %d)', ...
%     %             mean_Bperp, mean_thresh, i, j, b);
%     %         insarData(k).penetration = NaN(size(phzUnwrapped));
%     %         insarData(k).penetrationValid = false;
%     %     else
%     %         penetration = phase_to_penetration(phzUnwrapped - referencePhase, params.lambda, Bperp, rslant, inc);
%     %         insarData(k).penetration = penetration;
%     %         insarData(k).penetrationValid = true;
%     %     end
%     %     insarData(k).meanBperp = mean_Bperp;
%     %     insarData(k).minRequiredBperp = mean_thresh;
%     %     insarData(k).sensitivityTarget = target_sensitivity;
%     % end
%     % fprintf('Processed (pair %d-%d burst %d). Elapsed time: %.2f seconds\n', ...
%     %     i, j, b, toc');
% end
% end
% % function insarData = process_interferometric_phase(sarData, geomData, crResult, unwrapOpts)
% % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % %   Inputs:
% % %       sarData     - Struct of SLC data and metadata
% % %       geomData    - Struct of geometry data including lookMask
% % %       crResult    - CR calibration structure (from CRpower)
% % %       unwrapOpts  - Struct with fields:
% % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % %           .filterSize      : Interferogram filter window (default = 5)
% % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % %           .alpha           : Coherence weighting exponent (default = 2)
% % %           .qualityThresh   : Mask threshold (default = 0.3)
% %
% % if nargin < 4, unwrapOpts = struct(); end
% % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
% %
% % numDirs = numel(sarData);
% %
% % insarData = struct([]);
% % k = 0;
% %
% % for i = 1:numDirs
% %     for j = i+1:numDirs
% %         % Ensure bursts match
% %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% %             continue;
% %         end
% %
% %         for b = 1:numel(sarData(i).slc)
% %             slc1 = sarData(i).slc{b};
% %             slc2 = sarData(j).slc{b};
% %
% %             % Apply lookMask using geomData
% %             if isfield(geomData, "slcGeomIndex") && isfield(geomData, "sarGeometry")
% %                 idx1 = geomData.slcGeomIndex(i, b).idx;
% %                 idx2 = geomData.slcGeomIndex(j, b).idx;
% %                 if ~isnan(idx1) && ~isnan(idx2) && idx1 == idx2
% %                     mask1 = geomData.sarGeometry{idx1}.lookMask;
% %                     mask = mask1;
% %                     slc1(~mask) = NaN;
% %                     slc2(~mask) = NaN;
% %                 end
% %             end
% %
% %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% %                 unwrapOpts.sigma, unwrapOpts.alpha);
% %
% %             % Inline unwrapping logic here
% %             switch lower(unwrapOpts.method)
% %                 case 'multiseed'
% %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% %                 case 'goldstein'
% %                     unwrapped = insar.unwrap_goldstein(phzWrapped, unwrapOpts.alpha, unwrapOpts.filterSize);
% %                 otherwise
% %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% %             end
% %
% %             % Inline reference logic
% %             switch lower(unwrapOpts.referenceMethod)
% %                 case 'crweightedmean'
% %                     refPhases = [];
% %                     weights = [];
% %                     for c = 1:numel(crResult)
% %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% %                         wEntry = crResult(c).Weights{i,j};
% %                         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% %                         px = wEntry.indices;
% %                         wt = wEntry.value;
% %                         if any(isnan(unwrapped(px))), continue; end
% %                         refPhases(end+1) = sum(wt .* unwrapped(px));
% %                         weights(end+1) = sum(wt);
% %                     end
% %                     if ~isempty(refPhases)
% %                         referencePhase = sum(refPhases .* weights) / sum(weights);
% %                     else
% %                         referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% %                     end
% %
% %                 case 'crmean'
% %                     refPx = [];
% %                     for c = 1:numel(crResult)
% %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% %                         wEntry = crResult(c).Weights{i,j};
% %                         if ~isfield(wEntry, 'indices'), continue; end
% %                         refPx = [refPx; wEntry.indices(:)];
% %                     end
% %                     refPx = unique(refPx);
% %                     referencePhase = mean(unwrapped(refPx), 'omitnan');
% %
% %                 case 'scenemean'
% %                     referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% %
% %                 case 'none'
% %                     referencePhase = 0;
% %
% %                 otherwise
% %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% %             end
% %
% %             k = k + 1;
% %             insarData(k).pair = [i j];
% %             insarData(k).burst = b;
% %             insarData(k).phzWrapped = phzWrapped;
% %             insarData(k).coherence = cor;
% %             insarData(k).phzUnwrapped = unwrapped;
% %             insarData(k).phzReferenced = unwrapped - referencePhase;
% %             insarData(k).refCR = unwrapOpts.referenceMethod;
% %             insarData(k).refValue = referencePhase;
% %             insarData(k).unwrapMethod = unwrapOpts.method;
% %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% %         end
% %     end
% % end
% %
% % end
% %
% % % function insarData = process_interferometric_phase(sarData, geomData, crResult, unwrapOpts)
% % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % %   Inputs:
% % % %       sarData     - Struct of SLC data and metadata
% % % %       geomData    - Struct of geometry data including lookMask
% % % %       crResult    - CR calibration structure (from CRpower)
% % % %       unwrapOpts  - Struct with fields:
% % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % %           .filterSize      : Interferogram filter window (default = 5)
% % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % %           .qualityThresh   : Mask threshold (default = 0.75)
% % %
% % % if nargin < 4, unwrapOpts = struct(); end
% % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.75; end
% % %
% % % numDirs = numel(sarData);
% % %
% % % insarData = struct([]);
% % % k = 0;
% % %
% % % for i = 1:numDirs
% % %     for j = i+1:numDirs
% % %         % Ensure bursts match
% % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % %             continue;
% % %         end
% % %
% % %         for b = 1:numel(sarData(i).slc)
% % %             slc1 = sarData(i).slc{b};
% % %             slc2 = sarData(j).slc{b};
% % %             mask = []; % default
% % %             if isfield(geomData(i).lookMask, 'data')
% % %                 mask = geomData(i).lookMask(b).data & geomData(j).lookMask(b).data;
% % %                 slc1(~mask) = NaN;
% % %                 slc2(~mask) = NaN;
% % %             end
% % %
% % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % %
% % %             % Inline unwrapping logic here
% % %             switch lower(unwrapOpts.method)
% % %                 case 'multiseed'
% % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % %                 case 'goldstein'
% % %                     unwrapped = insar.unwrap_goldstein(phzWrapped);
% % %                 otherwise
% % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % %             end
% % %
% % %             % Inline reference logic
% % %             switch lower(unwrapOpts.referenceMethod)
% % %                 case 'crweightedmean'
% % %                     refPhases = [];
% % %                     weights = [];
% % %                     for c = 1:numel(crResult)
% % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % %                         wEntry = crResult(c).Weights{i,j};
% % %                         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% % %                         px = wEntry.indices;
% % %                         wt = wEntry.value;
% % %                         if any(isnan(unwrapped(px))), continue; end
% % %                         refPhases(end+1) = sum(wt .* unwrapped(px));
% % %                         weights(end+1) = sum(wt);
% % %                     end
% % %                     if ~isempty(refPhases)
% % %                         referencePhase = sum(refPhases .* weights) / sum(weights);
% % %                     else
% % %                         referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % %                     end
% % %
% % %                 case 'crmean'
% % %                     refPx = [];
% % %                     for c = 1:numel(crResult)
% % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % %                         wEntry = crResult(c).Weights{i,j};
% % %                         if ~isfield(wEntry, 'indices'), continue; end
% % %                         refPx = [refPx; wEntry.indices(:)];
% % %                     end
% % %                     refPx = unique(refPx);
% % %                     referencePhase = mean(unwrapped(refPx), 'omitnan');
% % %
% % %                 case 'scenemean'
% % %                     referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % %
% % %                 case 'none'
% % %                     referencePhase = 0;
% % %
% % %                 otherwise
% % %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% % %             end
% % %
% % %             k = k + 1;
% % %             insarData(k).pair = [i j];
% % %             insarData(k).burst = b;
% % %             insarData(k).phzWrapped = phzWrapped;
% % %             insarData(k).coherence = cor;
% % %             insarData(k).phzUnwrapped = unwrapped;
% % %             insarData(k).phzReferenced = unwrapped - referencePhase;
% % %             insarData(k).refCR = unwrapOpts.referenceMethod;
% % %             insarData(k).refValue = referencePhase;
% % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % %         end
% % %     end
% % % end
% % %
% % % end
% % %
% % % % function insarData = process_interferometric_phase(sarData, crResult, unwrapOpts)
% % % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % % %   Inputs:
% % % % %       sarData   - Struct of SLC data and metadata
% % % % %       crResult  - CR calibration structure (from CRpower)
% % % % %       unwrapOpts - Struct with fields:
% % % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % % %           .filterSize      : Interferogram filter window (default = 5)
% % % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % % %           .qualityThresh   : Mask threshold (default = 0.75)
% % % %
% % % % if nargin < 3, unwrapOpts = struct(); end
% % % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
% % % %
% % % % numDirs = numel(sarData);
% % % %
% % % % insarData = struct([]);
% % % % k = 0;
% % % %
% % % % for i = 1:numDirs
% % % %     for j = i+1:numDirs
% % % %         % Ensure bursts match
% % % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % % %             continue;
% % % %         end
% % % %
% % % %         for b = 1:numel(sarData(i).slc)
% % % %             slc1 = sarData(i).slc{b};
% % % %             slc2 = sarData(j).slc{b};
% % % %
% % % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % % %
% % % %             % Inline unwrapping logic here
% % % %             switch lower(unwrapOpts.method)
% % % %                 case 'multiseed'
% % % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % % %                 case 'goldstein'
% % % %                     unwrapped = insar.unwrap_goldstein(phzWrapped);
% % % %                 otherwise
% % % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % % %             end
% % % %
% % % %             % Inline reference logic
% % % %             switch lower(unwrapOpts.referenceMethod)
% % % %                 case 'crweightedmean'
% % % %                     refPhases = [];
% % % %                     weights = [];
% % % %                     for c = 1:numel(crResult)
% % % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % % %                         wEntry = crResult(c).Weights{i,j};
% % % %                         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% % % %                         px = wEntry.indices;
% % % %                         wt = wEntry.value;
% % % %                         if any(isnan(unwrapped(px))), continue; end
% % % %                         refPhases(end+1) = sum(wt .* unwrapped(px));
% % % %                         weights(end+1) = sum(wt);
% % % %                     end
% % % %                     if ~isempty(refPhases)
% % % %                         referencePhase = sum(refPhases .* weights) / sum(weights);
% % % %                     else
% % % %                         referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % % %                     end
% % % %
% % % %                 case 'crmean'
% % % %                     refPx = [];
% % % %                     for c = 1:numel(crResult)
% % % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % % %                         wEntry = crResult(c).Weights{i,j};
% % % %                         if ~isfield(wEntry, 'indices'), continue; end
% % % %                         refPx = [refPx; wEntry.indices(:)];
% % % %                     end
% % % %                     refPx = unique(refPx);
% % % %                     referencePhase = mean(unwrapped(refPx), 'omitnan');
% % % %
% % % %                 case 'scenemean'
% % % %                     referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % % %
% % % %                 case 'none'
% % % %                     referencePhase = 0;
% % % %
% % % %                 otherwise
% % % %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% % % %             end
% % % %
% % % %             k = k + 1;
% % % %             insarData(k).pair = [i j];
% % % %             insarData(k).burst = b;
% % % %             insarData(k).phzWrapped = phzWrapped;
% % % %             insarData(k).coherence = cor;
% % % %             insarData(k).phzUnwrapped = unwrapped;
% % % %             insarData(k).phzReferenced = unwrapped - referencePhase;
% % % %             insarData(k).refCR = unwrapOpts.referenceMethod;
% % % %             insarData(k).refValue = referencePhase;
% % % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % % %         end
% % % %     end
% % % % end
% % % %
% % % % end
% % % %
% % % % % function insarData = process_interferometric_phase(sarData, crResult, unwrapOpts)
% % % % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % % % %   Inputs:
% % % % % %       sarData   - Struct of SLC data and metadata
% % % % % %       crResult  - CR calibration structure (from CRpower)
% % % % % %       unwrapOpts - Struct with fields:
% % % % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % % % %           .filterSize      : Interferogram filter window (default = 5)
% % % % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % % % %           .qualityThresh   : Mask threshold (default = 0.75)
% % % % %
% % % % % if nargin < 3, unwrapOpts = struct(); end
% % % % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.75; end
% % % % %
% % % % % numDirs = numel(sarData);
% % % % %
% % % % % insarData = struct([]);
% % % % % k = 0;
% % % % %
% % % % % for i = 1:numDirs
% % % % %     for j = i+1:numDirs
% % % % %         % Ensure bursts match
% % % % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % % % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % % % %             continue;
% % % % %         end
% % % % %
% % % % %         for b = 1:numel(sarData(i).slc)
% % % % %             slc1 = sarData(i).slc{b};
% % % % %             slc2 = sarData(j).slc{b};
% % % % %
% % % % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % % % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % % % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % % % %
% % % % %             % Unwrapping
% % % % %             switch lower(unwrapOpts.method)
% % % % %                 case 'multiseed'
% % % % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % % % %                 case 'goldstein'
% % % % %                     unwrapped = insar.unwrap_goldstein(phzWrapped);
% % % % %                 otherwise
% % % % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % % % %             end
% % % % %             % Phase Referencing
% % % % %             [refPhz, refCR, refVal] = insar.reference_phase(unwrapped, crResult, unwrapOpts.referenceMethod);
% % % % %
% % % % %             k = k + 1;
% % % % %             insarData(k).pair = [i j];
% % % % %             insarData(k).burst = b;
% % % % %             insarData(k).phzWrapped = phzWrapped;
% % % % %             insarData(k).coherence = cor;
% % % % %             insarData(k).phzUnwrapped = unwrapped;
% % % % %             insarData(k).phzReferenced = refPhz;
% % % % %             insarData(k).refCR = refCR;
% % % % %             insarData(k).refValue = refVal;
% % % % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % % % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % % % %         end
% % % % %     end
% % % % % end
% % % % %
% % % % % end
% 
% % function [insarData, unwrapOpts] = process_interferometric_phase(sarData, geomData, crResult, params, unwrapOpts)
% % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % %   Includes multiseed unwrapping with optional post-processing upgrades:
% % %       - Region bridging (cross-region stitching)
% % %       - Phase gradient blending
% % %       - MST-based merging
% % %       - Fallback for tiny regions
% %
% % if nargin < 5, unwrapOpts = struct(); end
% % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % if ~isfield(unwrapOpts, 'useCoherenceFilter'), unwrapOpts.useCoherenceFilter = true; end
% % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize ./ 2; end
% % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
% % if ~isfield(unwrapOpts, 'postprocessRegions'), unwrapOpts.postprocessRegions = false; end
% % if ~isfield(unwrapOpts, 'stitchMethod'), unwrapOpts.stitchMethod = 'mst'; end
% % if ~isfield(unwrapOpts, 'tinyRegionSize'), unwrapOpts.tinyRegionSize = 30; end
% % if ~isfield(unwrapOpts, 'stitchAlpha'), unwrapOpts.stitchAlpha = 0.5; end
% % if ~isfield(unwrapOpts,'crSlcPolicy'), unwrapOpts.crSlcPolicy = 'burst'; end
% %
% %
% % numDirs = numel(sarData);
% % insarData = struct([]);
% % k = 0;
% %
% % for i = 1:numDirs
% %     for j = i+1:numDirs
% %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% %             continue;
% %         end
% %
% %         for b = 1:numel(sarData(i).slc)
% %             tic
% %             slc1 = sarData(i).slc{b};
% %             slc2 = sarData(j).slc{b};
% %
% %             % Apply lookMask using geomData
% %             if isfield(geomData, "slcGeomIndex") && isfield(geomData, "sarGeometry")
% %                 idx1 = geomData.slcGeomIndex(i, b).idx;
% %                 idx2 = geomData.slcGeomIndex(j, b).idx;
% %                 if ~isnan(idx1) && ~isnan(idx2) && idx1 == idx2
% %                     mask = geomData.sarGeometry{idx1}.lookMask;
% %                     slc1(~mask) = NaN;
% %                     slc2(~mask) = NaN;
% %                 end
% %             end
% %
% %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, unwrapOpts.useCoherenceFilter, ...
% %                 unwrapOpts.sigma, unwrapOpts.alpha);
% %
% %             % Compute unwrapped phase
% %             switch lower(unwrapOpts.method)
% %                 case 'multiseed'
% %                     quality = cor;
% %                     phzUnwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, quality, unwrapOpts.qualityThresh);
% %                     % unwrapped = phzWrapped;
% %                     if unwrapOpts.postprocessRegions
% %                         [phzCorrected, cycleMap, correctionMap, regionLabels] = insar.correct_phase_cycles_regionwise(phzUnwrapped, cor, unwrapOpts.tinyRegionSize);
% %                         % [phzCorrected, cycleMap, correctionMap] = insar.correct_phase_cycles_regionwise(phzUnwrapped, cor, unwrapOpts.tinyRegionSize);
% %                         % [phzCorrected, cycleMap] = insar.fix_integer_phase_offsets(phzUnwrapped);
% %                         % unwrapped = insar.correct_nan_isolated_regions(unwrapped, unwrapOpts.filterSize);
% %                         % unwrapped = insar.postprocess_unwrap_regions(unwrapped, quality,...
% %                         %     unwrapOpts);
% %                     end
% %
% %                 case 'goldstein'
% %                     phzUnwrapped = insar.unwrap_goldstein(phzWrapped);
% %                 otherwise
% %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% %             end
% %
% %             % Reference phase
% %             switch lower(unwrapOpts.referenceMethod)
% %                 case 'crweightedmean'
% %                     refNum  = 0; refDen  = 0;     % unwrapped weighted mean
% %                     refNumW = 0; refDenW = 0;     % wrapped circular mean
% %
% %                     for c = 1:numel(crResult)
% %                         wEntry = utils.pickCRWeights(crResult(c), i, b, unwrapOpts.crSlcPolicy);
% %                         if isempty(wEntry), continue; end
% %
% %                         px = wEntry.indices(:);
% %                         wt = wEntry.value(:);
% %
% %                         % Unwrapped: weighted arithmetic mean
% %                         vU = isfinite(phzUnwrapped(px));
% %                         if any(vU)
% %                             refNum = refNum + sum(wt(vU) .* phzUnwrapped(px(vU)));
% %                             refDen = refDen + sum(wt(vU));
% %                         end
% %
% %                         % Wrapped: weighted circular mean
% %                         vW = isfinite(phzWrapped(px));
% %                         if any(vW)
% %                             refNumW = refNumW + sum(wt(vW) .* exp(1i * phzWrapped(px(vW))));
% %                             refDenW = refDenW + sum(wt(vW));
% %                         end
% %                     end
% %
% %                     if refDen > 0
% %                         referencePhase = refNum / refDen;
% %                     else
% %                         referencePhase = mean(phzUnwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% %                     end
% %
% %                     if refDenW > 0
% %                         referencePhaseWrapped = angle(refNumW / refDenW);
% %                     else
% %                         mask = cor >= unwrapOpts.qualityThresh;
% %                         referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(mask)), 'omitnan'));
% %                     end
% %
% %                     % case 'crweightedmean'
% %                     %     refPhasesUnwrapped = []; refPhasesWrapped = [];
% %                     %     weights = [];
% %                     %
% %                     %     for c = 1:numel(crResult)
% %                     %         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% %                     %         wEntry = crResult(c).Weights{i,j};
% %                     %         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% %                     %         px = wEntry.indices;
% %                     %         wt = wEntry.value;
% %                     %
% %                     %         % Skip if unwrapped is invalid
% %                     %         if any(isnan(phzUnwrapped(px))), continue; end
% %                     %
% %                     %         % Weighted mean phase (unwrapped)
% %                     %         refPhasesUnwrapped(end+1) = sum(wt .* phzUnwrapped(px));
% %                     %
% %                     %         % Weighted mean phase (wrapped)
% %                     %         phzWrap = phzWrapped(px);
% %                     %         refPhasesWrapped(end+1) = angle(sum(wt .* exp(1i * phzWrap)));
% %                     %
% %                     %         weights(end+1) = sum(wt);
% %                     %     end
% %                     %
% %                     %     if ~isempty(refPhasesUnwrapped)
% %                     %         referencePhase        = sum(refPhasesUnwrapped .* weights) / sum(weights);
% %                     %         referencePhaseWrapped = angle(sum(exp(1i * refPhasesWrapped) .* weights));
% %                     %     else
% %                     %         referencePhase        = mean(phzUnwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% %                     %         referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(cor >= unwrapOpts.qualityThresh)), 'omitnan'));
% %                     %     end
% %                 case 'crmean'
% %                     refPx = [];
% %                     for c = 1:numel(crResult)
% %                         wEntry = utils.pickCRWeights(crResult(c), i, b, unwrapOpts.crSlcPolicy);
% %                         if isempty(wEntry) || ~isfield(wEntry,'indices'), continue; end
% %                         refPx = [refPx; wEntry.indices(:)];
% %                     end
% %                     refPx = unique(refPx);
% %
% %                     referencePhase        = mean(phzUnwrapped(refPx), 'omitnan');
% %                     referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(refPx)), 'omitnan'));
% %
% %                     % case 'crmean'
% %                     %     refPx = [];
% %                     %     for c = 1:numel(crResult)
% %                 %         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% %                 %         wEntry = crResult(c).Weights{i,j};
% %                 %         if ~isfield(wEntry, 'indices'), continue; end
% %                 %         refPx = [refPx; wEntry.indices(:)];
% %                 %     end
% %                 %     refPx = unique(refPx);
% %                 %
% %                 %     referencePhase        = mean(phzUnwrapped(refPx), 'omitnan');
% %                 %     referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(refPx)), 'omitnan'));
% %
% %                 case 'scenemean'
% %                     mask = cor >= unwrapOpts.qualityThresh;
% %                     referencePhase        = mean(phzUnwrapped(mask), 'omitnan');
% %                     referencePhaseWrapped = angle(mean(exp(1i * phzWrapped(mask)), 'omitnan'));
% %
% %                 case 'none'
% %                     referencePhase        = 0;
% %                     referencePhaseWrapped = 0;
% %
% %                 otherwise
% %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% %             end
% %
% %             % Reference Wrapped Phase
% %             % referencePhaseWrapped = angle(exp(1i * referencePhase));
% %             phzReferenced = phzUnwrapped - referencePhase;
% %             phzWrappedReferenced = angle(exp(1i * (phzWrapped - referencePhaseWrapped)));
% %
% %
% %             k = k + 1;
% %             insarData(k).pair = [i j];
% %             insarData(k).burst = b;
% %             insarData(k).phzWrapped = phzWrapped;
% %             insarData(k).phzWrappedReferenced = phzWrappedReferenced;
% %             insarData(k).refValueWrapped = referencePhaseWrapped;
% %             insarData(k).coherence = cor;
% %             insarData(k).phzUnwrapped = phzUnwrapped;
% %             if unwrapOpts.postprocessRegions && exist("phzCorrected","var")
% %             insarData(k).phzCorrected = phzCorrected;
% %             insarData(k).phzCycleMap = cycleMap;
% %             insarData(k).phzCorrectionMap = correctionMap;
% %             insarData(k).regionLabels = regionLabels;
% %             insarData(k).phzReferenced = phzReferenced;
% %             else
% %             insarData(k).phzReferenced = phzReferenced;
% %             end
% %             insarData(k).refCR = unwrapOpts.referenceMethod;
% %             insarData(k).refValue = referencePhase;
% %             insarData(k).unwrapMethod = unwrapOpts.method;
% %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% %
% %             % -- Estimate signal penetration using phase and geometry --
% %             if isfield(geomData, "slcGeomIndex") && isfield(geomData, "sarGeometry")
% %                 idx = geomData.slcGeomIndex(i, b).idx;
% %                 geom = geomData.sarGeometry{idx};
% %
% %                 rslant = geom.slant;
% %                 inc = geom.incidence;
% %                 Bperp = geom.Bperp;
% %
% %                 target_sensitivity = 2; % m per radian
% %
% %                 % Convert incidence to radians if needed
% %                 if any(inc(:) > pi), inc = deg2rad(inc); end
% %                 % Guard against zero or invalid sin(θ)
% %                 sin_theta = sin(inc);
% %                 sin_theta(~isfinite(sin_theta) | abs(sin_theta) < 1e-6) = eps;
% %
% %                 % Avoid division by zero in sensitivity
% %                 if target_sensitivity <= 0
% %                     error('Invalid target_sensitivity: must be > 0');
% %                 end
% %
% %                 % Compute threshold Bperp per pixel
% %                 Bperp_thresh = (params.lambda.*rslant) ./ (4.*pi.* target_sensitivity.*sin_theta);
% %
% %                 % Use average Bperp across scene
% %                 mean_thresh = median(Bperp_thresh(:), 'omitnan');
% %                 mean_Bperp  = median(abs(Bperp(:)), 'omitnan');
% %
% %                 if mean_Bperp < mean_thresh
% %                     warning('Mean Bperp %.2f m < min required %.2f m (pair %d-%d burst %d)', ...
% %                         mean_Bperp, mean_thresh, i, j, b);
% %                     insarData(k).penetration = NaN(size(phzUnwrapped));
% %                     insarData(k).penetrationValid = false;
% %                 else
% %                     penetration = phase_to_penetration(phzUnwrapped - referencePhase, params.lambda, Bperp, rslant, inc);
% %                     insarData(k).penetration = penetration;
% %                     insarData(k).penetrationValid = true;
% %                 end
% %                 insarData(k).meanBperp = mean_Bperp;
% %                 insarData(k).minRequiredBperp = mean_thresh;
% %                 insarData(k).sensitivityTarget = target_sensitivity;
% %             end
% %             fprintf('Processed (pair %d-%d burst %d). Elapsed time: %.2f seconds\n', ...
% %                 i, j, b, toc');
% %         end
% %     end
% % end
% % end
% %
% % % function insarData = process_interferometric_phase(sarData, geomData, crResult, unwrapOpts)
% % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % %   Inputs:
% % % %       sarData     - Struct of SLC data and metadata
% % % %       geomData    - Struct of geometry data including lookMask
% % % %       crResult    - CR calibration structure (from CRpower)
% % % %       unwrapOpts  - Struct with fields:
% % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % %           .filterSize      : Interferogram filter window (default = 5)
% % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % %           .qualityThresh   : Mask threshold (default = 0.3)
% % %
% % % if nargin < 4, unwrapOpts = struct(); end
% % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
% % %
% % % numDirs = numel(sarData);
% % %
% % % insarData = struct([]);
% % % k = 0;
% % %
% % % for i = 1:numDirs
% % %     for j = i+1:numDirs
% % %         % Ensure bursts match
% % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % %             continue;
% % %         end
% % %
% % %         for b = 1:numel(sarData(i).slc)
% % %             slc1 = sarData(i).slc{b};
% % %             slc2 = sarData(j).slc{b};
% % %
% % %             % Apply lookMask using geomData
% % %             if isfield(geomData, "slcGeomIndex") && isfield(geomData, "sarGeometry")
% % %                 idx1 = geomData.slcGeomIndex(i, b).idx;
% % %                 idx2 = geomData.slcGeomIndex(j, b).idx;
% % %                 if ~isnan(idx1) && ~isnan(idx2) && idx1 == idx2
% % %                     mask1 = geomData.sarGeometry{idx1}.lookMask;
% % %                     mask = mask1;
% % %                     slc1(~mask) = NaN;
% % %                     slc2(~mask) = NaN;
% % %                 end
% % %             end
% % %
% % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % %
% % %             % Inline unwrapping logic here
% % %             switch lower(unwrapOpts.method)
% % %                 case 'multiseed'
% % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % %                 case 'goldstein'
% % %                     unwrapped = insar.unwrap_goldstein(phzWrapped, unwrapOpts.alpha, unwrapOpts.filterSize);
% % %                 otherwise
% % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % %             end
% % %
% % %             % Inline reference logic
% % %             switch lower(unwrapOpts.referenceMethod)
% % %                 case 'crweightedmean'
% % %                     refPhases = [];
% % %                     weights = [];
% % %                     for c = 1:numel(crResult)
% % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % %                         wEntry = crResult(c).Weights{i,j};
% % %                         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% % %                         px = wEntry.indices;
% % %                         wt = wEntry.value;
% % %                         if any(isnan(unwrapped(px))), continue; end
% % %                         refPhases(end+1) = sum(wt .* unwrapped(px));
% % %                         weights(end+1) = sum(wt);
% % %                     end
% % %                     if ~isempty(refPhases)
% % %                         referencePhase = sum(refPhases .* weights) / sum(weights);
% % %                     else
% % %                         referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % %                     end
% % %
% % %                 case 'crmean'
% % %                     refPx = [];
% % %                     for c = 1:numel(crResult)
% % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % %                         wEntry = crResult(c).Weights{i,j};
% % %                         if ~isfield(wEntry, 'indices'), continue; end
% % %                         refPx = [refPx; wEntry.indices(:)];
% % %                     end
% % %                     refPx = unique(refPx);
% % %                     referencePhase = mean(unwrapped(refPx), 'omitnan');
% % %
% % %                 case 'scenemean'
% % %                     referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % %
% % %                 case 'none'
% % %                     referencePhase = 0;
% % %
% % %                 otherwise
% % %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% % %             end
% % %
% % %             k = k + 1;
% % %             insarData(k).pair = [i j];
% % %             insarData(k).burst = b;
% % %             insarData(k).phzWrapped = phzWrapped;
% % %             insarData(k).coherence = cor;
% % %             insarData(k).phzUnwrapped = unwrapped;
% % %             insarData(k).phzReferenced = unwrapped - referencePhase;
% % %             insarData(k).refCR = unwrapOpts.referenceMethod;
% % %             insarData(k).refValue = referencePhase;
% % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % %         end
% % %     end
% % % end
% % %
% % % end
% % %
% % % % function insarData = process_interferometric_phase(sarData, geomData, crResult, unwrapOpts)
% % % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % % %   Inputs:
% % % % %       sarData     - Struct of SLC data and metadata
% % % % %       geomData    - Struct of geometry data including lookMask
% % % % %       crResult    - CR calibration structure (from CRpower)
% % % % %       unwrapOpts  - Struct with fields:
% % % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % % %           .filterSize      : Interferogram filter window (default = 5)
% % % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % % %           .qualityThresh   : Mask threshold (default = 0.75)
% % % %
% % % % if nargin < 4, unwrapOpts = struct(); end
% % % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.75; end
% % % %
% % % % numDirs = numel(sarData);
% % % %
% % % % insarData = struct([]);
% % % % k = 0;
% % % %
% % % % for i = 1:numDirs
% % % %     for j = i+1:numDirs
% % % %         % Ensure bursts match
% % % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % % %             continue;
% % % %         end
% % % %
% % % %         for b = 1:numel(sarData(i).slc)
% % % %             slc1 = sarData(i).slc{b};
% % % %             slc2 = sarData(j).slc{b};
% % % %             mask = []; % default
% % % %             if isfield(geomData(i).lookMask, 'data')
% % % %                 mask = geomData(i).lookMask(b).data & geomData(j).lookMask(b).data;
% % % %                 slc1(~mask) = NaN;
% % % %                 slc2(~mask) = NaN;
% % % %             end
% % % %
% % % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % % %
% % % %             % Inline unwrapping logic here
% % % %             switch lower(unwrapOpts.method)
% % % %                 case 'multiseed'
% % % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % % %                 case 'goldstein'
% % % %                     unwrapped = insar.unwrap_goldstein(phzWrapped);
% % % %                 otherwise
% % % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % % %             end
% % % %
% % % %             % Inline reference logic
% % % %             switch lower(unwrapOpts.referenceMethod)
% % % %                 case 'crweightedmean'
% % % %                     refPhases = [];
% % % %                     weights = [];
% % % %                     for c = 1:numel(crResult)
% % % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % % %                         wEntry = crResult(c).Weights{i,j};
% % % %                         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% % % %                         px = wEntry.indices;
% % % %                         wt = wEntry.value;
% % % %                         if any(isnan(unwrapped(px))), continue; end
% % % %                         refPhases(end+1) = sum(wt .* unwrapped(px));
% % % %                         weights(end+1) = sum(wt);
% % % %                     end
% % % %                     if ~isempty(refPhases)
% % % %                         referencePhase = sum(refPhases .* weights) / sum(weights);
% % % %                     else
% % % %                         referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % % %                     end
% % % %
% % % %                 case 'crmean'
% % % %                     refPx = [];
% % % %                     for c = 1:numel(crResult)
% % % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % % %                         wEntry = crResult(c).Weights{i,j};
% % % %                         if ~isfield(wEntry, 'indices'), continue; end
% % % %                         refPx = [refPx; wEntry.indices(:)];
% % % %                     end
% % % %                     refPx = unique(refPx);
% % % %                     referencePhase = mean(unwrapped(refPx), 'omitnan');
% % % %
% % % %                 case 'scenemean'
% % % %                     referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % % %
% % % %                 case 'none'
% % % %                     referencePhase = 0;
% % % %
% % % %                 otherwise
% % % %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% % % %             end
% % % %
% % % %             k = k + 1;
% % % %             insarData(k).pair = [i j];
% % % %             insarData(k).burst = b;
% % % %             insarData(k).phzWrapped = phzWrapped;
% % % %             insarData(k).coherence = cor;
% % % %             insarData(k).phzUnwrapped = unwrapped;
% % % %             insarData(k).phzReferenced = unwrapped - referencePhase;
% % % %             insarData(k).refCR = unwrapOpts.referenceMethod;
% % % %             insarData(k).refValue = referencePhase;
% % % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % % %         end
% % % %     end
% % % % end
% % % %
% % % % end
% % % %
% % % % % function insarData = process_interferometric_phase(sarData, crResult, unwrapOpts)
% % % % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % % % %   Inputs:
% % % % % %       sarData   - Struct of SLC data and metadata
% % % % % %       crResult  - CR calibration structure (from CRpower)
% % % % % %       unwrapOpts - Struct with fields:
% % % % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % % % %           .filterSize      : Interferogram filter window (default = 5)
% % % % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % % % %           .qualityThresh   : Mask threshold (default = 0.75)
% % % % %
% % % % % if nargin < 3, unwrapOpts = struct(); end
% % % % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.3; end
% % % % %
% % % % % numDirs = numel(sarData);
% % % % %
% % % % % insarData = struct([]);
% % % % % k = 0;
% % % % %
% % % % % for i = 1:numDirs
% % % % %     for j = i+1:numDirs
% % % % %         % Ensure bursts match
% % % % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % % % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % % % %             continue;
% % % % %         end
% % % % %
% % % % %         for b = 1:numel(sarData(i).slc)
% % % % %             slc1 = sarData(i).slc{b};
% % % % %             slc2 = sarData(j).slc{b};
% % % % %
% % % % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % % % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % % % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % % % %
% % % % %             % Inline unwrapping logic here
% % % % %             switch lower(unwrapOpts.method)
% % % % %                 case 'multiseed'
% % % % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % % % %                 case 'goldstein'
% % % % %                     unwrapped = insar.unwrap_goldstein(phzWrapped);
% % % % %                 otherwise
% % % % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % % % %             end
% % % % %
% % % % %             % Inline reference logic
% % % % %             switch lower(unwrapOpts.referenceMethod)
% % % % %                 case 'crweightedmean'
% % % % %                     refPhases = [];
% % % % %                     weights = [];
% % % % %                     for c = 1:numel(crResult)
% % % % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % % % %                         wEntry = crResult(c).Weights{i,j};
% % % % %                         if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
% % % % %                         px = wEntry.indices;
% % % % %                         wt = wEntry.value;
% % % % %                         if any(isnan(unwrapped(px))), continue; end
% % % % %                         refPhases(end+1) = sum(wt .* unwrapped(px));
% % % % %                         weights(end+1) = sum(wt);
% % % % %                     end
% % % % %                     if ~isempty(refPhases)
% % % % %                         referencePhase = sum(refPhases .* weights) / sum(weights);
% % % % %                     else
% % % % %                         referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % % % %                     end
% % % % %
% % % % %                 case 'crmean'
% % % % %                     refPx = [];
% % % % %                     for c = 1:numel(crResult)
% % % % %                         if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
% % % % %                         wEntry = crResult(c).Weights{i,j};
% % % % %                         if ~isfield(wEntry, 'indices'), continue; end
% % % % %                         refPx = [refPx; wEntry.indices(:)];
% % % % %                     end
% % % % %                     refPx = unique(refPx);
% % % % %                     referencePhase = mean(unwrapped(refPx), 'omitnan');
% % % % %
% % % % %                 case 'scenemean'
% % % % %                     referencePhase = mean(unwrapped(cor >= unwrapOpts.qualityThresh), 'omitnan');
% % % % %
% % % % %                 case 'none'
% % % % %                     referencePhase = 0;
% % % % %
% % % % %                 otherwise
% % % % %                     error('Unknown reference method: %s', unwrapOpts.referenceMethod);
% % % % %             end
% % % % %
% % % % %             k = k + 1;
% % % % %             insarData(k).pair = [i j];
% % % % %             insarData(k).burst = b;
% % % % %             insarData(k).phzWrapped = phzWrapped;
% % % % %             insarData(k).coherence = cor;
% % % % %             insarData(k).phzUnwrapped = unwrapped;
% % % % %             insarData(k).phzReferenced = unwrapped - referencePhase;
% % % % %             insarData(k).refCR = unwrapOpts.referenceMethod;
% % % % %             insarData(k).refValue = referencePhase;
% % % % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % % % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % % % %         end
% % % % %     end
% % % % % end
% % % % %
% % % % % end
% % % % %
% % % % % % function insarData = process_interferometric_phase(sarData, crResult, unwrapOpts)
% % % % % % %PROCESS_INTERFEROMETRIC_PHASE Computes wrapped, unwrapped, and referenced phase.
% % % % % % %   Inputs:
% % % % % % %       sarData   - Struct of SLC data and metadata
% % % % % % %       crResult  - CR calibration structure (from CRpower)
% % % % % % %       unwrapOpts - Struct with fields:
% % % % % % %           .method          : 'multiseed' (default), 'goldstein', etc.
% % % % % % %           .referenceMethod : 'crWeightedMean' (default), 'crMean', 'sceneMean', 'none'
% % % % % % %           .filterSize      : Interferogram filter window (default = 5)
% % % % % % %           .sigma           : Gaussian std dev (default = filterSize / 2)
% % % % % % %           .alpha           : Coherence weighting exponent (default = 2)
% % % % % % %           .qualityThresh   : Mask threshold (default = 0.75)
% % % % % %
% % % % % % if nargin < 3, unwrapOpts = struct(); end
% % % % % % if ~isfield(unwrapOpts, 'method'), unwrapOpts.method = 'multiseed'; end
% % % % % % if ~isfield(unwrapOpts, 'referenceMethod'), unwrapOpts.referenceMethod = 'crWeightedMean'; end
% % % % % % if ~isfield(unwrapOpts, 'filterSize'), unwrapOpts.filterSize = 5; end
% % % % % % if ~isfield(unwrapOpts, 'sigma'), unwrapOpts.sigma = unwrapOpts.filterSize / 2; end
% % % % % % if ~isfield(unwrapOpts, 'alpha'), unwrapOpts.alpha = 2; end
% % % % % % if ~isfield(unwrapOpts, 'qualityThresh'), unwrapOpts.qualityThresh = 0.75; end
% % % % % %
% % % % % % numDirs = numel(sarData);
% % % % % %
% % % % % % insarData = struct([]);
% % % % % % k = 0;
% % % % % %
% % % % % % for i = 1:numDirs
% % % % % %     for j = i+1:numDirs
% % % % % %         % Ensure bursts match
% % % % % %         if numel(sarData(i).slc) ~= numel(sarData(j).slc)
% % % % % %             warning('Mismatched burst count between dir %d and dir %d', i, j);
% % % % % %             continue;
% % % % % %         end
% % % % % %
% % % % % %         for b = 1:numel(sarData(i).slc)
% % % % % %             slc1 = sarData(i).slc{b};
% % % % % %             slc2 = sarData(j).slc{b};
% % % % % %
% % % % % %             [phzWrapped, cor] = insar.compute_interferogram(slc1, slc2, ...
% % % % % %                 unwrapOpts.qualityThresh, unwrapOpts.filterSize, true, ...
% % % % % %                 unwrapOpts.sigma, unwrapOpts.alpha);
% % % % % %
% % % % % %             % Unwrapping
% % % % % %             switch lower(unwrapOpts.method)
% % % % % %                 case 'multiseed'
% % % % % %                     unwrapped = insar.region_grow_unwrap_multiseed(phzWrapped, cor, unwrapOpts.qualityThresh);
% % % % % %                 case 'goldstein'
% % % % % %                     unwrapped = insar.unwrap_goldstein(phzWrapped);
% % % % % %                 otherwise
% % % % % %                     error('Unknown unwrap method: %s', unwrapOpts.method);
% % % % % %             end
% % % % % %             % Phase Referencing
% % % % % %             [refPhz, refCR, refVal] = insar.reference_phase(unwrapped, crResult, unwrapOpts.referenceMethod);
% % % % % %
% % % % % %             k = k + 1;
% % % % % %             insarData(k).pair = [i j];
% % % % % %             insarData(k).burst = b;
% % % % % %             insarData(k).phzWrapped = phzWrapped;
% % % % % %             insarData(k).coherence = cor;
% % % % % %             insarData(k).phzUnwrapped = unwrapped;
% % % % % %             insarData(k).phzReferenced = refPhz;
% % % % % %             insarData(k).refCR = refCR;
% % % % % %             insarData(k).refValue = refVal;
% % % % % %             insarData(k).unwrapMethod = unwrapOpts.method;
% % % % % %             insarData(k).refMethod = unwrapOpts.referenceMethod;
% % % % % %         end
% % % % % %     end
% % % % % % end
% % % % % %
% % % % % % end

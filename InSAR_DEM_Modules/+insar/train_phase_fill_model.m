function model = train_phase_fill_model(feat, opts)
%TRAIN_PHASE_FILL_MODEL Train lightweight model for dense phase completion.

if nargin < 2
    opts = struct();
end

if ~isfield(opts, 'trainCohThresh') || isempty(opts.trainCohThresh)
    opts.trainCohThresh = 0.65;
end
if ~isfield(opts, 'maxTrainSamples') || isempty(opts.maxTrainSamples)
    opts.maxTrainSamples = 150000;
end
if ~isfield(opts, 'modelType') || isempty(opts.modelType)
    opts.modelType = 'bag';
end
if ~isfield(opts, 'numTrees') || isempty(opts.numTrees)
    opts.numTrees = 100;
end
if ~isfield(opts, 'minLeafSize') || isempty(opts.minLeafSize)
    opts.minLeafSize = 10;
end
if ~isfield(opts, 'minValidFrac') || isempty(opts.minValidFrac)
    opts.minValidFrac = 0.25;   % keep feature if >=25% finite over candidate training rows
end
if ~isfield(opts, 'learnRate') || isempty(opts.learnRate)
    opts.learnRate = 0.05;
end
if ~isfield(opts, 'selectedFeatures') || isempty(opts.selectedFeatures)
    opts.selectedFeatures = {'Z','V','slope','aspect_sin','aspect_cos','power','cor','incidence'};
end

% Build matrix only from usable features
[Xmat, names, keepMask] = make_feature_matrix_safe(feat, opts);

y = feat.phz(:);

% Prefer stricter ML training mask when provided
if isfield(feat, 'trainMask') && ~isempty(feat.trainMask)
    baseTrainMask = feat.trainMask(:) & feat.validFeatureMask(:) & isfinite(y);
else
    baseTrainMask = feat.trustedMask(:) & feat.validFeatureMask(:) & isfinite(y);
end

if nnz(baseTrainMask) < 1000
    warning('Phase fill training mask too small; falling back to coherence-only training mask.');
    baseTrainMask = isfinite(y) & feat.validFeatureMask(:) & isfinite(feat.cor(:)) & ...
                    feat.cor(:) >= opts.trainCohThresh;
end

if isempty(Xmat)
    error('train_phase_fill_model:NoUsableFeatures', ...
        'No usable feature columns remained after filtering NaN-heavy predictors.');
end

trainMask = baseTrainMask & all(isfinite(Xmat), 2);

Xtrain = Xmat(trainMask, :);
ytrain = y(trainMask);

fprintf('Phase fill training candidates: %d\n', nnz(baseTrainMask));
fprintf('Phase fill retained rows:       %d\n', size(Xtrain,1));
fprintf('Phase fill retained features:   %d\n', size(Xtrain,2));

if isempty(Xtrain) || isempty(ytrain)
    error('train_phase_fill_model:EmptyTrainingSet', ...
        ['Phase fill training set is empty. ', ...
         'Try lowering trainCohThresh, checking feature coverage, ', ...
         'or reducing required optional predictors.']);
end

% Downsample if too many points
nTrain = size(Xtrain, 1);
if nTrain > opts.maxTrainSamples
    rng(1);
    idx = randperm(nTrain, opts.maxTrainSamples);
    Xtrain = Xtrain(idx, :);
    ytrain = ytrain(idx);
end

switch lower(opts.modelType)
    case 'bag'
        mdl = fitrensemble(Xtrain, ytrain, ...
            'Method', 'Bag', ...
            'NumLearningCycles', opts.numTrees, ...
            'Learners', templateTree('MinLeafSize', opts.minLeafSize));

    case 'lsboost'
        mdl = fitrensemble(Xtrain, ytrain, ...
            'Method', 'LSBoost', ...
            'NumLearningCycles', opts.numTrees, ...
            'LearnRate', opts.learnRate, ...
            'Learners', templateTree('MinLeafSize', opts.minLeafSize));

    otherwise
        error('Unknown modelType: %s', opts.modelType);
end

importance = [];
try
    importance = predictorImportance(mdl);
catch
end

model = struct();
model.mdl = mdl;
model.featureNames = names;
model.featureKeepMask = keepMask;
model.opts = opts;
model.importance = importance;
model.nTrain = size(Xtrain, 1);

end

% -------------------------------------------------------------------------
function [Xmat, names, keepMask] = make_feature_matrix_safe(feat, opts)

allNames = {'Z','V', 'slope','aspect_sin','aspect_cos','curvature','power','cor','localCoh','incidence','Bperp','slant'};
allCols = { ...
    feat.Z(:), feat.V(:), feat.slope(:), feat.aspect_sin(:), feat.aspect_cos(:), feat.curvature(:), ...
    feat.power(:), feat.cor(:), feat.localCoh(:), feat.incidence(:), feat.Bperp(:), feat.slant(:)};

keep = ismember(allNames, opts.selectedFeatures);
rawNames = allNames(keep);
rawCols = allCols(keep);

y = feat.phz(:);

if isfield(feat, 'trainMask') && ~isempty(feat.trainMask)
    candidateMask = feat.trainMask(:) & feat.validFeatureMask(:) & isfinite(y);
else
    candidateMask = feat.trustedMask(:) & feat.validFeatureMask(:) & isfinite(y);
end

if ~any(candidateMask)
    error('train_phase_fill_model:EmptyCandidateMask', ...
        'No valid training pixels remain after applying trainMask/trustedMask.');
end

keepMask = false(numel(rawCols),1);

for k = 1:numel(rawCols)
    col = rawCols{k};
    fracFinite = mean(isfinite(col(candidateMask)));
    if fracFinite >= opts.minValidFrac
        keepMask(k) = true;
    end
end

% Always keep core coordinates/terrain/coherence if possible
coreNames = {'Z','slope','aspect_sin','aspect_cos', ...
             'cor','power'};
for k = 1:numel(rawNames)
    if ismember(rawNames{k}, coreNames)
        if any(isfinite(rawCols{k}(candidateMask)))
            keepMask(k) = true;
        end
    end
end

names = rawNames(keepMask);
cols = rawCols(keepMask);

Xmat = [];
for k = 1:numel(cols)
    Xmat = [Xmat, cols{k}]; %#ok<AGROW>
end

end
% function model = train_phase_fill_model(feat, opts)
% %TRAIN_PHASE_FILL_MODEL Train lightweight model for dense phase completion.
% %
% % Target is referenced phase residual itself:
% %   y = feat.phz
% %
% % The screen has already been removed upstream, so phz is the referenced
% % residual phase we want to fill.
% 
% if nargin < 2
%     opts = struct();
% end
% 
% if ~isfield(opts, 'trainCohThresh') || isempty(opts.trainCohThresh)
%     opts.trainCohThresh = 0.65;
% end
% if ~isfield(opts, 'maxTrainSamples') || isempty(opts.maxTrainSamples)
%     opts.maxTrainSamples = 150000;
% end
% if ~isfield(opts, 'modelType') || isempty(opts.modelType)
%     opts.modelType = 'bag';
% end
% if ~isfield(opts, 'numTrees') || isempty(opts.numTrees)
%     opts.numTrees = 100;
% end
% if ~isfield(opts, 'minLeafSize') || isempty(opts.minLeafSize)
%     opts.minLeafSize = 10;
% end
% 
% [Xmat, names] = make_feature_matrix(feat);
% 
% y = feat.phz(:);
% trainMask = feat.trustedMask(:) & feat.validFeatureMask(:) & all(isfinite(Xmat), 2) & isfinite(y);
% 
% Xtrain = Xmat(trainMask, :);
% ytrain = y(trainMask);
% 
% % Downsample if too many points
% nTrain = size(Xtrain, 1);
% if nTrain > opts.maxTrainSamples
%     rng(1);
%     idx = randperm(nTrain, opts.maxTrainSamples);
%     Xtrain = Xtrain(idx, :);
%     ytrain = ytrain(idx);
% end
% 
% switch lower(opts.modelType)
%     case 'bag'
%         mdl = fitrensemble(Xtrain, ytrain, ...
%             'Method', 'Bag', ...
%             'NumLearningCycles', opts.numTrees, ...
%             'Learners', templateTree('MinLeafSize', opts.minLeafSize));
% 
%     case 'lsboost'
%         mdl = fitrensemble(Xtrain, ytrain, ...
%             'Method', 'LSBoost', ...
%             'NumLearningCycles', opts.numTrees, ...
%             'Learners', templateTree('MinLeafSize', opts.minLeafSize));
% 
%     otherwise
%         error('Unknown modelType: %s', opts.modelType);
% end
% 
% % OOB importance if available
% importance = [];
% try
%     importance = predictorImportance(mdl);
% catch
% end
% 
% model = struct();
% model.mdl = mdl;
% model.featureNames = names;
% model.opts = opts;
% model.importance = importance;
% model.nTrain = size(Xtrain, 1);
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
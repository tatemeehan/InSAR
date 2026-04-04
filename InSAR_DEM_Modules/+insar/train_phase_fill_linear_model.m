function model = train_phase_fill_linear_model(feat, featureNames)
%TRAIN_PHASE_FILL_LINEAR_MODEL
% Fit simple linear regression on trusted pixels.

if nargin < 2 || isempty(featureNames)
    featureNames = {'Z','slope','aspect_sin','aspect_cos','power','cor'};
end

[Xmat, names] = make_feature_matrix_selected(feat, featureNames);
y = feat.phz(:);

trainMask = feat.trustedMask(:) & feat.validFeatureMask(:) & isfinite(y) & all(isfinite(Xmat),2);

Xtrain = Xmat(trainMask,:);
ytrain = y(trainMask);

% add intercept
Xaug = [ones(size(Xtrain,1),1), Xtrain];

beta = Xaug \ ytrain;

model = struct();
model.beta = beta;
model.featureNames = names;

end

% -------------------------------------------------------------------------
function [Xmat, names] = make_feature_matrix_selected(feat, featureNames)

allNames = {'Z','slope','aspect_sin','aspect_cos','power','cor','curvature','localCoh'};
allCols = { ...
    feat.Z(:), feat.slope(:), feat.aspect_sin(:), feat.aspect_cos(:), ...
    feat.power(:), feat.cor(:), feat.curvature(:), feat.localCoh(:)};

keep = ismember(allNames, featureNames);
names = allNames(keep);
cols = allCols(keep);

Xmat = [];
for k = 1:numel(cols)
    Xmat = [Xmat, cols{k}]; %#ok<AGROW>
end

end
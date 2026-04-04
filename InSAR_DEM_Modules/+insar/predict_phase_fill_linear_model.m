function pred = predict_phase_fill_linear_model(feat, model)
%PREDICT_PHASE_FILL_LINEAR_MODEL
% Predict referenced phase using linear model.

[Xmat, ~] = make_feature_matrix_selected(feat, model.featureNames);

valid = all(isfinite(Xmat),2);
yhat = nan(size(feat.phz(:)));

Xaug = [ones(nnz(valid),1), Xmat(valid,:)];
yhat(valid) = Xaug * model.beta;

pred = reshape(yhat, size(feat.phz));

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
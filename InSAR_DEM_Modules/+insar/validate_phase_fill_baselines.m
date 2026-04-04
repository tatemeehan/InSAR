function out = validate_phase_fill_baselines(feat, holdoutMask, featureNames, sigma)
%VALIDATE_PHASE_FILL_BASELINES
% Evaluate linear and Gaussian interpolation baselines on a given holdout.

if nargin < 3 || isempty(featureNames)
    featureNames = {'Z','slope','aspect_sin','aspect_cos','power','cor'};
end
if nargin < 4 || isempty(sigma)
    sigma = 8;
end

trustedMask0 = logical(feat.trustedMask);
trainMask = trustedMask0 & ~holdoutMask;

featTrain = feat;
featTrain.trustedMask = trainMask;

yTrue = feat.phz(holdoutMask);

% --- Linear model
linModel = insar.train_phase_fill_linear_model(featTrain, featureNames);
linPredMap = insar.predict_phase_fill_linear_model(feat, linModel);
yPredLin = linPredMap(holdoutMask);
validLin = isfinite(yTrue) & isfinite(yPredLin);
linMetrics = compute_metrics(yTrue(validLin), yPredLin(validLin));

% --- Gaussian interpolation
gaussPredMap = insar.phase_fill_gaussian_interp(feat.phz, trainMask, sigma);
yPredGauss = gaussPredMap(holdoutMask);
validGauss = isfinite(yTrue) & isfinite(yPredGauss);
gaussMetrics = compute_metrics(yTrue(validGauss), yPredGauss(validGauss));

out = struct();
out.trainMask = trainMask;
out.holdoutMask = holdoutMask;
out.linear.predMap = linPredMap;
out.linear.metrics = linMetrics;
out.gaussian.predMap = gaussPredMap;
out.gaussian.metrics = gaussMetrics;

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
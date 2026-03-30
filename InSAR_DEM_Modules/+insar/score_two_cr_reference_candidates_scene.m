function out = score_two_cr_reference_candidates_scene(obsV, pred, phzUnwrapped, cor, opts)
% Compare crAzRamp vs crRangeRamp using whole-scene residual nuisance score.

out = struct();
out.bestMethod = 'crAzRamp';
out.azScore = inf;
out.rangeScore = inf;
out.azDetail = struct();
out.rangeDetail = struct();
out.reason = '';

if numel(obsV) ~= 2
    out.reason = 'Requires exactly 2 valid CRs.';
    return;
end

% Fit both candidate models
[refAz, ~, coeffAz] = insar.fit_cr_sensor_ramp(obsV, pred.azCoord, 'az');
[refRg, ~, coeffRg] = insar.fit_cr_sensor_ramp(obsV, pred.rangeCoord, 'range');
refAvg = 0.5 * (refAz + refRg);

% Residual scenes
resAz = phzUnwrapped - refAz;
resRg = phzUnwrapped - refRg;
resAvg = phzUnwrapped - refAvg;

% Score both residuals
[out.azScore, out.azDetail] = local_scene_score(resAz, pred, cor, opts);
[out.rangeScore, out.rangeDetail] = local_scene_score(resRg, pred, cor, opts);
[out.avgScore, out.avgDetail] = local_scene_score(resAvg, pred, cor, opts);

scores = [out.azScore,out.rangeScore,out.avgScore];
% Pick winner
if all(isfinite(scores))
    [~, minIx] = min(scores);
    if minIx == 1
        out.bestMethod = 'crAzRamp';
    elseif minIx == 2
        out.bestMethod = 'crRangeRamp';
    elseif minIx == 3
        out.bestMethod = 'crAvgRamp';
    else
        out.bestMethod = local_dominant_choice(obsV);
        out.reason = 'Scene scoring unavailable; fell back to dominant.';
        return;
    end
    out.reason = 'Chosen by whole-scene residual nuisance score.';
else
    out.bestMethod = local_dominant_choice(obsV);
    out.reason = 'Scene scoring unavailable; fell back to dominant.';
    return;
end

% Near-tie fallback    
[~, minIx] = mink(scores,2);
relDiff = abs(scores(minIx(1)) - scores(minIx(2))) / max(min(minIx(1), minIx(2)), eps);
if relDiff < 0.05
    out.bestMethod = local_dominant_choice(obsV);
    out.reason = 'Scene scores nearly tied; fell back to dominant.';
end

end

function [score, detail] = local_scene_score(res, pred, cor, opts)
detail = struct('corrAz', NaN, 'corrRg', NaN, 'planeRMS', NaN, 'nPix', 0);

mask = pred.mask & isfinite(res) & isfinite(pred.azCoord) & isfinite(pred.rangeCoord);

if isfield(opts, 'qualityThresh') && ~isempty(cor)
    mask = mask & cor >= opts.qualityThresh;
end

idx = find(mask);
detail.nPix = numel(idx);

if numel(idx) < 20
    score = inf;
    return;
end

phi = res(idx);
az  = pred.azCoord(idx);
rg  = pred.rangeCoord(idx);

% Remove mean first
phi = phi - mean(phi, 'omitnan');

% Remaining correlation with sensor axes
detail.corrAz = abs(local_corr(phi, az));
detail.corrRg = abs(local_corr(phi, rg));

% Residual plane energy in sensor coords
A = [ones(numel(phi),1), az, rg];
beta = A \ phi;
phiPlane = A * beta;
detail.planeRMS = sqrt(mean(phiPlane.^2, 'omitnan'));

% Composite score
w1 = 1.0;
w2 = 1.0;
w3 = 1.0;
score = w1 * detail.corrAz + w2 * detail.corrRg + w3 * detail.planeRMS;

end

function r = local_corr(x, y)
x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);

if numel(x) < 3
    r = NaN;
    return;
end

x = x - mean(x, 'omitnan');
y = y - mean(y, 'omitnan');

den = sqrt(sum(x.^2) * sum(y.^2));
if den < eps
    r = NaN;
else
    r = sum(x .* y) / den;
end
end

function method = local_dominant_choice(obsV)
azVals = [obsV.az];
rgVals = [obsV.rg];

daz = max(azVals) - min(azVals);
drg = max(rgVals) - min(rgVals);

if abs(daz) >= abs(drg)
    method = 'crAzRamp';
else
    method = 'crRangeRamp';
end
end
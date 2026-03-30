function out = score_two_cr_reference_candidates(obsV, pred, phzUnwrapped, opts)
% Compare crAzRamp vs crRangeRamp using CR support pixels.
%
% Returns:
%   out.bestMethod
%   out.azScore
%   out.rangeScore
%   out.azCoeffs
%   out.rangeCoeffs
%   out.reason

out = struct();
out.bestMethod = 'crAzRamp';
out.azScore = inf;
out.rangeScore = inf;
out.azCoeffs = [];
out.rangeCoeffs = [];
out.reason = '';

if numel(obsV) ~= 2
    out.reason = 'Requires exactly 2 valid CRs.';
    return;
end

% Fit both candidate models
[~, ~, coeffAz] = insar.fit_cr_sensor_ramp(obsV, pred.azCoord, 'az');
[~, ~, coeffRg] = insar.fit_cr_sensor_ramp(obsV, pred.rangeCoord, 'range');

out.azCoeffs = coeffAz;
out.rangeCoeffs = coeffRg;

% Score each candidate on CR support pixels
out.azScore = local_score_model(obsV, phzUnwrapped, pred.azCoord, coeffAz);
out.rangeScore = local_score_model(obsV, phzUnwrapped, pred.rangeCoord, coeffRg);

% Choose lower score
if isfinite(out.azScore) && isfinite(out.rangeScore)
    if out.azScore <= out.rangeScore
        out.bestMethod = 'crAzRamp';
    else
        out.bestMethod = 'crRangeRamp';
    end
    out.reason = 'Chosen by lower weighted CR-support residual.';
else
    % fallback: leverage rule
    azVals = [obsV.az];
    rgVals = [obsV.rg];
    daz = max(azVals) - min(azVals);
    drg = max(rgVals) - min(rgVals);

    if abs(daz) >= abs(drg)
        out.bestMethod = 'crAzRamp';
    else
        out.bestMethod = 'crRangeRamp';
    end
    out.reason = 'Residual scoring unavailable; fell back to leverage rule.';
end

% If nearly tied, also fall back to leverage rule for stability
if isfinite(out.azScore) && isfinite(out.rangeScore)
    relDiff = abs(out.azScore - out.rangeScore) / max(min(out.azScore, out.rangeScore), eps);
    if relDiff < 0.05
        azVals = [obsV.az];
        rgVals = [obsV.rg];
        daz = max(azVals) - min(azVals);
        drg = max(rgVals) - min(rgVals);

        if abs(daz) >= abs(drg)
            out.bestMethod = 'crAzRamp';
        else
            out.bestMethod = 'crRangeRamp';
        end
        out.reason = 'Scores nearly tied; fell back to leverage rule.';
    end
end

end

function score = local_score_model(obsV, phzUnwrapped, coordGrid, coeffs)
scoreNum = 0;
scoreDen = 0;

for k = 1:numel(obsV)
    px = obsV(k).supportPx(:);
    wt = obsV(k).supportWt(:);

    good = isfinite(px) & px >= 1 & px <= numel(phzUnwrapped);
    px = px(good);
    wt = wt(good);

    if isempty(px)
        continue;
    end

    good2 = isfinite(phzUnwrapped(px)) & isfinite(coordGrid(px));
    px = px(good2);
    wt = wt(good2);

    if isempty(px)
        continue;
    end

    wt = wt / max(sum(wt), eps);

    phiPred = coeffs(1) + coeffs(2) * coordGrid(px);
    resid = phzUnwrapped(px) - phiPred;

    % weighted RMS-like score
    scoreNum = scoreNum + sum(wt .* resid.^2);
    scoreDen = scoreDen + sum(wt);
end

if scoreDen > 0
    score = sqrt(scoreNum / scoreDen);
else
    score = inf;
end
end
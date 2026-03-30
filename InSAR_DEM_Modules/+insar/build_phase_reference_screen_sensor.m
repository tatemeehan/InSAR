function [refScreen, refScreenWrapped, meta] = build_phase_reference_screen_sensor( ...
    phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, i, b1, j, b2, opts)

sz = size(phzWrapped);
refScreen = zeros(sz);
refScreenWrapped = zeros(sz);

meta = struct();
meta.methodRequested = opts.referenceMethod;
meta.methodUsed = 'none';
meta.numValidCR = 0;
meta.validCRNames = {};
meta.coeffs = [];
meta.residuals = [];
meta.scalarEquivalent = 0;
meta.scalarEquivalentWrapped = 0;
meta.reason = '';
meta.twoCRDecision = '';
meta.obs = [];

obs = insar.extract_cr_phase_observations_sensor( ...
    phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, i, b1, j, b2, opts);

valid = [obs.isValid];
obsV = obs(valid);
meta.numValidCR = numel(obsV);
meta.validCRNames = {obsV.name};
meta.obs = obsV;

requested = lower(opts.referenceMethod);

if strcmp(requested, 'auto')
    if meta.numValidCR >= 3
        requested = 'crsensorplane';
    elseif meta.numValidCR == 2
        if isfield(opts, 'referenceAutoTwoCRMode') && strcmpi(opts.referenceAutoTwoCRMode, 'score')
            % twoCR = insar.score_two_cr_reference_candidates(obsV, pred, phzUnwrapped, opts);
            twoCR = insar.score_two_cr_reference_candidates_scene(obsV, pred, phzUnwrapped, cor, opts);
            requested = lower(twoCR.bestMethod);
            meta.twoCRDecision = requested;
            meta.twoCRScore = twoCR;
        else
            requested = local_pick_two_cr_mode(obsV, opts);
            meta.twoCRDecision = requested;
        end
    elseif meta.numValidCR == 1
        requested = 'crconstant';
    else
        requested = lower(opts.referenceFallback);
    end
end

switch requested
    case 'none'
        meta.methodUsed = 'none';
        return;

    case 'scenemean'
        [refScreen, refScreenWrapped, meta] = local_scene_mean(phzWrapped, phzUnwrapped, cor, opts, meta);
        return;

    case 'crconstant'
        if meta.numValidCR < 1
            [refScreen, refScreenWrapped, meta] = local_scene_mean(phzWrapped, phzUnwrapped, cor, opts, meta);
            meta.reason = 'No valid CRs for crConstant';
            return;
        end

        w = [obsV.weight]';
        phiU = [obsV.phiUnwrapped]';
        phiW = [obsV.phiWrapped]';

        sU = sum(w .* phiU) / max(sum(w), eps);
        sW = angle(sum(w .* exp(1i * phiW)) / max(sum(w), eps));

        refScreen(:) = sU;
        refScreenWrapped(:) = sW;

        meta.methodUsed = 'crConstant';
        meta.scalarEquivalent = sU;
        meta.scalarEquivalentWrapped = sW;
        return;

    case 'crazramp'
        if meta.numValidCR < 2
            opts2 = opts; opts2.referenceMethod = 'crConstant';
            [refScreen, refScreenWrapped, meta] = insar.build_phase_reference_screen_sensor( ...
                phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, i, b1, j, b2, opts2);
            meta.reason = 'Fell back from crAzRamp to crConstant';
            return;
        end

        [refScreen, refScreenWrapped, coeffs, residuals] = ...
            insar.fit_cr_sensor_ramp(obsV, pred.azCoord, 'az');

        meta.methodUsed = 'crAzRamp';
        meta.coeffs = coeffs;
        meta.residuals = residuals;
        meta.scalarEquivalent = coeffs(1);
        meta.scalarEquivalentWrapped = angle(mean(exp(1i * refScreenWrapped(:)), 'omitnan'));
        return;

    case 'crrangeramp'
        if meta.numValidCR < 2
            opts2 = opts; opts2.referenceMethod = 'crConstant';
            [refScreen, refScreenWrapped, meta] = insar.build_phase_reference_screen_sensor( ...
                phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, i, b1, j, b2, opts2);
            meta.reason = 'Fell back from crRangeRamp to crConstant';
            return;
        end

        [refScreen, refScreenWrapped, coeffs, residuals] = ...
            insar.fit_cr_sensor_ramp(obsV, pred.rangeCoord, 'range');

        meta.methodUsed = 'crRangeRamp';
        meta.coeffs = coeffs;
        meta.residuals = residuals;
        meta.scalarEquivalent = coeffs(1);
        meta.scalarEquivalentWrapped = angle(mean(exp(1i * refScreenWrapped(:)), 'omitnan'));
        return;

    case 'cravgramp'
        [refAz, refAzW, coeffsAz, residualsAz] = insar.fit_cr_sensor_ramp(obsV, pred.azCoord, 'az');
        [refRg, refRgW, coeffsRg, residualsRg] = insar.fit_cr_sensor_ramp(obsV, pred.rangeCoord, 'range');
        refScreen = 0.5 * (refAz + refRg);
        refScreenWrapped = angle(exp(1i * 0.5 * (angle(exp(1i*refAzW)) + angle(exp(1i*refRgW)))));

        meta.methodUsed = 'crAvgRamp';
        meta.coeffs = mean([coeffsAz,coeffsRg]);
        meta.residuals = mean([residualsAz,residualsRg]);
        meta.scalarEquivalent = mean([coeffsAz(1),coeffsRg(1)]);
        meta.scalarEquivalentWrapped = angle(mean(exp(1i * refScreenWrapped(:)), 'omitnan'));

    case 'crsensorplane'
        if meta.numValidCR < 3
            if meta.numValidCR == 2
                opts2 = opts;
                opts2.referenceMethod = local_pick_two_cr_mode(obsV, opts);
            elseif meta.numValidCR == 1
                opts2 = opts;
                opts2.referenceMethod = 'crConstant';
            else
                opts2 = opts;
                opts2.referenceMethod = lower(opts.referenceFallback);
            end

            [refScreen, refScreenWrapped, meta] = insar.build_phase_reference_screen_sensor( ...
                phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, i, b1, j, b2, opts2);
            meta.reason = 'Fell back from crSensorPlane';
            return;
        end

        [refScreen, refScreenWrapped, coeffs, residuals] = ...
            insar.fit_cr_sensor_plane(obsV, pred.azCoord, pred.rangeCoord, opts.referenceUseRobustFit);

        meta.methodUsed = 'crSensorPlane';
        meta.coeffs = coeffs;
        meta.residuals = residuals;
        meta.scalarEquivalent = coeffs(1);
        meta.scalarEquivalentWrapped = angle(mean(exp(1i * refScreenWrapped(:)), 'omitnan'));
        return;

    otherwise
        error('Unknown reference method: %s', opts.referenceMethod);
end
end

function method = local_pick_two_cr_mode(obsV, opts)
mode = lower(opts.referenceAutoTwoCRMode);

switch mode
    case 'az'
        method = 'crAzRamp';

    case 'range'
        method = 'crRangeRamp';

    case 'dominant'
        azVals = [obsV.az];
        rgVals = [obsV.rg];

        daz = max(azVals) - min(azVals);
        drg = max(rgVals) - min(rgVals);

        if abs(daz) >= abs(drg)
            method = 'crAzRamp';
        else
            method = 'crRangeRamp';
        end

    otherwise
        error('Unknown referenceAutoTwoCRMode: %s', opts.referenceAutoTwoCRMode);
end
method = lower(method);
end

function [refScreen, refScreenWrapped, meta] = local_scene_mean(phzWrapped, phzUnwrapped, cor, opts, meta)
maskU = cor >= opts.qualityThresh & isfinite(phzUnwrapped);
maskW = cor >= opts.qualityThresh & isfinite(phzWrapped);

sU = mean(phzUnwrapped(maskU), 'omitnan');
sW = angle(mean(exp(1i * phzWrapped(maskW)), 'omitnan'));

refScreen = sU * ones(size(phzWrapped));
refScreenWrapped = sW * ones(size(phzWrapped));

meta.methodUsed = 'sceneMean';
meta.scalarEquivalent = sU;
meta.scalarEquivalentWrapped = sW;
end
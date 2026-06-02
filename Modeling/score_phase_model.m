function stats = score_phase_model(phiData, phiModel, fitMask, evalMask)
%SCORE_PHASE_MODEL Score wrapped model phase against wrapped data phase.
%
% This reports both raw residual metrics and scalar-bias-corrected metrics.
%
% Inputs
%   phiData  : observed wrapped phase [rad]
%   phiModel : modeled wrapped phase [rad]
%   fitMask  : mask used to estimate scalar circular bias
%   evalMask : mask used to evaluate residual statistics
%
% Output
%   stats : structure with bias, raw metrics, bias-corrected metrics,
%           residual rasters, and bias-corrected model phase.

wrap = @(x) angle(exp(1i .* double(x)));

phiData  = double(phiData);
phiModel = double(phiModel);

fitMask  = logical(fitMask);
evalMask = logical(evalMask);

if ~isequal(size(phiData), size(phiModel))
    error('phiData and phiModel must have the same size.');
end

if ~isequal(size(fitMask), size(phiData))
    error('fitMask must have the same size as phiData.');
end

if ~isequal(size(evalMask), size(phiData))
    error('evalMask must have the same size as phiData.');
end

valid0 = isfinite(phiData) & isfinite(phiModel);

% -------------------------
% Raw residual
% -------------------------
resRaw = wrap(phiData - phiModel);

evalValidRaw = evalMask & valid0 & isfinite(resRaw);
rvRaw = resRaw(evalValidRaw);

Rraw = abs(mean(exp(1i .* rvRaw), 'omitnan'));
circStdRaw = sqrt(-2 .* log(Rraw));
medAbsRaw = median(abs(rvRaw), 'omitnan');

% Complex phasor correlation / coherence-like score
zData_raw  = exp(1i .* phiData(evalValidRaw));
zModel_raw = exp(1i .* phiModel(evalValidRaw));

phaseCorr_complex_raw = abs(sum(zData_raw .* conj(zModel_raw), 'omitnan')) ./ ...
    sqrt(sum(abs(zData_raw).^2, 'omitnan') .* sum(abs(zModel_raw).^2, 'omitnan'));

% -------------------------
% Scalar circular bias fit
% -------------------------
fitValid = fitMask & valid0 & isfinite(resRaw);

if nnz(fitValid) == 0
    bias = NaN;
else
    bias = angle(mean(exp(1i .* resRaw(fitValid)), 'omitnan'));
end

phiModelBias = wrap(phiModel + bias);
resBias = wrap(phiData - phiModelBias);

evalValidBias = evalMask & valid0 & isfinite(resBias);
rvBias = resBias(evalValidBias);

Rbias = abs(mean(exp(1i .* rvBias), 'omitnan'));
circStdBias = sqrt(-2 .* log(Rbias));
medAbsBias = median(abs(rvBias), 'omitnan');

% Bias-corrected complex phasor correlation
zData_bias  = exp(1i .* phiData(evalValidBias));
zModel_bias = exp(1i .* phiModelBias(evalValidBias));

phaseCorr_complex_bias = abs(sum(zData_bias .* conj(zModel_bias), 'omitnan')) ./ ...
    sqrt(sum(abs(zData_bias).^2, 'omitnan') .* sum(abs(zModel_bias).^2, 'omitnan'));

% -------------------------
% Package
% -------------------------
stats = struct();

stats.bias = bias;

stats.R_raw = Rraw;
stats.circStd_raw = circStdRaw;
stats.medAbs_raw = medAbsRaw;
stats.phaseCorr_complex_raw = phaseCorr_complex_raw;


stats.R_bias = Rbias;
stats.circStd_bias = circStdBias;
stats.medAbs_bias = medAbsBias;
stats.phaseCorr_complex_bias = phaseCorr_complex_bias;


stats.nFit = nnz(fitValid);
stats.nEval = nnz(evalValidBias);

stats.resRaw = single(resRaw);
stats.resBias = single(resBias);
stats.phiModelBias = single(phiModelBias);

end
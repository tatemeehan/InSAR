function stats = phase_model_bias_stats(phiData, phiModel, fitMask, evalMask)
%PHASE_MODEL_BIAS_STATS Scalar-align model phase to data and compute residual stats.
%
% Inputs
%   phiData  : observed/reference phase [rad]
%   phiModel : modeled phase [rad]
%   fitMask  : logical mask used to estimate scalar circular bias
%   evalMask : logical mask used to evaluate residual statistics
%
% Output
%   stats.bias          : scalar circular phase correction added to model
%   stats.R             : circular residual concentration
%   stats.circStd       : circular standard deviation [rad]
%   stats.medAbs        : median absolute wrapped residual [rad]
%   stats.nFit          : number of pixels used for bias fit
%   stats.nEval         : number of pixels used for evaluation
%   stats.phiModelBias  : bias-corrected model phase
%   stats.residual      : wrapped residual, data - bias-corrected model

wrap = @(x) angle(exp(1i .* double(x)));

phiData  = double(phiData);
phiModel = double(phiModel);

if nargin < 3 || isempty(fitMask)
    fitMask = isfinite(phiData) & isfinite(phiModel);
end

if nargin < 4 || isempty(evalMask)
    evalMask = fitMask;
end

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

% Initial residual
res0 = wrap(phiData - phiModel);

% Fit scalar circular bias over fit mask
fitValid = fitMask & isfinite(res0);

if nnz(fitValid) == 0
    warning('No valid pixels in fitMask.');
    bias = NaN;
else
    bias = angle(mean(exp(1i .* res0(fitValid)), 'omitnan'));
end

% Apply model bias
phiModelBias = wrap(phiModel + bias);

% Final residual
resBias = wrap(phiData - phiModelBias);

% Evaluate residuals
evalValid = evalMask & isfinite(resBias);

rv = resBias(evalValid);

if isempty(rv)
    R = NaN;
    circStd = NaN;
    medAbs = NaN;
else
    R = abs(mean(exp(1i .* rv), 'omitnan'));
    circStd = sqrt(-2 .* log(R));
    medAbs = median(abs(rv), 'omitnan');
end

stats = struct();
stats.bias = bias;
stats.R = R;
stats.circStd = circStd;
stats.medAbs = medAbs;
stats.nFit = nnz(fitValid);
stats.nEval = nnz(evalValid);
stats.phiModelBias = single(phiModelBias);
stats.residual = single(resBias);

end
function stats = pearson_phase_bias_corrected(phiData, phiModel, fitMask, evalMask)

wrap = @(x) angle(exp(1i .* x));

fitValid = fitMask & isfinite(phiData) & isfinite(phiModel);
evalValid = evalMask & isfinite(phiData) & isfinite(phiModel);

resFit = wrap(phiData - phiModel);
bias = angle(mean(exp(1i .* resFit(fitValid)), 'omitnan'));

phiModelBias = wrap(phiModel + bias);

x = double(phiData(evalValid));
y = double(phiModelBias(evalValid));

% Direct Pearson on bias-corrected wrapped phase
Rxy = corrcoef(x, y, 'Rows', 'complete');
r_wrapped_bias = Rxy(1,2);

% Safer Pearson on unit-circle components
cu = cos(x); su = sin(x);
cm = cos(y); sm = sin(y);

Rc = corrcoef(cu, cm, 'Rows', 'complete');
Rs = corrcoef(su, sm, 'Rows', 'complete');

stats.bias = bias;
stats.r_wrapped_bias = r_wrapped_bias;
stats.r_cos_bias = Rc(1,2);
stats.r_sin_bias = Rs(1,2);
stats.r_component_mean_bias = mean([stats.r_cos_bias, stats.r_sin_bias], 'omitnan');

end
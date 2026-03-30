function [refScreen, refScreenWrapped, coeffs, residuals] = fit_cr_sensor_plane(obs, azGrid, rgGrid, useRobustFit)
%FIT_CR_SENSOR_PLANE Fit phi = a + b*az + c*rg

if numel(obs) < 3
    error('fit_cr_sensor_plane requires at least 3 CR observations.');
end
if numel(obs) <= 4
    useRobustFit = false;
end

az = [obs.az]';
rg = [obs.rg]';
phiU = [obs.phiUnwrapped]';
phiW = [obs.phiWrapped]';
w    = [obs.weight]';

A = [ones(numel(az),1), az, rg];

% unwrapped fit
if useRobustFit && exist('robustfit', 'file') == 2
    bU = robustfit([az, rg], phiU);
else
    W = diag(w / max(sum(w), eps));
    bU = (A' * W * A) \ (A' * W * phiU);
end

% wrapped fit relative to highest-weight anchor
[~, ia] = max(w);
phiW_rel = phiW(ia) + angle(exp(1i * (phiW - phiW(ia))));

if useRobustFit && exist('robustfit', 'file') == 2
    bW = robustfit([az, rg], phiW_rel);
else
    W = diag(w / max(sum(w), eps));
    bW = (A' * W * A) \ (A' * W * phiW_rel);
end

coeffs = bU;
refScreen = bU(1) + bU(2) * azGrid + bU(3) * rgGrid;
refScreenWrapped = angle(exp(1i * (bW(1) + bW(2) * azGrid + bW(3) * rgGrid)));

pred = A * bU;
residuals = phiU - pred;
end
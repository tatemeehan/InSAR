function [refScreen, refScreenWrapped, coeffs, residuals] = fit_cr_sensor_ramp(obs, coordGrid, mode)
%FIT_CR_SENSOR_RAMP Fit phi = a + b*u, where u is azimuth or range predictor.
%
% obs must contain fields:
%   phiUnwrapped, phiWrapped, weight, az, rg

if numel(obs) < 2
    error('fit_cr_sensor_ramp requires at least 2 CR observations.');
end

switch lower(mode)
    case 'az'
        u = [obs.az]';
    case 'range'
        u = [obs.rg]';
    otherwise
        error('Unknown mode: %s', mode);
end

phiU = [obs.phiUnwrapped]';
phiW = [obs.phiWrapped]';
w    = [obs.weight]';

A = [ones(numel(u),1), u];

% weighted LS for unwrapped phase
W = diag(w / max(sum(w), eps));
coeffs = (A' * W * A) \ (A' * W * phiU);

% wrapped: unwrap relative to best CR anchor
[~, ia] = max(w);
phiW_rel = phiW(ia) + angle(exp(1i * (phiW - phiW(ia))));
coeffsW = (A' * W * A) \ (A' * W * phiW_rel);

refScreen = coeffs(1) + coeffs(2) * coordGrid;
refScreenWrapped = angle(exp(1i * (coeffsW(1) + coeffsW(2) * coordGrid)));

pred = A * coeffs;
residuals = phiU - pred;
end
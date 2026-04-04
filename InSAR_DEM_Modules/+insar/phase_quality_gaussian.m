function [quality, supportFrac] = phase_quality_gaussian(wrapped_phase, sigma, minValidFrac, validMask)
%PHASE_QUALITY_GAUSSIAN
% NaN-/mask-aware local wrapped-phase quality using Gaussian normalized convolution.
%
% Inputs
%   wrapped_phase : wrapped phase [rad]
%   validMask     : logical mask of pixels that belong to the interferogram
%   sigma         : Gaussian sigma in pixels (e.g. 4-6)
%   minValidFrac  : minimum local support fraction required (e.g. 0.35)
%
% Outputs
%   quality       : [0,1] local phase concentration, NaN outside valid support
%   supportFrac   : local valid-support fraction under the kernel

if nargin < 4 || isempty(validMask)
    validMask = isfinite(wrapped_phase);
end
if nargin < 2 || isempty(sigma)
    sigma = 5;
end
if nargin < 3 || isempty(minValidFrac)
    minValidFrac = 0.35;
end

phi = double(wrapped_phase);
validMask = logical(validMask) & isfinite(phi);

% Set invalid pixels to NaN explicitly
phi(~validMask) = NaN;

% Complex phase
z = exp(1i * phi);

% Gaussian kernel
halfSize = max(2, ceil(3*sigma));
x = -halfSize:halfSize;
[X, Y] = meshgrid(x, x);
K = exp(-(X.^2 + Y.^2) / (2*sigma^2));
K = K / sum(K(:));

% Zero-fill invalid pixels for normalized convolution
zr = real(z); zr(~validMask) = 0;
zi = imag(z); zi(~validMask) = 0;

num_r = conv2(zr, K, 'same');
num_i = conv2(zi, K, 'same');

% Since K sums to 1, this is local valid fraction in [0,1]
supportFrac = conv2(double(validMask), K, 'same');

quality = nan(size(phi));

% Only compute where:
%   1) center pixel is valid
%   2) neighborhood has enough valid support
good = validMask & (supportFrac >= minValidFrac);

mean_r = nan(size(phi));
mean_i = nan(size(phi));
mean_r(good) = num_r(good) ./ supportFrac(good);
mean_i(good) = num_i(good) ./ supportFrac(good);

quality(good) = hypot(mean_r(good), mean_i(good));

% Numerical clamp
quality = max(0, min(1, quality));

% Hard-mask invalid center pixels and weak-support pixels
quality(~validMask) = NaN;
quality(supportFrac < minValidFrac) = NaN;

end
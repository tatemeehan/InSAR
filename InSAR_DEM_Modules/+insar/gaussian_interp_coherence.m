function cohFilled = gaussian_interp_coherence(coh, trainMask, sigma)
%GAUSSIAN_INTERP_COHERENCE
% Gaussian normalized interpolation of coherence from trusted pixels.

if nargin < 3 || isempty(sigma)
    sigma = 2;
end

C = double(coh);
trainMask = logical(trainMask) & isfinite(C);

C0 = C;
C0(~trainMask) = 0;

halfSize = max(2, ceil(3*sigma));
x = -halfSize:halfSize;
[X, Y] = meshgrid(x, x);
K = exp(-(X.^2 + Y.^2) / (2*sigma^2));
K = K / sum(K(:));

num = conv2(C0, K, 'same');
den = conv2(double(trainMask), K, 'same');

cohFilled = nan(size(C));
good = den > 1e-6;
cohFilled(good) = num(good) ./ den(good);

cohFilled = max(0, min(1, cohFilled));

% keep invalid center pixels as NaN
cohFilled(~isfinite(C) & ~trainMask) = nan;

end
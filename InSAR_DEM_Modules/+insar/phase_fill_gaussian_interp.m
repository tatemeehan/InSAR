function pred = phase_fill_gaussian_interp(phz, trainMask, sigma)
%PHASE_FILL_GAUSSIAN_INTERP
% Gaussian normalized interpolation from trusted pixels only.

if nargin < 3 || isempty(sigma)
    sigma = 8;
end

phi = double(phz);
trainMask = logical(trainMask) & isfinite(phi);

phi0 = phi;
phi0(~trainMask) = 0;

halfSize = max(2, ceil(3*sigma));
x = -halfSize:halfSize;
[X, Y] = meshgrid(x, x);
K = exp(-(X.^2 + Y.^2)/(2*sigma^2));
K = K / sum(K(:));

num = conv2(phi0, K, 'same');
den = conv2(double(trainMask), K, 'same');

pred = nan(size(phi));
good = den > 1e-6;
pred(good) = num(good) ./ den(good);

end
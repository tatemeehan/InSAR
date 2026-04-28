function out = nan_gauss_smooth(in, sigma)

if nargin < 2 || isempty(sigma)
    sigma = 2;
end

valid = ~isnan(in);
x = in;
x(~valid) = 0;

num = imgaussfilt(x, sigma, 'FilterDomain', 'spatial');
den = imgaussfilt(double(valid), sigma, 'FilterDomain', 'spatial');

out = num ./ den;
out(den < eps) = NaN;

end
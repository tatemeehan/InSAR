function cmap = warpColormapPiecewise(ctrl, n, splitVal, fracLow)
% ctrl     : control RGB values, size m x 3
% n        : final number of colors
% splitVal : where to split the map in [0,1], e.g. 0.85
% fracLow  : fraction of colors assigned below splitVal, e.g. 0.5

if nargin < 2 || isempty(n), n = 256; end
if nargin < 3 || isempty(splitVal), splitVal = 0.85; end
if nargin < 4 || isempty(fracLow), fracLow = 0.5; end

nLow  = round(fracLow * n);
nHigh = n - nLow;

x = linspace(0,1,size(ctrl,1));

xiLow  = linspace(0, splitVal, nLow);
xiHigh = linspace(splitVal, 1, nHigh+1);
xiHigh = xiHigh(2:end);   % avoid duplicate split point

xi = [xiLow, xiHigh];

cmap = interp1(x, ctrl, xi, 'pchip');
cmap = max(min(cmap,1),0);
end
function cmap = terrainPermCmapCool(n)

if nargin < 1
    n = 256;
end

ctrl = [
    0.96 0.97 0.98
    0.84 0.88 0.91
    0.67 0.75 0.80
    0.53 0.60 0.63
    0.45 0.40 0.36
    0.20 0.19 0.18
];

x  = linspace(0,1,size(ctrl,1));
xi = linspace(0,1,n);

cmap = interp1(x, ctrl, xi, 'pchip');
cmap = max(min(cmap,1),0);
end
function cmap = terrainPermCmap(n)
%TERRAINPERMCMAP Custom muted terrain-inspired colormap for poster figures.
%
% Usage:
%   colormap(terrainPermCmap)
%   colormap(terrainPermCmap(256))

if nargin < 1
    n = 256;
end

% Control colors: light snow -> icy blue -> tan -> brown -> dark earth
ctrl = [
    0.95 0.96 0.97  % very light snow-gray
    0.82 0.87 0.90  % pale icy blue
    0.67 0.74 0.78  % blue-gray
    0.82 0.75 0.63  % light tan
    0.63 0.52 0.40  % brown
    0.30 0.25 0.22  % dark earth
];

x  = linspace(0,1,size(ctrl,1));
% xi = linspace(0,.85,n);xi = [xi,linspace(0.85,1,n)];
xi = linspace(0,1,n);


cmap = interp1(x, ctrl, xi, 'pchip');
cmap = max(min(cmap,1),0);
end
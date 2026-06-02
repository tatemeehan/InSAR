function lut = build_soil_VWC0_constraint_LUT(soilprofBase, f_Hz, theta_i_deg, Lg_default, varargin)
%BUILD_SOIL_VWC0_CONSTRAINT_LUT
% Build inverse lookup table from modeled real soil permittivity to VWC0.
%
% Purpose
%   Given a fixed 1-D soil profile parameterization, vary only VWC0 and
%   compute the modeled real soil permittivity target. This allows a
%   GPR/ML-derived near-surface eps' raster to be converted into an
%   effective VWC0 raster.
%
% Inputs
%   soilprofBase : baseline soil profile struct
%   f_Hz         : radar frequency [Hz]
%   theta_i_deg  : incidence angle [deg], scalar for LUT construction
%   Lg_default   : soil weighting length used by forward model
%
% Name-value options
%   'VWC0Grid'    : candidate VWC0 values, default linspace(0, SWC, 1000)
%   'Metric'      : 'surface', 'shallow', 'boxcar', or 'energy'
%   'AnchorDepth' : depth scale [m] for shallow/boxcar metrics, default 0.05
%   'Verbose'     : true/false
%
% Output
%   lut : structure containing forward and inverse LUTs

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;

p.addParameter('VWC0Grid', [], @isnumeric);
p.addParameter('Metric', 'surface', @(x) ischar(x) || isstring(x));
p.addParameter('AnchorDepth', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Verbose', true, @(x) islogical(x) || isnumeric(x));

p.parse(varargin{:});
opts = p.Results;

metric = lower(string(opts.Metric));

% -------------------------
% VWC0 search grid
% -------------------------
if isempty(opts.VWC0Grid)
    VWC0Grid = linspace(0, soilprofBase.SWC, 1000).';
else
    VWC0Grid = double(opts.VWC0Grid(:));
end

VWC0Grid = unique(VWC0Grid, 'sorted');

if any(VWC0Grid < 0) || any(VWC0Grid > soilprofBase.SWC)
    warning('VWC0Grid extends outside [0, soilprofBase.SWC]. Values will be clipped.');
    VWC0Grid = min(max(VWC0Grid, 0), soilprofBase.SWC);
    VWC0Grid = unique(VWC0Grid, 'sorted');
end

nV = numel(VWC0Grid);

% -------------------------
% Compute kx from incidence angle in air
% -------------------------
c = 299792458;
k0 = 2*pi*f_Hz/c;
kx = k0 * sind(theta_i_deg);

% -------------------------
% Allocate
% -------------------------
epsMetric  = nan(nV, 1);
epsSurface = nan(nV, 1);
IgGrid     = complex(nan(nV, 1), nan(nV, 1));

if opts.Verbose
    fprintf('Building soil VWC0 LUT: %d VWC0 nodes, metric = %s\n', ...
        nV, metric);
end

% -------------------------
% Evaluate profile model over VWC0
% -------------------------
for ii = 1:nV

    prof_i = soilprofBase;
    prof_i.enable = true;
    prof_i.VWC0 = VWC0Grid(ii);

    [Ig_i, diag_i] = soil_profile_integral( ...
        f_Hz, ...
        kx, ...
        prof_i, ...
        Lg_default);

    z = diag_i.z(:);
    e = real(diag_i.eps(:));

    epsSurface(ii) = e(1);
    IgGrid(ii) = Ig_i;

    switch metric

        case "surface"
            % Match the modeled z = 0 real permittivity.
            epsMetric(ii) = e(1);

        case "shallow"
            % Exponentially weighted shallow average.
            w = exp(-z ./ opts.AnchorDepth);
            epsMetric(ii) = trapz(z, e .* w) ./ trapz(z, w);

        case "boxcar"
            % Simple average over top AnchorDepth.
            use = z <= opts.AnchorDepth;

            if ~any(use)
                error('No soil nodes found within AnchorDepth = %.4f m.', opts.AnchorDepth);
            end

            epsMetric(ii) = mean(e(use), 'omitnan');

        case "energy"
            % Energy-weighted average using soil integral diagnostic.
            W = diag_i.W_energy(:);
            den = trapz(z, W);

            if den > 0
                epsMetric(ii) = trapz(z, e .* W) ./ den;
            else
                epsMetric(ii) = NaN;
            end

        otherwise
            error('Unknown Metric: %s', metric);
    end
end

% -------------------------
% Build inverse mapping epsReal -> VWC0
% -------------------------
good = isfinite(epsMetric) & isfinite(VWC0Grid);

epsGood = epsMetric(good);
vwcGood = VWC0Grid(good);

% Sort by eps for inverse interpolation
[epsSorted, idx] = sort(epsGood);
vwcSorted = vwcGood(idx);

% Remove duplicate eps values
[epsUnique, ia] = unique(epsSorted, 'stable');
vwcUnique = vwcSorted(ia);

% Diagnostics: monotonicity
dE = diff(epsMetric);
fracIncreasing = mean(dE > 0, 'omitnan');

if opts.Verbose
    fprintf('  epsMetric range: %.3f to %.3f\n', ...
        min(epsMetric, [], 'omitnan'), max(epsMetric, [], 'omitnan'));
    fprintf('  VWC0 range: %.3f to %.3f\n', min(VWC0Grid), max(VWC0Grid));
    fprintf('  Fraction increasing: %.3f\n', fracIncreasing);
end

if fracIncreasing < 0.95
    warning(['epsMetric(VWC0) is not strongly monotonic. ', ...
             'Inspect the LUT before trusting the inversion.']);
end

% -------------------------
% Package output
% -------------------------
lut = struct();

lut.VWC0Grid = VWC0Grid;
lut.epsMetric = epsMetric;
lut.epsSurface = epsSurface;
lut.IgGrid = IgGrid;

lut.epsUnique = epsUnique;
lut.VWC0Unique = vwcUnique;

lut.epsMin = min(epsUnique);
lut.epsMax = max(epsUnique);

lut.metric = metric;
lut.anchorDepth = opts.AnchorDepth;
lut.theta_i_deg = theta_i_deg;
lut.f_Hz = f_Hz;
lut.kx = kx;
lut.Lg_default = Lg_default;
lut.SWC = soilprofBase.SWC;
lut.soilprofBase = soilprofBase;

lut.fracIncreasing = fracIncreasing;

end
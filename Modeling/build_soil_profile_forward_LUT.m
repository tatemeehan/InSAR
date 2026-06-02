function soilFwdLUT = build_soil_profile_forward_LUT(soilprofBase, f_Hz, thetaGrid_deg, VWC0Grid, Lg_default)
%BUILD_SOIL_PROFILE_FORWARD_LUT Build LUT for soil profile response.
%
% Inputs
%   soilprofBase  : fixed soil profile struct
%   f_Hz          : radar frequency [Hz]
%   thetaGrid_deg : incidence-angle grid [deg]
%   VWC0Grid      : effective surface VWC grid
%   Lg_default    : soil volume weighting length
%
% Outputs
%   soilFwdLUT.Ig        : complex Ig(VWC0, theta)
%   soilFwdLUT.epsSurf   : complex eps(z=0), mostly VWC0-dependent
%   soilFwdLUT.FIgReal   : gridded interpolant
%   soilFwdLUT.FIgImag   : gridded interpolant
%   soilFwdLUT.FepsReal  : 1-D interpolant
%   soilFwdLUT.FepsImag  : 1-D interpolant

VWC0Grid = double(VWC0Grid(:));
thetaGrid_deg = double(thetaGrid_deg(:));

nV = numel(VWC0Grid);
nT = numel(thetaGrid_deg);

IgTable = complex(nan(nV, nT), nan(nV, nT));
epsSurfTable = complex(nan(nV, 1), nan(nV, 1));

c = 299792458;
k0 = 2*pi*f_Hz/c;

fprintf('Building soil profile forward LUT: %d VWC0 x %d theta = %d profile calls\n', ...
    nV, nT, nV*nT);

for iV = 1:nV

    prof_i = soilprofBase;
    prof_i.enable = true;
    prof_i.VWC0 = VWC0Grid(iV);

    for iT = 1:nT

        theta_i = thetaGrid_deg(iT);
        kx = k0 * sind(theta_i);

        [Ig_i, diag_i] = soil_profile_integral( ...
            f_Hz, ...
            kx, ...
            prof_i, ...
            Lg_default);

        IgTable(iV, iT) = Ig_i;

        if iT == 1
            epsSurfTable(iV) = diag_i.eps(1);
        end
    end
end

soilFwdLUT = struct();

soilFwdLUT.VWC0Grid = VWC0Grid;
soilFwdLUT.thetaGrid_deg = thetaGrid_deg;
soilFwdLUT.f_Hz = f_Hz;
soilFwdLUT.Lg_default = Lg_default;

soilFwdLUT.IgTable = IgTable;
soilFwdLUT.epsSurfTable = epsSurfTable;

soilFwdLUT.FIgReal = griddedInterpolant( ...
    {VWC0Grid, thetaGrid_deg}, ...
    real(IgTable), ...
    'linear', ...
    'nearest');

soilFwdLUT.FIgImag = griddedInterpolant( ...
    {VWC0Grid, thetaGrid_deg}, ...
    imag(IgTable), ...
    'linear', ...
    'nearest');

soilFwdLUT.FepsReal = griddedInterpolant( ...
    {VWC0Grid}, ...
    real(epsSurfTable), ...
    'linear', ...
    'nearest');

soilFwdLUT.FepsImag = griddedInterpolant( ...
    {VWC0Grid}, ...
    imag(epsSurfTable), ...
    'linear', ...
    'nearest');

soilFwdLUT.VWC0Min = min(VWC0Grid);
soilFwdLUT.VWC0Max = max(VWC0Grid);
soilFwdLUT.thetaMin = min(thetaGrid_deg);
soilFwdLUT.thetaMax = max(thetaGrid_deg);

end

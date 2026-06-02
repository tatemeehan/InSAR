function [IgRaster, epsSurfRaster] = eval_soil_profile_forward_LUT( ...
    VWC0Raster, thetaRaster_deg, soilFwdLUT, validMask, blockSize)
%EVAL_SOIL_PROFILE_FORWARD_LUT Evaluate soil profile LUT on raster grids.
%
% Inputs
%   VWC0Raster      : effective VWC0 raster
%   thetaRaster_deg : local incidence angle raster [deg]
%   soilFwdLUT      : output from build_soil_profile_forward_LUT
%   validMask       : optional logical valid mask
%   blockSize       : number of valid pixels per block, default 1e6
%
% Outputs
%   IgRaster        : complex soil profile integral raster
%   epsSurfRaster   : complex soil surface permittivity raster, eps(z=0)

if nargin < 4 || isempty(validMask)
    validMask = isfinite(VWC0Raster) & isfinite(thetaRaster_deg);
end

if nargin < 5 || isempty(blockSize)
    blockSize = 1e6;
end

if ~isequal(size(VWC0Raster), size(thetaRaster_deg))
    error('VWC0Raster and thetaRaster_deg must have the same size.');
end

valid = validMask & isfinite(VWC0Raster) & isfinite(thetaRaster_deg);

IgRaster = complex( ...
    nan(size(VWC0Raster), 'single'), ...
    nan(size(VWC0Raster), 'single'));

epsSurfRaster = complex( ...
    nan(size(VWC0Raster), 'single'), ...
    nan(size(VWC0Raster), 'single'));

idx = find(valid);
n = numel(idx);
nBlocks = ceil(n / blockSize);

fprintf('Evaluating soil forward LUT for %d pixels in %d blocks\n', ...
    n, nBlocks);

for ib = 1:nBlocks

    i1 = (ib-1)*blockSize + 1;
    i2 = min(ib*blockSize, n);

    idb = idx(i1:i2);

    VWC0Eval = double(VWC0Raster(idb));
    thetaEval = double(thetaRaster_deg(idb));

    % Clip to LUT bounds
    VWC0Eval = min(max(VWC0Eval, soilFwdLUT.VWC0Min), soilFwdLUT.VWC0Max);
    thetaEval = min(max(thetaEval, soilFwdLUT.thetaMin), soilFwdLUT.thetaMax);

    IgReal = soilFwdLUT.FIgReal(VWC0Eval, thetaEval);
    IgImag = soilFwdLUT.FIgImag(VWC0Eval, thetaEval);

    epsReal = soilFwdLUT.FepsReal(VWC0Eval);
    epsImag = soilFwdLUT.FepsImag(VWC0Eval);

    IgRaster(idb) = complex(single(IgReal), single(IgImag));
    epsSurfRaster(idb) = complex(single(epsReal), single(epsImag));

    fprintf('  Soil LUT block %d / %d complete\n', ib, nBlocks);
end

end
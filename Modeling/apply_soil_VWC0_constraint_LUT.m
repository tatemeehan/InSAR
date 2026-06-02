function [VWC0Raster, diag] = apply_soil_VWC0_constraint_LUT(epsRealRaster, lut, validMask, varargin)
%APPLY_SOIL_VWC0_CONSTRAINT_LUT Convert real-permittivity raster to VWC0 raster.
%
% Inputs
%   epsRealRaster : GPR/ML-derived real soil permittivity raster
%   lut           : output from build_soil_VWC0_constraint_LUT
%   validMask     : optional valid mask
%
% Name-value options
%   'OutputClass' : 'double' or 'single', default 'single'
%   'ClipInput'   : true/false, default true
%
% Outputs
%   VWC0Raster : effective VWC0 raster
%   diag       : diagnostics

p = inputParser;
p.addParameter('OutputClass', 'single', @(x) ischar(x) || isstring(x));
p.addParameter('ClipInput', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

opts = p.Results;
outputClass = lower(string(opts.OutputClass));

if nargin < 3 || isempty(validMask)
    validMask = isfinite(epsRealRaster);
end

if ~isequal(size(validMask), size(epsRealRaster))
    error('validMask must be the same size as epsRealRaster.');
end

valid = validMask & isfinite(epsRealRaster);

VWC0Raster = nan(size(epsRealRaster));

epsEval = double(epsRealRaster(valid));

nBelow = nnz(epsEval < lut.epsMin);
nAbove = nnz(epsEval > lut.epsMax);

if logical(opts.ClipInput)
    epsEval = min(max(epsEval, lut.epsMin), lut.epsMax);
else
    bad = epsEval < lut.epsMin | epsEval > lut.epsMax;
    epsEval(bad) = NaN;
end

% Inverse LUT: epsReal -> VWC0
VWC0Eval = interp1( ...
    lut.epsUnique, ...
    lut.VWC0Unique, ...
    epsEval, ...
    'linear', ...
    NaN);

% Physical clipping
VWC0Eval = min(max(VWC0Eval, 0), lut.SWC);

VWC0Raster(valid) = VWC0Eval;

switch outputClass
    case "single"
        VWC0Raster = single(VWC0Raster);
    case "double"
        % leave as double
    otherwise
        error('OutputClass must be "single" or "double".');
end

diag = struct();
diag.nValid = nnz(valid);
diag.nBelowLUT = nBelow;
diag.nAboveLUT = nAbove;
diag.epsMinLUT = lut.epsMin;
diag.epsMaxLUT = lut.epsMax;
diag.VWC0Min = min(VWC0Eval, [], 'omitnan');
diag.VWC0Max = max(VWC0Eval, [], 'omitnan');
diag.VWC0Mean = mean(VWC0Eval, 'omitnan');
diag.VWC0Median = median(VWC0Eval, 'omitnan');

end
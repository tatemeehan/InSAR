function val = validate_soil_VWC0_constraint(epsRealRaster, VWC0Raster, lut, validMask)
%VALIDATE_SOIL_VWC0_CONSTRAINT Check LUT inversion consistency.
%
% This does not rerun the full soil profile integral per pixel. It checks
% the inverse LUT consistency:
%
%   epsRealRaster -> VWC0Raster -> epsRecon

if nargin < 4 || isempty(validMask)
    validMask = isfinite(epsRealRaster) & isfinite(VWC0Raster);
end

valid = validMask & isfinite(epsRealRaster) & isfinite(VWC0Raster);

epsTarget = double(epsRealRaster(valid));
VWC0Eval  = double(VWC0Raster(valid));

epsRecon = interp1( ...
    lut.VWC0Grid, ...
    lut.epsMetric, ...
    VWC0Eval, ...
    'linear', ...
    NaN);

err = epsRecon - epsTarget;

val = struct();
val.nValid = nnz(valid);
val.errMean = mean(err, 'omitnan');
val.errMedian = median(err, 'omitnan');
val.errRMSE = sqrt(mean(err.^2, 'omitnan'));
val.errMAE = mean(abs(err), 'omitnan');
val.errP05 = prctile(err, 5);
val.errP95 = prctile(err, 95);

fprintf('\nSoil VWC0 LUT inversion validation\n');
fprintf('  nValid     = %d\n', val.nValid);
fprintf('  mean err   = %.4f eps units\n', val.errMean);
fprintf('  median err = %.4f eps units\n', val.errMedian);
fprintf('  RMSE       = %.4f eps units\n', val.errRMSE);
fprintf('  MAE        = %.4f eps units\n', val.errMAE);
fprintf('  P05/P95    = %.4f / %.4f eps units\n', val.errP05, val.errP95);

end
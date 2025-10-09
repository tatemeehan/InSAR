function [amp, pow, db, meta] = compute_sli(slc, P, a, lambda, pixelArea, incidence, slope, slantRange, correction)
%COMPUTE_SLI Compute Single Look Geocoded Image with optional RCS, terrain, and range correction.
%
%   [amp, pow, db, meta] = compute_sli(slc, P, a, lambda, pixelArea, incidence, slope, correction, slantRange)
%
%   Inputs:
%       slc         - complex Single Look Complex image (matrix)
%       P           - Returned power from reference target (linear power)
%       a           - Side length of trihedral CR (m), default = 1
%       lambda      - Radar wavelength (m), default = 0.3/1.3 (L-band)
%       pixelArea   - (optional) pixel area raster (m^2), for area normalization
%       incidence   - (optional) incidence angle raster (degrees)
%       slope       - (optional) terrain slope raster (degrees)
%       slantRange  - (optional) slant range raster (m) for range correction (R^2)
%       correction  - (optional) string or cell array: {'gamma0','beta0','pixelarea','rangecorrect'}
%
%   Outputs:
%       amp   - Backscatter amplitude image
%       pow   - Calibrated backscatter intensity (linear)
%       db    - Backscatter intensity (dB scale)
%       meta  - Struct containing applied corrections and scale factors

% Defaults
if nargin < 3 || isempty(a),       a = 1; end
if nargin < 4 || isempty(lambda), lambda = 0.3/1.3; end
if nargin < 5, pixelArea = []; end
if nargin < 6, incidence = []; end
if nargin < 7, slope = []; end
if nargin < 8 || isempty(correction), correction = {}; end
if ischar(correction), correction = {correction}; end
if nargin < 9, slantRange = []; end

% Compute raw amplitude and power
amp = abs(slc);
pow = amp.^2;

meta = struct();
meta.corrections = correction;

% ----- RCS Calibration (if P is given) -----
if ~isempty(P)
    RCS = (4 * pi * a^4) / (3 * lambda^2); % Trihedral RCS
    scale = RCS / P;
    pow = pow * scale;
    meta.rcsScale = scale;
    meta.rcs = RCS;
else
    meta.rcsScale = 1;
    meta.rcs = NaN;
end

% ----- Sequential Corrections -----
for i = 1:length(correction)
    step = lower(correction{i});
    switch step
        case 'pixelarea'
            if isempty(pixelArea)
                warning('Pixel area raster not provided. Skipping pixelarea normalization.');
                meta.pixelAreaNorm = false;
            else
                pow = pow ./ pixelArea;
                meta.pixelAreaNorm = true;
            end

        case 'gamma0'
            if isempty(incidence)
                error('Incidence angle required for gamma0 correction.');
            end
            pow = pow ./ cosd(incidence);
            meta.terrainCorrectionGamma0 = true;

        case 'beta0'
            if isempty(slope)
                error('Slope required for beta0 correction.');
            end
            pow = pow ./ cosd(slope);
            meta.terrainCorrectionBeta0 = true;

        case {'rangecorrect', 'range2'}
            if isempty(slantRange)
                warning('Slant range raster not provided. Skipping range correction.');
                meta.rangeCorrection = false;
            else
                pow = pow .* (slantRange.^2);  % Standard R^2 correction
                meta.rangeCorrection = 'R^2';
            end

        case 'range4'
            if isempty(slantRange)
                warning('Slant range raster not provided. Skipping R^4 correction.');
                meta.rangeCorrection = false;
            else
                pow = pow .* (slantRange.^4);  % Full R^4 correction
                meta.rangeCorrection = 'R^4';
            end

        otherwise
            warning('Unknown correction: %s. Skipped.', step);
    end
end

% Convert to dB
% db = 10 * log10(max(pow, eps));
% db = 10 * log10(pow);
% db(~isfinite(pow)) = NaN;
db = nan(size(pow));
valid = pow > 0 & isfinite(pow);
db(valid) = 10 * log10(pow(valid));
end
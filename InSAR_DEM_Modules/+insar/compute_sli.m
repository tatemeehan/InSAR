function [amp, pow, db, meta] = compute_sli(slc, P, a, lambda, pixelArea, incidence, slope, slantRange, correction)
%COMPUTE_SLI Compute Single Look Geocoded Image with optional RCS, terrain, and range correction.

% Defaults
if nargin < 3 || isempty(a),       a = 1; end
if nargin < 4 || isempty(lambda),  lambda = 0.3/1.3; end
if nargin < 5, pixelArea = []; end
if nargin < 6, incidence = []; end
if nargin < 7, slope = []; end
if nargin < 8, slantRange = []; end
if nargin < 9 || isempty(correction), correction = {}; end
if ischar(correction), correction = {correction}; end

% Work in single for raster products
if ~isa(slc, 'single')
    slc = single(slc);
end
if ~isempty(pixelArea) && ~isa(pixelArea, 'single'),   pixelArea = single(pixelArea); end
if ~isempty(incidence) && ~isa(incidence, 'single'),   incidence = single(incidence); end
if ~isempty(slope) && ~isa(slope, 'single'),           slope = single(slope); end
if ~isempty(slantRange) && ~isa(slantRange, 'single'), slantRange = single(slantRange); end

% Compute raw amplitude and power
amp = abs(slc);       % single
pow = amp .* amp;     % single

meta = struct();
meta.corrections = correction;

% ----- RCS Calibration -----
if ~isempty(P)
    % Scalars can stay double; cast final scale for raster math
    RCS = (4 * pi * a^4) / (3 * lambda^2);
    scale = single(RCS / P);
    pow = pow .* scale;

    meta.rcsScale = double(scale);
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
                pow = pow .* (slantRange .* slantRange);
                meta.rangeCorrection = 'R^2';
            end

        case 'range4'
            if isempty(slantRange)
                warning('Slant range raster not provided. Skipping R^4 correction.');
                meta.rangeCorrection = false;
            else
                sr2 = slantRange .* slantRange;
                pow = pow .* (sr2 .* sr2);
                meta.rangeCorrection = 'R^4';
            end

        otherwise
            warning('Unknown correction: %s. Skipped.', step);
    end
end

% Convert to dB as single
db = nan(size(pow), 'single');
valid = (pow > 0) & isfinite(pow);
db(valid) = single(10) .* log10(pow(valid));

% Optional: ensure outputs are single even after any promotion
amp = single(amp);
pow = single(pow);
db  = single(db);
end
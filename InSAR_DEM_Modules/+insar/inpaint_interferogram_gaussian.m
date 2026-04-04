function out = inpaint_interferogram_gaussian(insarPair, opts)
%INPAINT_INTERFEROGRAM_GAUSSIAN
% Gaussian inpainting for referenced phase and optional coherence using a
% shared trusted/support mask definition.
%
% INPUTS
%   insarPair : struct with fields
%       .phzReferenced   referenced phase raster
%       .cor             coherence raster
%
%   opts : options struct
%       .trainCohThresh      default 0.75
%       .crSupportMask       logical mask
%       .directSolvedMask    logical mask
%       .maxTrainDistToCR    default inf
%       .sameWrapMask        optional logical mask
%       .validMask           optional logical mask
%
%       .phaseSigma          default 2
%       .cohSigma            default 2
%       .fillCoherence       default true
%
%       .qualitySigma        default 5
%       .qualityMinValidFrac default 0.35
%
% OUTPUT
%   out : struct with fields
%       .trustedMask
%       .crSupportMask
%       .directSolvedMask
%       .distToTrusted
%       .distToCRSupport
%       .distToDirectSolved
%       .phzFilled
%       .cohFilled
%       .phaseQuality
%       .qualitySupportFrac
%       .support

if nargin < 2
    opts = struct();
end

if ~isfield(opts, 'phaseSigma') || isempty(opts.phaseSigma)
    opts.phaseSigma = 2;
end
if ~isfield(opts, 'cohSigma') || isempty(opts.cohSigma)
    opts.cohSigma = opts.phaseSigma;
end
if ~isfield(opts, 'fillCoherence') || isempty(opts.fillCoherence)
    opts.fillCoherence = true;
end
if ~isfield(opts, 'qualitySigma') || isempty(opts.qualitySigma)
    opts.qualitySigma = 5;
end
if ~isfield(opts, 'qualityMinValidFrac') || isempty(opts.qualityMinValidFrac)
    opts.qualityMinValidFrac = 0.35;
end

phz = double(insarPair.phzReferenced);
cor = double(insarPair.cor);

% -------------------------------------------------------------------------
% Shared support logic
% -------------------------------------------------------------------------
support = insar.build_phase_fill_support_masks(phz, cor, opts);
trustedMask = support.trustedMask;

% -------------------------------------------------------------------------
% Phase interpolation
% -------------------------------------------------------------------------
phzFilled = insar.phase_fill_gaussian_interp(phz, trustedMask, opts.phaseSigma);

% Preserve trusted observed phase exactly
phzFilled(trustedMask) = phz(trustedMask);

% -------------------------------------------------------------------------
% Coherence interpolation (optional)
% -------------------------------------------------------------------------
if opts.fillCoherence
    cohFilled = insar.gaussian_interp_coherence(cor, trustedMask, opts.cohSigma);

    % Preserve trusted observed coherence exactly
    cohFilled(trustedMask) = cor(trustedMask);
else
    cohFilled = [];
end

% -------------------------------------------------------------------------
% Phase-quality map from restored phase
% -------------------------------------------------------------------------
validMaskForQuality = support.validMask & isfinite(phzFilled);
[phaseQuality, qualitySupportFrac] = insar.phase_quality_gaussian( ...
    phzFilled, opts.qualitySigma, opts.qualityMinValidFrac,validMaskForQuality);

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
out = struct();
out.trustedMask = trustedMask;
out.crSupportMask = support.crSupportMask;
out.directSolvedMask = support.directSolvedMask;
out.distToTrusted = support.distToTrusted;
out.distToCRSupport = support.distToCRSupport;
out.distToDirectSolved = support.distToDirectSolved;
out.support = support;

out.phzFilled = phzFilled;
out.cohFilled = cohFilled;
out.phaseQuality = phaseQuality;
out.qualitySupportFrac = qualitySupportFrac;

end
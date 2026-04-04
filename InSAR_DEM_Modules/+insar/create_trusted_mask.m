function support = create_trusted_mask(phzReferenced, cor, opts)
%CREATE_TRUSTED_MASK
% Build shared trusted/support masks for phase/coherence interpolation.
%
% Inputs
%   phzReferenced : referenced phase raster
%   cor           : coherence raster
%   opts          : struct
%
% Optional opts fields
%   .trainCohThresh    default 0.75
%   .crSupportMask     logical raster, default false(size(phzReferenced))
%   .directSolvedMask  logical raster, default true(size(phzReferenced))
%   .maxTrainDistToCR  default inf
%   .sameWrapMask      optional logical raster
%   .validMask         optional logical raster
%
% Output
%   support : struct with
%       .trustedMask
%       .validMask
%       .crSupportMask
%       .directSolvedMask
%       .distToTrusted
%       .distToCRSupport
%       .distToDirectSolved

if nargin < 3
    opts = struct();
end

if ~isfield(opts, 'trainCohThresh') || isempty(opts.trainCohThresh)
    opts.trainCohThresh = 0.75;
end
if ~isfield(opts, 'crSupportMask') || isempty(opts.crSupportMask)
    opts.crSupportMask = false(size(phzReferenced));
end
if ~isfield(opts, 'directSolvedMask') || isempty(opts.directSolvedMask)
    opts.directSolvedMask = true(size(phzReferenced));
end
if ~isfield(opts, 'maxTrainDistToCR') || isempty(opts.maxTrainDistToCR)
    opts.maxTrainDistToCR = inf;
end

phz = double(phzReferenced);
cor = double(cor);
[nr, nc] = size(phz);

validMask = isfinite(phz) & isfinite(cor);
if isfield(opts, 'validMask') && ~isempty(opts.validMask) && isequal(size(opts.validMask), [nr nc])
    validMask = validMask & logical(opts.validMask);
end

crSupportMask = logical(opts.crSupportMask);
if ~isequal(size(crSupportMask), [nr nc])
    crSupportMask = false(nr, nc);
end

directSolvedMask = logical(opts.directSolvedMask);
if ~isequal(size(directSolvedMask), [nr nc])
    directSolvedMask = true(nr, nc);
end

baseTrustedSeed = validMask & (cor >= opts.trainCohThresh);

distToTrusted = bwdist(baseTrustedSeed);
distToCRSupport = bwdist(crSupportMask);
distToDirectSolved = bwdist(directSolvedMask & validMask);

useCRGate = any(crSupportMask(:)) && isfinite(opts.maxTrainDistToCR);
if useCRGate
    crGate = distToCRSupport <= opts.maxTrainDistToCR;
else
    crGate = true(nr, nc);
end

trustedMask = validMask & ...
              (cor >= opts.trainCohThresh) & ...
              directSolvedMask & ...
              crGate;

if isfield(opts, 'sameWrapMask') && ~isempty(opts.sameWrapMask) && isequal(size(opts.sameWrapMask), [nr nc])
    trustedMask = trustedMask & logical(opts.sameWrapMask);
end

support = struct();
support.trustedMask = trustedMask;
support.validMask = validMask;
support.crSupportMask = crSupportMask;
support.directSolvedMask = directSolvedMask;
support.distToTrusted = distToTrusted;
support.distToCRSupport = distToCRSupport;
support.distToDirectSolved = distToDirectSolved;
support.baseTrustedSeed = baseTrustedSeed;

end
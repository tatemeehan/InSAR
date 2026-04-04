function support = build_phase_fill_support_masks(phzReferenced, cor, opts)
%BUILD_PHASE_FILL_SUPPORT_MASKS
% Centralized trusted/support mask builder for phase/coherence interpolation.

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

distToTrustedSeed  = bwdist(baseTrustedSeed);
distToCRSupport    = bwdist(crSupportMask);
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
support.validMask = validMask;
support.baseTrustedSeed = baseTrustedSeed;
support.trustedMask = trustedMask;
support.crSupportMask = crSupportMask;
support.directSolvedMask = directSolvedMask;
support.distToTrusted = distToTrustedSeed;
support.distToCRSupport = distToCRSupport;
support.distToDirectSolved = distToDirectSolved;
support.useCRGate = useCRGate;
support.trainCohThresh = opts.trainCohThresh;

end
function holdoutMask = make_random_holdout_mask(trustedMask, frac)
%MAKE_RANDOM_HOLDOUT_MASK Randomly hold out a fraction of trusted pixels.

if nargin < 2 || isempty(frac)
    frac = 0.15;
end

idx = find(trustedMask);
n = numel(idx);
nHold = max(1, round(frac * n));

rng(1);
perm = randperm(n, nHold);

holdoutMask = false(size(trustedMask));
holdoutMask(idx(perm)) = true;

end
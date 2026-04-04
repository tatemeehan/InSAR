function holdoutMask = make_block_holdout_mask(trustedMask, blockRadius, numBlocks)
%MAKE_BLOCK_HOLDOUT_MASK Hold out spatial blocks from the trusted region.

if nargin < 2 || isempty(blockRadius)
    blockRadius = 20;
end
if nargin < 3 || isempty(numBlocks)
    numBlocks = 25;
end

[nr, nc] = size(trustedMask);
holdoutMask = false(nr, nc);

[idxR, idxC] = find(trustedMask);
nCand = numel(idxR);

if nCand == 0
    return;
end

rng(1);
perm = randperm(nCand, min(numBlocks, nCand));

[RR, CC] = ndgrid(1:nr, 1:nc);

for ii = 1:numel(perm)
    r0 = idxR(perm(ii));
    c0 = idxC(perm(ii));

    patch = hypot(RR - r0, CC - c0) <= blockRadius;
    holdoutMask = holdoutMask | (patch & trustedMask);
end

end
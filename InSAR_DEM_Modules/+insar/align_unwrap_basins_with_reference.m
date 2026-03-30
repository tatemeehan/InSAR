function [phzAligned, basinMeta] = align_unwrap_basins_with_reference( ...
    phzWrapped, phzUnwrapped, cor, owner, crResult, pred, i, b1, j, b2, unwrapOpts)
%ALIGN_UNWRAP_BASINS_WITH_REFERENCE
% Resolve basin-to-basin integer 2*pi offsets after local multiseed unwrap.
%
% Strategy:
%   1) Build direct seam constraints from neighboring pixels where basins touch.
%   2) Build soft screen-based constraints across separated basins using a
%      smooth reference-like phase model.
%   3) Choose an anchor basin, preferably overlapping CR support.
%   4) Solve integer basin shifts and apply them.
%
% INPUTS
%   phzWrapped    - wrapped phase
%   phzUnwrapped  - locally unwrapped phase
%   cor           - coherence / quality
%   owner         - basin owner labels from region_grow_unwrap_bucket
%   crResult      - CR results struct used by pipeline
%   pred          - predictor struct from build_reference_predictors
%   i,b1,j,b2     - pair indices (kept for compatibility/debug)
%   unwrapOpts    - options
%
% Fields used from unwrapOpts
%   .basinMinEdgeSupport     default 5
%   .reconcileGapBasins      default true
%   .basinGapMaxPixels       default 25
%   .basinUseReferenceScreen default true
%   .basinScreenCohThresh    default max(0.4, qualityThresh)
%   .basinAnchorPolicy       default 'crFirst'
%
% OUTPUTS
%   phzAligned - globally aligned phase
%   basinMeta  - metadata

if nargin < 11
    unwrapOpts = struct();
end

if ~isfield(unwrapOpts,'basinMinEdgeSupport') || isempty(unwrapOpts.basinMinEdgeSupport)
    unwrapOpts.basinMinEdgeSupport = 5;
end
if ~isfield(unwrapOpts,'reconcileGapBasins') || isempty(unwrapOpts.reconcileGapBasins)
    unwrapOpts.reconcileGapBasins = true;
end
if ~isfield(unwrapOpts,'basinGapMaxPixels') || isempty(unwrapOpts.basinGapMaxPixels)
    unwrapOpts.basinGapMaxPixels = 25;
end
if ~isfield(unwrapOpts,'basinUseReferenceScreen') || isempty(unwrapOpts.basinUseReferenceScreen)
    unwrapOpts.basinUseReferenceScreen = true;
end
if ~isfield(unwrapOpts,'basinScreenCohThresh') || isempty(unwrapOpts.basinScreenCohThresh)
    if isfield(unwrapOpts,'qualityThresh') && ~isempty(unwrapOpts.qualityThresh)
        unwrapOpts.basinScreenCohThresh = max(0.4, unwrapOpts.qualityThresh);
    else
        unwrapOpts.basinScreenCohThresh = 0.4;
    end
end
if ~isfield(unwrapOpts,'basinAnchorPolicy') || isempty(unwrapOpts.basinAnchorPolicy)
    unwrapOpts.basinAnchorPolicy = 'crFirst';
end

if ~isfield(unwrapOpts,'basinGapMaxNeighbors') || isempty(unwrapOpts.basinGapMaxNeighbors)
    unwrapOpts.basinGapMaxNeighbors = 3;
end
if ~isfield(unwrapOpts,'basinGapMaxPasses') || isempty(unwrapOpts.basinGapMaxPasses)
    unwrapOpts.basinGapMaxPasses = 3;
end
if ~isfield(unwrapOpts,'basinMinGapArea') || isempty(unwrapOpts.basinMinGapArea)
    unwrapOpts.basinMinGapArea = 20;
end
phzAligned = phzUnwrapped;

numBasins = double(max(owner(:)));
if numBasins <= 1
    basinMeta = struct( ...
        'shift', 0, ...
        'solved', true, ...
        'anchorBasin', 1, ...
        'screenModel', [], ...
        'directEdges', struct('a',[],'b',[],'delta',[],'support',[],'weight',[]), ...
        'gapEdges', struct('a',[],'b',[],'delta',[],'support',[],'weight',[]), ...
        'numEdges', 0);
    return;
end

% % directEdges = build_direct_edges(owner, phzWrapped, phzUnwrapped, cor, unwrapOpts.basinMinEdgeSupport);
% % 
% % screenModel = [];
% % gapEdges = struct('a',[],'b',[],'delta',[],'support',[],'weight',[]);
% % 
% % if unwrapOpts.reconcileGapBasins && unwrapOpts.basinUseReferenceScreen
% %     screenModel = fit_basin_reference_screen(phzUnwrapped, cor, owner, pred, unwrapOpts);
% %     gapEdges = build_gap_edges(owner, phzUnwrapped, cor, screenModel, unwrapOpts);
% % end
% % 
% % allEdges = concat_edges(directEdges, gapEdges);
% % 
% % anchorBasin = choose_anchor_basin(owner, cor, crResult, unwrapOpts);
% % 
% % [shift, solved] = solve_basin_graph(numBasins, allEdges, anchorBasin);
% 
% directEdges = build_direct_edges(owner, phzWrapped, phzUnwrapped, cor, unwrapOpts.basinMinEdgeSupport);
% 
% anchorBasin = choose_anchor_basin(owner, cor, crResult, unwrapOpts);
% 
% % First pass: direct seam constraints only
% [shift, solved] = solve_basin_graph(numBasins, directEdges, anchorBasin);
% 
% screenModel = [];
% gapEdges = struct('a',[],'b',[],'delta',[],'support',[],'weight',[]);
% 
% % Second pass: only if some basins remain unsolved
% if unwrapOpts.reconcileGapBasins && unwrapOpts.basinUseReferenceScreen && any(~solved)
%     basinStats = compute_basin_stats(owner, phzUnwrapped, cor, unwrapOpts);
% 
%     screenModel = fit_basin_reference_screen(phzUnwrapped, cor, owner, pred, unwrapOpts);
% 
%     gapEdges = build_gap_edges_sparse( ...
%         owner, phzUnwrapped, cor, screenModel, basinStats, solved, shift, unwrapOpts);
% 
%     allEdges = concat_edges(directEdges, gapEdges);
%     [shift, solved] = solve_basin_graph(numBasins, allEdges, anchorBasin);
% else
%     allEdges = directEdges;
% end
% 
% for b = 1:numBasins
%     if solved(b)
%         phzAligned(owner == b) = phzAligned(owner == b) + 2*pi*shift(b);
%     end
% end
% 
% basinMeta = struct();
% basinMeta.shift = shift;
% basinMeta.solved = solved;
% basinMeta.anchorBasin = anchorBasin;
% basinMeta.screenModel = screenModel;
% basinMeta.directEdges = directEdges;
% basinMeta.gapEdges = gapEdges;
% basinMeta.numEdges = numel(allEdges.a);
% basinMeta.pair = [i j];
% basinMeta.burst = [b1 b2];
% 
% end

directEdges = build_direct_edges(owner, phzWrapped, phzUnwrapped, cor, unwrapOpts.basinMinEdgeSupport);

anchorBasin = choose_anchor_basin(owner, cor, crResult, unwrapOpts);

% Pass 0: direct seam solve only
[shift, solved] = solve_basin_graph(numBasins, directEdges, anchorBasin);

screenModel = [];
gapEdges = struct('a',[],'b',[],'delta',[],'support',[],'weight',[], ...
                  'pass',[]);

allEdges = directEdges;
nSolvedPerPass = nnz(solved);

if unwrapOpts.reconcileGapBasins && unwrapOpts.basinUseReferenceScreen && any(~solved)

    basinStats = compute_basin_stats(owner, phzUnwrapped, cor, unwrapOpts);
    screenModel = fit_basin_reference_screen(phzUnwrapped, cor, owner, pred, unwrapOpts);

    maxPasses = unwrapOpts.basinGapMaxPasses;

    for pass = 1:maxPasses
        solvedPrev = solved;

        gapEdgesPass = build_gap_edges_sparse( ...
            owner, phzUnwrapped, cor, screenModel, basinStats, solved, shift, unwrapOpts, pass);

        if isempty(gapEdgesPass.a)
            break;
        end

        allEdges = concat_edges(allEdges, gapEdgesPass);
        gapEdges = concat_edges(gapEdges, gapEdgesPass);

        [shift, solved] = solve_basin_graph(numBasins, allEdges, anchorBasin);
        nSolvedPerPass(end+1,1) = nnz(solved); %#ok<AGROW>

        if isequal(solved, solvedPrev)
            break;
        end

        if ~any(~solved)
            break;
        end
    end
end

for b = 1:numBasins
    if solved(b)
        phzAligned(owner == b) = phzAligned(owner == b) + 2*pi*shift(b);
    end
end

basinMeta = struct();
basinMeta.shift = shift;
basinMeta.solved = solved;
basinMeta.anchorBasin = anchorBasin;
basinMeta.screenModel = screenModel;
basinMeta.directEdges = directEdges;
basinMeta.gapEdges = gapEdges;
basinMeta.numEdges = numel(allEdges.a);
basinMeta.nSolvedPerPass = nSolvedPerPass;
basinMeta.pair = [i j];
basinMeta.burst = [b1 b2];
basinMeta.nSolvedFinal = nnz(solved);
basinMeta.nBasinsTotal = numBasins;
end

% -------------------------------------------------------------------------
function edges = build_direct_edges(owner, phzWrapped, phzUnwrapped, cor, minSupport)

pairA = [];
pairB = [];
pairK = [];
pairW = [];

% vertical neighbors
oa = owner(1:end-1,:);
ob = owner(2:end,:);
mask = oa > 0 & ob > 0 & oa ~= ob & ...
       isfinite(phzWrapped(1:end-1,:)) & isfinite(phzWrapped(2:end,:)) & ...
       isfinite(phzUnwrapped(1:end-1,:)) & isfinite(phzUnwrapped(2:end,:));

if any(mask(:))
    wa = phzWrapped(1:end-1,:);
    wb = phzWrapped(2:end,:);
    ua = phzUnwrapped(1:end-1,:);
    ub = phzUnwrapped(2:end,:);
    qa = cor(1:end-1,:);
    qb = cor(2:end,:);

    dphi = wrap_to_pi_local(wb(mask) - wa(mask));
    kest = round((ua(mask) + dphi - ub(mask)) ./ (2*pi));
    wgt  = 0.5 * (double(qa(mask)) + double(qb(mask)));

    a = double(oa(mask));
    b = double(ob(mask));
    [a, b, kest] = canonicalize_pairs(a, b, kest);

    pairA = [pairA; a];
    pairB = [pairB; b];
    pairK = [pairK; kest];
    pairW = [pairW; wgt];
end

% horizontal neighbors
oa = owner(:,1:end-1);
ob = owner(:,2:end);
mask = oa > 0 & ob > 0 & oa ~= ob & ...
       isfinite(phzWrapped(:,1:end-1)) & isfinite(phzWrapped(:,2:end)) & ...
       isfinite(phzUnwrapped(:,1:end-1)) & isfinite(phzUnwrapped(:,2:end));

if any(mask(:))
    wa = phzWrapped(:,1:end-1);
    wb = phzWrapped(:,2:end);
    ua = phzUnwrapped(:,1:end-1);
    ub = phzUnwrapped(:,2:end);
    qa = cor(:,1:end-1);
    qb = cor(:,2:end);

    dphi = wrap_to_pi_local(wb(mask) - wa(mask));
    kest = round((ua(mask) + dphi - ub(mask)) ./ (2*pi));
    wgt  = 0.5 * (double(qa(mask)) + double(qb(mask)));

    a = double(oa(mask));
    b = double(ob(mask));
    [a, b, kest] = canonicalize_pairs(a, b, kest);

    pairA = [pairA; a];
    pairB = [pairB; b];
    pairK = [pairK; kest];
    pairW = [pairW; wgt];
end

edges = aggregate_pair_votes(pairA, pairB, pairK, pairW, minSupport);

end

% -------------------------------------------------------------------------
function screenModel = fit_basin_reference_screen(phzUnwrapped, cor, owner, pred, unwrapOpts)

mask = isfinite(phzUnwrapped) & isfinite(cor) & ...
       (cor >= unwrapOpts.basinScreenCohThresh) & owner > 0;

screenModel = struct();
screenModel.kind = 'none';
screenModel.phase = zeros(size(phzUnwrapped), 'like', phzUnwrapped);
screenModel.coeff = [];

if nnz(mask) < 20
    return;
end

% Prefer sensor-plane predictors when available
hasAz    = isfield(pred,'azCoord')    && ~isempty(pred.azCoord);
hasRange = isfield(pred,'rangeCoord') && ~isempty(pred.rangeCoord);

if hasAz && hasRange
    x1 = double(pred.azCoord(mask));
    x2 = double(pred.rangeCoord(mask));
    y  = double(phzUnwrapped(mask));

    X = [ones(size(x1)), x1, x2];
    w = double(cor(mask));
    w(~isfinite(w)) = 0;
    w = max(w, eps);

    beta = lscov(X, y, w);

    phase = nan(size(phzUnwrapped));
    x1f = double(pred.azCoord(:));
    x2f = double(pred.rangeCoord(:));
    validPred = isfinite(x1f) & isfinite(x2f);

    tmp = nan(numel(phzUnwrapped),1);
    tmp(validPred) = [ones(nnz(validPred),1), x1f(validPred), x2f(validPred)] * beta;
    phase(:) = tmp;

    screenModel.kind = 'sensorPlane';
    screenModel.phase = phase;
    screenModel.coeff = beta;
    return;
end

% Fallback to image-plane plane fit
[rr, cc] = ndgrid(1:size(phzUnwrapped,1), 1:size(phzUnwrapped,2));
x1 = double(rr(mask));
x2 = double(cc(mask));
y  = double(phzUnwrapped(mask));
X = [ones(size(x1)), x1, x2];

w = double(cor(mask));
w(~isfinite(w)) = 0;
w = max(w, eps);

beta = lscov(X, y, w);

phase = nan(size(phzUnwrapped));
phase(:) = [ones(numel(rr),1), double(rr(:)), double(cc(:))] * beta;

screenModel.kind = 'imagePlane';
screenModel.phase = phase;
screenModel.coeff = beta;

end

% -------------------------------------------------------------------------
function edges = build_gap_edges(owner, phzUnwrapped, cor, screenModel, unwrapOpts)

edges = struct('a',[],'b',[],'delta',[],'support',[],'weight',[]);

if ~isfield(screenModel,'phase') || isempty(screenModel.phase)
    return;
end

numBasins = double(max(owner(:)));
if numBasins <= 1
    return;
end

screen = screenModel.phase;
resid = phzUnwrapped - screen;

cohMask = isfinite(cor) & cor >= unwrapOpts.basinScreenCohThresh;
gapPix  = max(1, round(double(unwrapOpts.basinGapMaxPixels)));

% Precompute per-basin stats without storing full dilated masks
rowMin   = nan(numBasins,1);
rowMax   = nan(numBasins,1);
colMin   = nan(numBasins,1);
colMax   = nan(numBasins,1);
centRow  = nan(numBasins,1);
centCol  = nan(numBasins,1);
medResid = nan(numBasins,1);
madResid = nan(numBasins,1);
meanCoh  = nan(numBasins,1);
npix     = zeros(numBasins,1);

for b = 1:numBasins
    idx = find(owner == b);
    if isempty(idx)
        continue;
    end

    [rr, cc] = ind2sub(size(owner), idx);

    rowMin(b) = min(rr);
    rowMax(b) = max(rr);
    colMin(b) = min(cc);
    colMax(b) = max(cc);
    centRow(b) = mean(rr);
    centCol(b) = mean(cc);
    npix(b) = numel(idx);

    good = idx(cohMask(idx) & isfinite(resid(idx)));
    if ~isempty(good)
        vals = double(resid(good));
        medResid(b) = median(vals, 'omitnan');
        madResid(b) = median(abs(vals - medResid(b)), 'omitnan');
        meanCoh(b)  = mean(double(cor(good)), 'omitnan');
    end
end

pairA = [];
pairB = [];
pairK = [];
pairW = [];

for a = 1:numBasins-1
    if ~isfinite(medResid(a))
        continue;
    end

    % expanded bbox for basin a
    a_r0 = rowMin(a) - gapPix;
    a_r1 = rowMax(a) + gapPix;
    a_c0 = colMin(a) - gapPix;
    a_c1 = colMax(a) + gapPix;

    for b = a+1:numBasins
        if ~isfinite(medResid(b))
            continue;
        end

        % quick bbox intersection test
        overlapRows = ~(rowMax(b) < a_r0 || rowMin(b) > a_r1);
        overlapCols = ~(colMax(b) < a_c0 || colMin(b) > a_c1);

        if ~(overlapRows && overlapCols)
            continue;
        end

        % optional centroid sanity check
        dCent = hypot(centRow(a) - centRow(b), centCol(a) - centCol(b));
        if dCent > 4 * gapPix && min(npix([a b])) < 20
            continue;
        end

        dk = round((medResid(a) - medResid(b)) / (2*pi));

        spread = 0;
        if isfinite(madResid(a)), spread = spread + madResid(a); end
        if isfinite(madResid(b)), spread = spread + madResid(b); end
        spread = max(spread, 1e-3);

        support = min(npix(a), npix(b));
        support = max(1, min(support, 100));  % cap so huge basins don't dominate

        cohPair = mean([meanCoh(a), meanCoh(b)], 'omitnan');
        if ~isfinite(cohPair)
            cohPair = 0;
        end

        w = support * cohPair / spread;

        pairA = [pairA; a]; %#ok<AGROW>
        pairB = [pairB; b]; %#ok<AGROW>
        pairK = [pairK; dk]; %#ok<AGROW>
        pairW = [pairW; w]; %#ok<AGROW>
    end
end

edges = aggregate_pair_votes(pairA, pairB, pairK, pairW, 1);

end

% -------------------------------------------------------------------------
function anchorBasin = choose_anchor_basin(owner, cor, crResult, unwrapOpts)

numBasins = double(max(owner(:)));

switch lower(unwrapOpts.basinAnchorPolicy)
    case 'largest'
        counts = accumarray(double(owner(owner>0)), 1, [numBasins 1]);
        [~, anchorBasin] = max(counts);
        return;

    otherwise
        % try CR-first
end

anchorMask = false(size(owner));

% Generic search for usable CR support masks/locations
if isstruct(crResult) && ~isempty(crResult)
    if isfield(crResult,'referenceMask') && ~isempty(crResult.referenceMask)
        if isequal(size(crResult.referenceMask), size(owner))
            anchorMask = anchorMask | logical(crResult.referenceMask);
        end
    end

    if isfield(crResult,'validMask') && ~isempty(crResult.validMask)
        if isequal(size(crResult.validMask), size(owner))
            anchorMask = anchorMask | logical(crResult.validMask);
        end
    end

    if isfield(crResult,'mask') && ~isempty(crResult.mask)
        if isequal(size(crResult.mask), size(owner))
            anchorMask = anchorMask | logical(crResult.mask);
        end
    end

    % optional point-style CR locations
    if isfield(crResult,'row') && isfield(crResult,'col')
        rr = crResult.row(:);
        cc = crResult.col(:);
        good = isfinite(rr) & isfinite(cc);
        rr = round(rr(good));
        cc = round(cc(good));

        keep = rr >= 1 & rr <= size(owner,1) & cc >= 1 & cc <= size(owner,2);
        rr = rr(keep);
        cc = cc(keep);
        if ~isempty(rr)
            anchorMask(sub2ind(size(owner), rr, cc)) = true;
        end
    end
end

if any(anchorMask(:))
    score = zeros(numBasins,1);
    for b = 1:numBasins
        m = (owner == b) & anchorMask;
        if any(m(:))
            score(b) = nnz(m) + 1e-6 * sum(double(cor(m)), 'omitnan');
        end
    end
    if any(score > 0)
        [~, anchorBasin] = max(score);
        return;
    end
end

% fallback highest integrated coherence
score = zeros(numBasins,1);
for b = 1:numBasins
    m = owner == b;
    if any(m(:))
        score(b) = sum(double(cor(m)), 'omitnan');
    end
end
if any(score > 0)
    [~, anchorBasin] = max(score);
    return;
end

% final fallback: largest basin
counts = accumarray(double(owner(owner>0)), 1, [numBasins 1]);
[~, anchorBasin] = max(counts);

end

% -------------------------------------------------------------------------
function [shift, solved] = solve_basin_graph(numBasins, edges, anchorBasin)

shift = zeros(numBasins,1);
solved = false(numBasins,1);

if isempty(anchorBasin) || anchorBasin < 1 || anchorBasin > numBasins
    anchorBasin = 1;
end

solved(anchorBasin) = true;
shift(anchorBasin) = 0;

changed = true;
while changed
    changed = false;
    for e = 1:numel(edges.a)
        a = edges.a(e);
        b = edges.b(e);
        d = edges.delta(e);   % shift(b) - shift(a) = d

        if solved(a) && ~solved(b)
            shift(b) = shift(a) + d;
            solved(b) = true;
            changed = true;
        elseif solved(b) && ~solved(a)
            shift(a) = shift(b) - d;
            solved(a) = true;
            changed = true;
        end
    end
end

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function edges = concat_edges(edges1, edges2)

fields = {'a','b','delta','support','weight','pass'};
edges = struct();

for k = 1:numel(fields)
    f = fields{k};
    v1 = [];
    v2 = [];
    if isfield(edges1,f), v1 = edges1.(f); end
    if isfield(edges2,f), v2 = edges2.(f); end
    edges.(f) = [v1; v2];
end

end

% -------------------------------------------------------------------------
function edges = aggregate_pair_votes(pairA, pairB, pairK, pairW, minSupport)

edges = struct('a',[],'b',[],'delta',[],'support',[],'weight',[]);

if isempty(pairA)
    return;
end

pairs = [pairA(:), pairB(:)];
[uniquePairs, ~, ic] = unique(pairs, 'rows');

for p = 1:size(uniquePairs,1)
    sel = (ic == p);
    if nnz(sel) < minSupport
        continue;
    end

    kvals = pairK(sel);
    wvals = pairW(sel);

    d = weighted_mode_integer(kvals, wvals);

    edges.a(end+1,1) = uniquePairs(p,1); %#ok<AGROW>
    edges.b(end+1,1) = uniquePairs(p,2); %#ok<AGROW>
    edges.delta(end+1,1) = d; %#ok<AGROW>
    edges.support(end+1,1) = nnz(sel); %#ok<AGROW>
    edges.weight(end+1,1) = sum(wvals); %#ok<AGROW>
end

end

% -------------------------------------------------------------------------
function basinStats = compute_basin_stats(owner, phzUnwrapped, cor, unwrapOpts)

numBasins = double(max(owner(:)));

screenResidDummy = phzUnwrapped; %#ok<NASGU> % placeholder for consistent layout

rowMin   = nan(numBasins,1);
rowMax   = nan(numBasins,1);
colMin   = nan(numBasins,1);
colMax   = nan(numBasins,1);
centRow  = nan(numBasins,1);
centCol  = nan(numBasins,1);
npix     = zeros(numBasins,1);
meanCoh  = nan(numBasins,1);

linIdxByBasin = accumarray(double(owner(owner>0)), find(owner(owner>0)), [numBasins 1], @(v){v}, {[]});

for b = 1:numBasins
    idx = linIdxByBasin{b};
    if isempty(idx)
        continue;
    end

    [rr, cc] = ind2sub(size(owner), idx);

    rowMin(b) = min(rr);
    rowMax(b) = max(rr);
    colMin(b) = min(cc);
    colMax(b) = max(cc);
    centRow(b) = mean(rr);
    centCol(b) = mean(cc);
    npix(b) = numel(idx);
    meanCoh(b) = mean(double(cor(idx)), 'omitnan');
end

basinStats = struct();
basinStats.rowMin = rowMin;
basinStats.rowMax = rowMax;
basinStats.colMin = colMin;
basinStats.colMax = colMax;
basinStats.centRow = centRow;
basinStats.centCol = centCol;
basinStats.npix = npix;
basinStats.meanCoh = meanCoh;
basinStats.linIdx = linIdxByBasin;

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function edges = build_gap_edges_sparse( ...
    owner, phzUnwrapped, cor, screenModel, basinStats, solved, shift, unwrapOpts, passID)

edges = struct('a',[],'b',[],'delta',[],'support',[],'weight',[], ...
               'pass',[]);

if ~isfield(screenModel,'phase') || isempty(screenModel.phase)
    return;
end

numBasins = double(max(owner(:)));
if numBasins <= 1
    return;
end

screen = screenModel.phase;
resid = phzUnwrapped - screen;

cohMask = isfinite(cor) & cor >= unwrapOpts.basinScreenCohThresh;
gapPix  = max(1, round(double(unwrapOpts.basinGapMaxPixels)));
maxNbrs = unwrapOpts.basinGapMaxNeighbors;
minArea = unwrapOpts.basinMinGapArea;

medResid = nan(numBasins,1);
madResid = nan(numBasins,1);

for b = 1:numBasins
    idx = basinStats.linIdx{b};
    if isempty(idx)
        continue;
    end

    good = idx(cohMask(idx) & isfinite(resid(idx)));
    if isempty(good)
        continue;
    end

    vals = double(resid(good));
    medResid(b) = median(vals, 'omitnan');
    madResid(b) = median(abs(vals - medResid(b)), 'omitnan');
end

pairA = [];
pairB = [];
pairK = [];
pairW = [];
pairP = [];

unsolvedBasins = find(~solved);
solvedBasins   = find(solved);

for uu = 1:numel(unsolvedBasins)
    a = unsolvedBasins(uu);

    if basinStats.npix(a) < minArea
        continue;
    end
    if ~isfinite(medResid(a))
        continue;
    end

    a_r0 = basinStats.rowMin(a) - gapPix;
    a_r1 = basinStats.rowMax(a) + gapPix;
    a_c0 = basinStats.colMin(a) - gapPix;
    a_c1 = basinStats.colMax(a) + gapPix;

    candB = [];
    candScore = [];

    for ss = 1:numel(solvedBasins)
        b = solvedBasins(ss);

        if ~isfinite(medResid(b))
            continue;
        end

        overlapRows = ~(basinStats.rowMax(b) < a_r0 || basinStats.rowMin(b) > a_r1);
        overlapCols = ~(basinStats.colMax(b) < a_c0 || basinStats.colMin(b) > a_c1);

        if ~(overlapRows && overlapCols)
            continue;
        end

        dCent = hypot(basinStats.centRow(a) - basinStats.centRow(b), ...
                      basinStats.centCol(a) - basinStats.centCol(b));

        spread = 0;
        if isfinite(madResid(a)), spread = spread + madResid(a); end
        if isfinite(madResid(b)), spread = spread + madResid(b); end
        spread = max(spread, 1e-3);

        support = min(basinStats.npix(a), basinStats.npix(b));
        support = max(1, min(support, 100));

        cohPair = mean([basinStats.meanCoh(a), basinStats.meanCoh(b)], 'omitnan');
        if ~isfinite(cohPair), cohPair = 0; end

        % stronger score than before: prefer coherent, supported, nearby, low-spread basins
        score = cohPair * sqrt(double(support)) / ((1 + dCent) * (1 + spread));

        candB(end+1,1) = b; %#ok<AGROW>
        candScore(end+1,1) = score; %#ok<AGROW>
    end

    if isempty(candB)
        continue;
    end

    [~, ord] = sort(candScore, 'descend');
    ord = ord(1:min(numel(ord), maxNbrs));
    candB = candB(ord);

    for kk = 1:numel(candB)
        b = candB(kk);

        % shift(b) is already solved; infer shift(a)
        % medResid(a) + 2*pi*shift(a) ~= medResid(b) + 2*pi*shift(b)
        dk = round((medResid(a) - medResid(b)) / (2*pi));
        shiftA_from_b = shift(b) - dk;

        % convert to graph edge form: shift(b) - shift(a) = d
        d = shift(b) - shiftA_from_b;

        spread = 0;
        if isfinite(madResid(a)), spread = spread + madResid(a); end
        if isfinite(madResid(b)), spread = spread + madResid(b); end
        spread = max(spread, 1e-3);

        support = min(basinStats.npix(a), basinStats.npix(b));
        support = max(1, min(support, 100));

        cohPair = mean([basinStats.meanCoh(a), basinStats.meanCoh(b)], 'omitnan');
        if ~isfinite(cohPair), cohPair = 0; end

        w = support * cohPair / spread;

        aa = a;
        bb = b;
        dd = d;

        if aa > bb
            tmp = aa; aa = bb; bb = tmp;
            dd = -dd;
        end

        pairA(end+1,1) = aa; %#ok<AGROW>
        pairB(end+1,1) = bb; %#ok<AGROW>
        pairK(end+1,1) = dd; %#ok<AGROW>
        pairW(end+1,1) = w;  %#ok<AGROW>
        pairP(end+1,1) = passID; %#ok<AGROW>
    end
end

edges = aggregate_pair_votes(pairA, pairB, pairK, pairW, 1);

if ~isempty(edges.a)
    % preserve pass provenance for debugging
    edges.pass = repmat(passID, numel(edges.a), 1);
else
    edges.pass = [];
end

end

% -------------------------------------------------------------------------
function [a, b, k] = canonicalize_pairs(a, b, k)
flipMask = a > b;
tmp = a(flipMask);
a(flipMask) = b(flipMask);
b(flipMask) = tmp;
k(flipMask) = -k(flipMask);
end

% -------------------------------------------------------------------------
function m = weighted_mode_integer(kvals, wvals)
uk = unique(kvals(:));
score = zeros(size(uk));
for ii = 1:numel(uk)
    score(ii) = sum(wvals(kvals == uk(ii)));
end
[~, idx] = max(score);
m = uk(idx);
end

% -------------------------------------------------------------------------
function y = wrap_to_pi_local(x)
y = mod(x + pi, 2*pi) - pi;
end
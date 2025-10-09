function [SNR, NESg0, diagOut] = snr_from_gamma0_knn(gamma0, rslant, sAlong, lookMask, opts)
% Local low-percentile noise estimate using KNN / ball search in (r,s) space
% — robust for geocoded SLCs (no striping).
%
% Inputs
%   gamma0   [ny,nx]  : calibrated backscatter (linear units)
%   rslant   [ny,nx]  : slant range map (this SLC/pass)
%   sAlong   [ny,nx]  : along-track index / nearest-trajectory sample
%   lookMask [ny,nx]  : valid mask (logical)
%   opts (all optional)
%     .percentile = 2          % low percentile for NESγ0
%     .K          = 300        % target neighbors (KNN)
%     .minK       = 80         % minimum neighbors to accept
%     .rRadius    = 40         % max radius in range (m)   (used with rangesearch)
%     .sRadius    = 120        % max radius in along-track *index* units
%     .scale      = [20 60]    % feature scaling for [r, s] before KD-tree (m, index)
%     .batchN     = 2e5        % query batch size to limit memory
%     .winsorHi   = 99         % winsorize gamma0 within neighborhood (cap bright pts)
%
% Outputs
%   SNR     [ny,nx]  : gamma0 ./ NESγ0
%   NESg0   [ny,nx]  : estimated noise-equivalent gamma0 (linear)
%   diagOut struct   : diagnostics (neighbor counts, usedK, etc.)
%
% Requires Statistics and Machine Learning Toolbox (KDTreeSearcher/rangesearch)

if nargin<5, opts = struct; end
p  = getOpt(opts,'percentile',2);
K  = getOpt(opts,'K',300);
minK = getOpt(opts,'minK',80);
rRad = getOpt(opts,'rRadius',40);
sRad = getOpt(opts,'sRadius',120);
sc   = getOpt(opts,'scale',[20 60]);     % scale r by sc(1), s by sc(2)
batchN = getOpt(opts,'batchN',2e5);
winsHi = getOpt(opts,'winsorHi',99);

gamma0   = double(gamma0);
rslant   = double(rslant);
sAlong   = double(sAlong);
if nargin < 4 || isempty(lookMask), lookMask = true(size(gamma0)); end
valid = isfinite(gamma0) & isfinite(rslant) & isfinite(sAlong) & lookMask;

NESg0 = nan(size(gamma0));
nnUsed = zeros(size(gamma0));

if ~any(valid(:))
    SNR = NESg0;
    diagOut = struct('nnUsed',nnUsed,'opts',opts);
    return
end

% --- Build feature matrices (r,s) with scaling so KD-tree distance is balanced
r  = rslant(valid);   s = sAlong(valid);
F  = [r./sc(1), s./sc(2)];               % reference set (valid pixels)
Z  = gamma0(valid);

% KD-tree on reference set
Mdl = KDTreeSearcher(F,'Distance','euclidean');

% Query points = ALL valid pixels (same set), but we can also evaluate a subset if desired
[ny,nx] = size(gamma0);
[qI, qJ] = find(valid);
rq = rslant(valid);  sq = sAlong(valid);
Q  = [rq./sc(1), sq./sc(2)];

% Radii in scaled space
rad = hypot(rRad./sc(1), sRad./sc(2));

% Batch the queries to limit memory
Nq = size(Q,1);
b0 = 1;
while b0 <= Nq
    be = min(Nq, b0 + batchN - 1);
    Qb = Q(b0:be,:);
    idxCell = rangesearch(Mdl, Qb, rad);         % neighbors within radius
    % Ensure at least K neighbors by KNN fallback where needed
    needK = cellfun(@numel, idxCell) < minK;
    if any(needK)
        [idxK,~] = knnsearch(Mdl, Qb(needK,:), 'K', min(K, size(F,1)));
        % convert rows to cell
        addCell = mat2cell(idxK, ones(sum(needK),1), size(idxK,2));
        idxCell(needK) = addCell;
    end

    % Compute low-percentile per query (winsorize to avoid bright-point leakage)
    for t = 1:numel(idxCell)
        ids = idxCell{t};
        if isempty(ids), continue; end
        zz  = Z(ids);
        if ~isempty(winsHi) && isfinite(winsHi)
            cap = prctile(zz, winsHi);
            zz(zz>cap) = cap;
        end
        NESq = prctile(zz, p);
        NESg0(valid(b0 + t - 1)) = NESq;
        nnUsed(valid(b0 + t - 1)) = numel(ids);
    end

    b0 = be + 1;
end

% Final SNR
SNR = gamma0 ./ NESg0;

diagOut = struct('nnUsed',nnUsed,'opts',opts,'scaleUsed',sc,'radiusScaled',rad);
end

function v = getOpt(S,f,def), if isfield(S,f)&&~isempty(S.(f)), v=S.(f); else, v=def; end, end

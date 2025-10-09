function [SNR, NESg0, diagOut] = snr_from_gamma0_knn_fast(gamma0, rslant, sAlong, lookMask, opts)
% Local low-percentile NESγ0 via KNN / rangesearch in (r,s) space.
% Robust and faster: decimates queries and guarantees neighbors.
%
% Inputs
%   gamma0   [ny,nx]  linear backscatter
%   rslant   [ny,nx]  slant range
%   sAlong   [ny,nx]  along-track index (nearest traj sample, or meters)
%   lookMask [ny,nx]
%   opts (optional)
%     .percentile = 2.5
%     .K          = 300
%     .minK       = 80
%     .rRadius    = 40      % meters
%     .sRadius    = 120     % index units or meters (match sAlong units)
%     .scale      = [20 60] % scale r,s before KD-tree
%     .batchN     = 2e5
%     .winsorHi   = 99
%     .decimQuery = 2       % evaluate every Nth valid pixel (>=1), then NN-fill
%
% Outputs
%   SNR, NESg0, diagOut

if nargin<5, opts = struct; end
p   = getOpt(opts,'percentile',2.5);
K   = getOpt(opts,'K',300);
minK= getOpt(opts,'minK',80);
rRad= getOpt(opts,'rRadius',40);
sRad= getOpt(opts,'sRadius',120);
sc  = getOpt(opts,'scale',[20 60]);
batchN = getOpt(opts,'batchN',2e5);
winsHi = getOpt(opts,'winsorHi',99);
decQ = max(1, getOpt(opts,'decimQuery',50));  % 1 = no decimation

gamma0 = double(gamma0); rslant = double(rslant); sAlong = double(sAlong);
if nargin<4 || isempty(lookMask), lookMask = true(size(gamma0)); end
valid = isfinite(gamma0) & isfinite(rslant) & isfinite(sAlong) & lookMask;
NESg0  = nan(size(gamma0)); nnUsed = zeros(size(gamma0));

if ~any(valid(:))
    SNR = NESg0; diagOut=struct('reason','no valid pixels'); return
end

% Reference (r,s) & values
linIdx = find(valid);             % <-- correct writeback indices
r = rslant(linIdx); s = sAlong(linIdx); z = gamma0(linIdx);
F = [r./sc(1), s./sc(2)];         % KD-tree features
Mdl = KDTreeSearcher(F,'Distance','euclidean');

% Query subset (decimate to speed up)
Qidx = linIdx(1:decQ:end);
Qr   = rslant(Qidx); Qs = sAlong(Qidx);
Q    = [Qr./sc(1), Qs./sc(2)];
rad  = hypot(rRad./sc(1), sRad./sc(2));

% Compute NES only at decimated queries
NESsub = nan(numel(Qidx),1);
nnSub  = zeros(numel(Qidx),1);

b0 = 1; Nq = size(Q,1);
while b0 <= Nq
    be = min(Nq, b0+batchN-1);
    Qb = Q(b0:be,:);
    % radius search
    idxCell = rangesearch(Mdl, Qb, rad);

    % Ensure neighbors: fallback to KNN for empties
    empties = cellfun(@isempty, idxCell);
    if any(empties)
        [idxK,~] = knnsearch(Mdl, Qb(empties,:), 'K', min(K, size(F,1)));
        idxCell(empties) = mat2cell(idxK, ones(sum(empties),1), size(idxK,2));
    end

    % If still too small, second chance with 1.5x radius
    tooSmall = cellfun(@numel, idxCell) < minK;
    if any(tooSmall)
        idxCell2 = rangesearch(Mdl, Qb(tooSmall,:), 1.5*rad);
        % merge
        for t = 1:sum(tooSmall)
            if ~isempty(idxCell2{t})
                idxCell{find(tooSmall,1,'first')-1+t} = idxCell2{t};
            end
        end
        % final KNN fallback
        still = cellfun(@numel, idxCell) < minK;
        if any(still)
            [idxK2,~] = knnsearch(Mdl, Qb(still,:), 'K', min(K, size(F,1)));
            idxCell(still) = mat2cell(idxK2, ones(sum(still),1), size(idxK2,2));
        end
    end

    % Low-percentile with winsorization
    for t = 1:numel(idxCell)
        ids = idxCell{t};
        if isempty(ids), continue; end
        zz  = z(ids);
        if isfinite(winsHi)
            cap = prctile(zz, winsHi);
            zz(zz>cap) = cap;
        end
        NESsub(b0+t-1) = prctile(zz, p);
        nnSub(b0+t-1)  = numel(ids);
    end

    b0 = be + 1;
end
tiny = 1e-9;
NESsub(~isfinite(NESsub)) = NaN;        % mark bad now; filled later
NESsub(NESsub < tiny & isfinite(NESsub)) = tiny;

% Write back decimated NES, then nearest-neighbor fill to all valid pixels
NESg0(Qidx) = NESsub;

% if decQ > 1
%     ok = isfinite(NESsub);
%     if any(ok)
%         % KD-tree on decimated (r,s) for which we actually have NES values
%         Qr_ok = Qr(ok); Qs_ok = Qs(ok);
%         NES_ok = NESsub(ok);
%         nn_ok  = nnSub(ok);
% 
%         MdlQ = KDTreeSearcher([Qr_ok./sc(1), Qs_ok./sc(2)]);
%         [idxNear,~] = knnsearch(MdlQ, [r./sc(1), s./sc(2)]);   % size = numel(linIdx)
% 
%         % Map back to full valid set
%         NESg0(linIdx) = NES_ok(idxNear);
%         nnUsed(linIdx) = nn_ok(idxNear);
%     else
%         % no valid sub NES; bail gracefully
%         SNR = NESg0; 
%         diagOut = struct('reason','no valid NESsub after decimation');
%         return
%     end
% else
%     % no decimation: store counts for the computed subset
%     nnUsed(Qidx) = nnSub;
% end
if decQ > 1

ok = isfinite(NESsub);
if any(ok)
    Qr_ok = Qr(ok); Qs_ok = Qs(ok);
    NES_ok = NESsub(ok);
    nn_ok  = nnSub(ok);

    % 2) If ok is sparse, grow neighbors for bad queries (one-shot rescue)
    if nnz(ok) < numel(NESsub)
        bad = ~ok;
        if any(bad)
            % Try a larger KNN on the bad subset to compute a fallback NES
            [idxKbad,~] = knnsearch(Mdl, Q(bad,:), 'K', min(3*K, size(F,1)));
            for t = 1:sum(bad)
                ids = idxKbad(t,:);
                zz  = z(ids);
                cap = prctile(zz, winsHi);
                zz(zz>cap) = cap;
                NESsub(find(bad,1,'first')-1+t) = max(tiny, prctile(zz, p));
            end
            % refresh ok-set after rescue
            ok = isfinite(NESsub);
            Qr_ok = Qr(ok); Qs_ok = Qs(ok);
            NES_ok = NESsub(ok);
            nn_ok  = nnSub(ok);
        end
    end

    % 3) Nearest fill from the OK subset (guaranteed finite, tiny-floored)
    MdlQ = KDTreeSearcher([Qr_ok./sc(1), Qs_ok./sc(2)]);
    [idxNear,~] = knnsearch(MdlQ, [r./sc(1), s./sc(2)]);
    NESg0(linIdx)  = NES_ok(idxNear);
    nnUsed(linIdx) = nn_ok(idxNear);
else
    % as a last resort, set NES to small floor on all valid pixels
    NESg0(linIdx) = tiny;
    nnUsed(linIdx) = 0;
end
end
% Final clamp and SNR
%% === Robust NES cleanup + SNR computation (patch) ======================
validAll = isfinite(gamma0) & lookMask;
Gv = gamma0(validAll);

% 1) Scene-derived floor/ceiling (adaptive to this image)
nesFloorScene = max(prctile(Gv, 1), 1e-8);     % linear power
nesCeilScene  = max(prctile(Gv, 70), 1e-4);    % soft upper cap for NES

% Identify unresolved NES inside the swath after KNN fill
badNES = ~isfinite(NESg0) & lookMask;

% 2a) Local-window fallback (fills holes without far-away artifacts)
if any(badNES(:))
    win   = getOpt(opts,'localWin',25);         % odd int recommended
    half  = floor(win/2);
    need  = find(badNES);                        % linear indices

    [ny,nx] = size(gamma0);
    for k = 1:numel(need)
        [yy,xx] = ind2sub([ny,nx], need(k));
        r1 = max(1, yy-half); r2 = min(ny, yy+half);
        c1 = max(1, xx-half); c2 = min(nx, xx+half);

        patch = gamma0(r1:r2, c1:c2);
        mpatch= lookMask(r1:r2, c1:c2);
        patch = patch(isfinite(patch) & mpatch);

        if numel(patch) >= getOpt(opts,'localMin',50)
            NESg0(yy,xx) = max(prctile(patch, getOpt(opts,'percentile',2.5)), nesFloorScene);
        end
    end
end

% 2b) Optional slant-range fallback for any remaining holes
stillBad = ~isfinite(NESg0) & lookMask;
if any(stillBad(:)) && isfield(opts,'NES_slant') && ~isempty(opts.NES_slant)
    NESg0(stillBad) = opts.NES_slant(stillBad);
end

% 2c) Absolute last resort: scene floor
stillBad = ~isfinite(NESg0) & lookMask;
if any(stillBad(:))
    NESg0(stillBad) = nesFloorScene;
end

% 3) Clamp NES to sane range (prevents huge SNR)
NESg0(~lookMask) = NaN;
inside = isfinite(NESg0) & lookMask;
NESg0(inside) = min(max(NESg0(inside), nesFloorScene), nesCeilScene);

% 4) Compute SNR and clamp to a physical maximum (e.g., 30 dB ~ 1000 lin)
SNR = NaN(size(NESg0));
SNR(inside) = gamma0(inside) ./ NESg0(inside);

snrMaxLin = getOpt(opts,'snrMaxLin',1e3);       % ~30 dB
SNR(inside) = min(SNR(inside), snrMaxLin);
%% ======================================================================
% NESg0(~lookMask) = NaN;
% NESg0(isfinite(NESg0) & NESg0 < tiny) = tiny;
% SNR = gamma0 ./ NESg0;

% % Write back decimated NES, then nearest-neighbor fill to all valid pixels
% NESg0(Qidx) = NESsub;
% % Fill remaining valid pixels by nearest decimated query in (r,s)
% if decQ > 1
%     % KD-tree on decimated (r,s)
%     ok = isfinite(NESsub);
%     if any(ok)
%         MdlQ = KDTreeSearcher([Qr(ok)./sc(1), Qs(ok)./sc(2)]);
%         [idxNear,~] = knnsearch(MdlQ, [r./sc(1), s./sc(2)]);
%         NESg0(linIdx) = NESsub(ok(idxNear));   % nearest assignment
%         nnUsed(Qidx(ok)) = nnSub(ok);
%     else
%         % no valid sub NES; bail
%         SNR = NESg0; diagOut = struct('reason','no NESsub ok'); return
%     end
% else
%     % no decimation: store counts (optional)
%     nnUsed(Qidx) = nnSub;
% end

% % Mask out invalid
% NESg0(~lookMask) = NaN;
% SNR = gamma0 ./ NESg0;

diagOut = struct('nnUsed',nnUsed,'opts',opts,'radScaled',rad,'decimQuery',decQ);
end

function v = getOpt(S,f,def)
if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = def; end
end

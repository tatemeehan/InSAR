function [SNR, NESg0, prof] = snr_from_gamma0_by_slant2D(gamma0, rslant, sAlong, lookMask, opts)
% 2-D NESγ0 estimate in bins of (slant range r, along-track s), robust for geocoded SLCs.
% Inputs:
%   gamma0  : calibrated backscatter (linear), size [ny,nx]
%   rslant  : slant range map for THIS SLC, same size
%   sAlong  : along-track index/map (e.g., nearest-trajectory sample idx), same size
%   lookMask: logical valid mask
%   opts (all optional):
%     .dr        (default 5 m)       slant bin width
%     .ds        (default 30 idx)    along-track bin width (in sAlong units)
%     .prctLow   (default 2)         low percentile per bin
%     .minPerBin (default 400)       min samples per (r,s) bin
%     .smoothWin (default [5 5])     2-D moving median window in bins [Wr, Ws]
%     .clipNES   (default [1e-6 1])  clamp NESγ0
%
% Outputs:
%   SNR   : gamma0 ./ NESγ0
%   NESg0 : NESγ0 map
%   prof  : struct with axes and NES bin grid

if nargin<5, opts = struct; end
prc  = getOpt(opts,'prctLow',   2);
dr   = getOpt(opts,'dr',        11);
ds   = getOpt(opts,'ds',        11);
minN = getOpt(opts,'minPerBin', 5);
W    = getOpt(opts,'smoothWin', [5,5]);
CL   = getOpt(opts,'clipNES',   [1e-6, 1]);
% g = @(f,def) iff(isfield(opts,f)&&~isempty(opts.(f)), opts.(f), def);
% dr   = g('dr',5);
% ds   = g('ds',30);
% prc  = g('prctLow',2);
% minN = g('minPerBin',400);
% W    = g('smoothWin',[5 5]);
% CL   = g('clipNES',[1e-6 1]);

valid = isfinite(gamma0) & isfinite(rslant) & isfinite(sAlong) & lookMask;
SNR = nan(size(gamma0)); NESg0 = SNR; prof = struct();
if ~any(valid(:)), return; end

r = rslant(valid); s = sAlong(valid); z = gamma0(valid);
rmin = floor(min(r)); rmax = ceil(max(r));
smin = floor(min(s)); smax = ceil(max(s));
rEdges = rmin:dr:rmax;
sEdges = smin:ds:smax;
rCent  = rEdges(1:end-1) + dr/2;
sCent  = sEdges(1:end-1) + ds/2;

NES = nan(numel(rCent), numel(sCent));  % [Nr x Ns]
for ir = 1:numel(rCent)
    inR = r>=rEdges(ir) & r<rEdges(ir+1);
    if ~any(inR), continue; end
    sR  = s(inR); zR = z(inR);
    % bin along s
    [~,sb] = histc(sR, sEdges);
    for is = 1:numel(sCent)
        inBin = sb==is;
        if nnz(inBin) >= minN
            NES(ir,is) = prctile(zR(inBin), prc);
        end
    end
end

% fill gaps and smooth in 2-D bin space
% NESfill = fillmissing(fillmissing(NES,'nearest',2),'nearest',1);
NESfill = fillmissing2(NES,'natural');%,MissingLocations=isnan(NES));
Wr = max(3,W(1)); if mod(Wr,2)==0, Wr=Wr+1; end
Ws = max(3,W(2)); if mod(Ws,2)==0, Ws=Ws+1; end
NESsmooth = medfilt2(NESfill, [Wr Ws], 'symmetric');
NESsmooth = max(CL(1), min(CL(2), NESsmooth));

    % map back by 2-D interp
    NESg0 = interp2(sCent, rCent, NESsmooth, sAlong, rslant, 'makima', NaN);
NESg0(~lookMask) = NaN;
SNR = gamma0 ./ NESg0;

prof = struct('rCent',rCent,'sCent',sCent,'NES_raw',NES,'NES_smooth',NESsmooth);
end

function y = iff(c,a,b), if c, y=a; else, y=b; end
end

function val = getOpt(S, field, def)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        val = S.(field);
    else
        val = def;
    end
end


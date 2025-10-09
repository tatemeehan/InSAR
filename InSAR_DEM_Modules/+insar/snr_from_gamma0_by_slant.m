function [SNR, NESg0, prof] = snr_from_gamma0_by_slant(gamma0, rslant, lookMask, opts)
% Estimate SNR from calibrated gamma0 by deriving a noise-equivalent gamma0
% profile vs slant range (robust for geocoded images; no striping).
%
% Inputs
%   gamma0   : calibrated backscatter (linear units), [ny,nx]
%   rslant   : slant range map for THIS SLC (same size as gamma0)
%   lookMask : logical valid mask (same size). If empty, all-true.
%   opts     : struct with optional fields:
%       prctLow   (default 2)      low percentile per slant bin
%       dr        (default 5)      bin width in meters
%       minPerBin (default 500)    min pixels per bin
%       smoothWin (default 9)      smoothing window (bins, odd)
%       clipNES   (default [1e-6 1]) clamp NESγ0
%
% Outputs
%   SNR      : gamma0 ./ NESg0
%   NESg0    : NESγ0 map (same size as gamma0)
%   prof     : struct with profile diagnostics (centers, raw/smoothed NES)

    % ---- inputs & defaults
    if nargin < 3 || isempty(lookMask), lookMask = true(size(gamma0)); end
    if nargin < 4 || isempty(opts),     opts = struct();               end

    prc  = getOpt(opts,'prctLow',   2);
    dr   = getOpt(opts,'dr',        5);
    minN = getOpt(opts,'minPerBin', 500);
    W    = getOpt(opts,'smoothWin', 9);
    CL   = getOpt(opts,'clipNES',   [1e-6, 1]);

    gamma0  = double(gamma0);
    rslant  = double(rslant);
    lookMask = logical(lookMask);

    valid = isfinite(gamma0) & isfinite(rslant) & lookMask;
    if ~any(valid(:))
        SNR = nan(size(gamma0)); NESg0 = SNR; prof = struct(); return;
    end

    % ---- bin by slant range
    r = rslant(valid);
    z = gamma0(valid);

    rmin = floor(min(r)); rmax = ceil(max(r));
    edges = rmin:dr:rmax;
    if numel(edges) < 3
        % too few bins; fall back to global low-percentile
        nes = prctile(z, prc);
        NESg0 = nes * ones(size(gamma0), 'like', gamma0);
        SNR   = gamma0 ./ NESg0; NESg0(~lookMask) = NaN;
        prof = struct('centers',[], 'nes_raw',nes, 'nes_smooth',nes);
        return
    end

    cent  = edges(1:end-1) + dr/2;
    nes_raw = nan(size(cent));
    for k = 1:numel(cent)
        in = r>=edges(k) & r<edges(k+1);
        if nnz(in) >= minN
            nes_raw(k) = prctile(z(in), prc);
        end
    end

    % ---- fill gaps & smooth across bins
    ik = find(isfinite(nes_raw));
    if isempty(ik)
        SNR = nan(size(gamma0)); NESg0 = SNR;
        prof = struct('centers',cent,'nes_raw',nes_raw,'nes_smooth',nan(size(nes_raw)));
        return
    end

    nes_fill = interp1(cent(ik), nes_raw(ik), cent, 'linear','extrap');
    if mod(W,2)==0, W = W+1; end
    nes_smooth = smoothdata(nes_fill, 'movmedian', W);
    nes_smooth = max(CL(1), min(CL(2), nes_smooth));

    % ---- map back to pixels via slant range
    NESg0 = interp1(cent, nes_smooth, rslant, 'linear', 'extrap');
    NESg0(~lookMask) = NaN;

    SNR = gamma0 ./ NESg0;

    % ---- profile diagnostics
    prof = struct('centers',cent,'nes_raw',nes_raw,'nes_smooth',nes_smooth,'edges',edges);
end

function val = getOpt(S, field, def)
    if isstruct(S) && isfield(S, field) && ~isempty(S.(field))
        val = S.(field);
    else
        val = def;
    end
end

function [SNR, NESg0] = snr_from_gamma0_empirical(gamma0, lookMask, varargin)
% Estimate SNR map from calibrated gamma0 by deriving an NESγ0 profile
% from low-percentile statistics along range.
%
% gamma0   : calibrated backscatter (linear units)
% lookMask : valid pixels (logical)
% opts:
%   .percentile  (default 2)    -- low percentile per range to estimate noise floor
%   .minSamples  (default 200)  -- min valid pixels per column to trust percentile
%   .smoothWin   (default 101)  -- odd window for 1-D smoothing along range
%   .clipNESZ    (default [1e-6, 1]) -- clamp NESγ0 (in gamma0 units)
%
% Returns:
%   SNR   : gamma0 ./ NESγ0 (same size as input; NaN where invalid)
%   NESg0 : estimated NESγ0 map (same size as input)

opts = struct('percentile',2,'minSamples',200,'smoothWin',101,'clipNESZ',[1e-6,1]);
if ~isempty(varargin), opts = mergestruct(opts, varargin{:}); end

[ny,nx] = size(gamma0);
valid = isfinite(gamma0) & lookMask;

% Low-percentile per range column
nes_col = nan(1,nx);
for ix = 1:nx
    v = gamma0(valid(:,ix), ix);
    if numel(v) >= opts.minSamples
        nes_col(ix) = prctile(v, opts.percentile);
    end
end

% Smooth along range (replace gaps by nearest before smoothing)
idx = find(isfinite(nes_col));
if isempty(idx), NESg0 = nan(size(gamma0)); SNR = NESg0; return; end
nes_fill = interp1(idx, nes_col(idx), 1:nx, 'linear','extrap');

w = max(3, opts.smoothWin); if mod(w,2)==0, w=w+1; end
nes_smooth = smoothdata(nes_fill, 'movmean', w);

% Clamp to safe bounds
nes_smooth = max(opts.clipNESZ(1), min(opts.clipNESZ(2), nes_smooth));

% Broadcast to image and compute SNR
NESg0 = repmat(nes_smooth, ny, 1);
SNR   = nan(size(gamma0));
SNR(valid) = gamma0(valid) ./ NESg0(valid);

end

function s = mergestruct(s, s2)
f = fieldnames(s2);
for k=1:numel(f), s.(f{k}) = s2.(f{k}); end
end

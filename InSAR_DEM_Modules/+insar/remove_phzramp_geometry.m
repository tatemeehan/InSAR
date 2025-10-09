function [I_corr, model] = remove_phzramp_geometry(Ifilt, coh, geom, trajM, trajS, lambda, opts)
% Geometry ramp removal with FE->data scale calibration (wrap-safe).
% Uses your per-pixel closestIndex_master/slave (MxN), DEM (X,Y,Z), and trajectories.
%
% opts.sigmaFE   : Gaussian LP on φ_FE (default 60 px; use 30–100)
% opts.cohThresh : min coherence (default 0.30)
% opts.gammaMode : 'global' (default) or 'column' (per-column γ(j) smoothed)
% opts.gammaClip : [min max] clip for γ (default [0 2])
% opts.sigmaGamma: if column mode, 1D smooth along columns (default 80 px)

if nargin<7, opts = struct(); end
if ~isfield(opts,'sigmaFE'),    opts.sigmaFE    = 50;    end
if ~isfield(opts,'cohThresh'),  opts.cohThresh  = 0.30;  end
if ~isfield(opts,'gammaMode'),  opts.gammaMode  = 'global'; end
if ~isfield(opts,'gammaClip'),  opts.gammaClip  = [0 2]; end
if ~isfield(opts,'sigmaGamma'), opts.sigmaGamma = 80;    end

[M,N] = size(Ifilt); k = 4*pi/lambda;

X = geom.X; Y = geom.Y; Z = geom.Z;
idxM = double(geom.closestIndex_master);   % MxN (fractional/weighted ok)
idxS = double(geom.closestIndex_slave);

PM = trajM.pos;  PS = trajS.pos;  NsM = size(PM,1); NsS = size(PS,1);
tM = (1:NsM)';   tS = (1:NsS)';

% --- interpolate trajectories at fractional indices (assume 1-based fractional) ---
rmx = interp1(tM, PM(:,1), idxM, 'linear', NaN);
rmy = interp1(tM, PM(:,2), idxM, 'linear', NaN);
rmz = interp1(tM, PM(:,3), idxM, 'linear', NaN);
rsx = interp1(tS, PS(:,1), idxS, 'linear', NaN);
rsy = interp1(tS, PS(:,2), idxS, 'linear', NaN);
rsz = interp1(tS, PS(:,3), idxS, 'linear', NaN);

valid = isfinite(Ifilt) & isfinite(coh) & (coh >= opts.cohThresh) & ...
        isfinite(rmx)&isfinite(rmy)&isfinite(rmz)&isfinite(rsx)&isfinite(rsy)&isfinite(rsz);

% --- LOS and baseline per pixel ---
dx = X - rmx;  dy = Y - rmy;  dz = Z - rmz;

R  = sqrt(dx.^2 + dy.^2 + dz.^2);  R(R==0)=NaN;
sx = dx./R; sy = dy./R; sz = dz./R;

Bx = rsx - rmx;  By = rsy - rmy;  Bz = rsz - rmz;

% --- FE prediction (unwrapped), then gentle LP and gauge removal ---
phi_FE = -k * (Bx.*sx + By.*sy + Bz.*sz);
phi_FE(~valid) = NaN;
phi_FE = nanGaussLP(phi_FE, opts.sigmaFE);
phi_FE = phi_FE - median(phi_FE(valid),'omitnan');

% --- measured wrap-safe gradients from complex field ---
[ gx, gy, Wx, Wy ] = phaseGradAndWeights(Ifilt, coh);   % below

% --- predicted gradients from φ_FE ---
[ px, py ] = gradsFinite(phi_FE);                       % central/edge guarded

% --- WLS scale γ (global or per-column) ---
switch lower(opts.gammaMode)
    case 'global'
        maskX = isfinite(gx) & isfinite(px);  wX = Wx(maskX);
        maskY = isfinite(gy) & isfinite(py);  wY = Wy(maskY);
        num = nansum(wX .* px(maskX) .* gx(maskX)) + nansum(wY .* py(maskY) .* gy(maskY));
        den = nansum(wX .* px(maskX).^2)      + nansum(wY .* py(maskY).^2);
        gamma = num / max(den, eps).*5;
        gamma = min(max(gamma, opts.gammaClip(1)), opts.gammaClip(2));
        phi_pred = gamma * phi_FE;

    case 'column'
        gamma_j = nan(1,N);
        % horizontal edges (gx, px) live between j and j+1
        for j = 1:N
            numj = 0; denj = 0;
            if j <= N-1
                mx = isfinite(gx(:,j)) & isfinite(px(:,j));
                wx = Wx(mx,j);
                numj = numj + nansum(wx .* px(mx,j) .* gx(mx,j));
                denj = denj + nansum(wx .* px(mx,j).^2);
            end
            if j >= 2
                mx = isfinite(gx(:,j-1)) & isfinite(px(:,j-1));
                wx = Wx(mx,j-1);
                numj = numj + nansum(wx .* px(mx,j-1) .* gx(mx,j-1));
                denj = denj + nansum(wx .* px(mx,j-1).^2);
            end
            if denj > 0
                gamma_j(j) = numj / denj;
            end
        end
        % smooth γ along columns, clip, and expand to 2D
        gamma_j = smooth1D_nan(gamma_j(:), opts.sigmaGamma).';
        gamma_j = min(max(gamma_j, opts.gammaClip(1)), opts.gammaClip(2));
        phi_pred = phi_FE .* repmat(gamma_j, M, 1);

    otherwise
        error('Unknown gammaMode: %s', opts.gammaMode);
end

% --- subtract in complex domain ---
I_corr = Ifilt .* exp(-1i * phi_pred);

% pack
model.phi_FE   = phi_FE;
model.gamma    = strcmpi(opts.gammaMode,'global') * gamma + ...
                 strcmpi(opts.gammaMode,'column') * NaN;
if strcmpi(opts.gammaMode,'column'), model.gamma_j = gamma_j; end
model.phi_pred = phi_pred;
model.valid    = valid;
end

% ---- helpers ----
function [gx, gy, Wx, Wy] = phaseGradAndWeights(I, coh)
% I   : MxN complex interferogram
% coh : MxN coherence
% gx  : Mx(N-1) horizontal wrapped gradient (rad)
% gy  : (M-1)xN vertical wrapped gradient (rad)
% Wx,Wy : pair weights = sqrt(coh_pair1 * coh_pair2)

    [M,N] = size(I);
    U  = I ./ max(abs(I), eps);

    % Preallocate
    gx = NaN(M, N-1);  gy = NaN(M-1, N);
    Wx = NaN(M, N-1);  Wy = NaN(M-1, N);

    % ----- Horizontal pairs (j vs j+1) -----
    U1x = U(:,1:N-1);           U2x = U(:,2:N);
    C1x = coh(:,1:N-1);         C2x = coh(:,2:N);
    mx  = isfinite(U1x) & isfinite(U2x) & isfinite(C1x) & isfinite(C2x);

    tmpGx = U2x .* conj(U1x);           % complex ratio
    tmpWx = sqrt(max(C1x,0) .* max(C2x,0));

    gx(mx) = angle(tmpGx(mx));
    Wx(mx) = tmpWx(mx);

    % ----- Vertical pairs (i vs i+1) -----
    U1y = U(1:M-1,:);           U2y = U(2:M,:);
    C1y = coh(1:M-1,:);         C2y = coh(2:M,:);
    my  = isfinite(U1y) & isfinite(U2y) & isfinite(C1y) & isfinite(C2y);

    tmpGy = U2y .* conj(U1y);
    tmpWy = sqrt(max(C1y,0) .* max(C2y,0));

    gy(my) = angle(tmpGy(my));
    Wy(my) = tmpWy(my);
end


function [px, py] = gradsFinite(A)
% A   : MxN real field (e.g., phi_FE after LP)
% px  : Mx(N-1) horizontal diff A(:,j+1)-A(:,j)
% py  : (M-1)xN vertical   diff A(i+1,:)-A(i,:)

    [M,N] = size(A);

    % Horizontal
    A1x = A(:,1:N-1);  A2x = A(:,2:N);
    mx  = isfinite(A1x) & isfinite(A2x);
    Dx  = A2x - A1x;
    px  = NaN(M, N-1);
    px(mx) = Dx(mx);

    % Vertical
    A1y = A(1:M-1,:);  A2y = A(2:M,:);
    my  = isfinite(A1y) & isfinite(A2y);
    Dy  = A2y - A1y;
    py  = NaN(M-1, N);
    py(my) = Dy(my);
end


function Aout = nanGaussLP(A, sigma)
if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
L = max(3, ceil(5*sigma)); x = (-L:L);
g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
mask = isfinite(A); A(~mask)=0;
num = conv2(conv2(A, g, 'same'), g', 'same');
den = conv2(conv2(double(mask), g, 'same'), g', 'same');
Aout = num ./ max(den, eps); Aout(den==0)=NaN;
end

function b = smooth1D_nan(a, sigma)
if sigma<=0, b=a; return; end
L=max(3,ceil(5*sigma)); x=(-L:L); g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
num=conv(a2,g,'same'); den=conv(m2,g,'same');
b=num./max(den,eps); b(den==0)=NaN;
end

% function [I_corr, model] = remove_phzramp_geometry(Ifilt, coh, geom, trajM, trajS, lambda, opts)
% % Ifilt : MxN complex interferogram (filtered, pre-unwrapping)
% % coh   : MxN coherence (for masking only)
% % geom  : struct with fields X,Y,Z (MxN, meters, same frame as traj*),
% %         and closestIndex_master, closestIndex_slave (MxN, 1-based integers)
% % trajM.pos, trajS.pos : Ns x 3 positions (meters, lever arms applied) in SAME frame as X,Y,Z
% % lambda : wavelength (m), e.g., ~0.236 for L-band
% % opts.sigmaFE   : Gaussian LP on predicted phase (default 60 px, set 0 to skip)
% % opts.cohThresh : mask threshold (default 0.30)
% 
% if nargin<7, opts = struct(); end
% if ~isfield(opts,'sigmaFE'),   opts.sigmaFE   = 60;   end
% if ~isfield(opts,'cohThresh'), opts.cohThresh = 0.30; end
% 
% [M,N] = size(Ifilt);
% k = 4*pi/lambda;
% 
% % --- inputs ---
% X = geom.X; Y = geom.Y; Z = geom.Z;              % MxN
% idxM = double(geom.closestIndex_master);
% idxS = double(geom.closestIndex_slave);
% 
% % Trajectories (Ns x 3)
% PM = trajM.pos;   % [x y z], meters, NsM x 3
% PS = trajS.pos;   % [x y z], meters, NsS x 3
% NsM = size(PM,1); NsS = size(PS,1);
% 
% % Build two parameterizations:
% tM_idx = (1:NsM)';                   % sample-index parameter
% tS_idx = (1:NsS)';
% sM = [0; cumsum(vecnorm(diff(PM),2,2))];      % arc length (meters)
% sS = [0; cumsum(vecnorm(diff(PS),2,2))];
% 
% % If your indices are 0-based or use -1 for invalid, uncomment:
% % idxM(idxM<=0) = NaN; idxS(idxS<=0) = NaN; idxM = idxM+1; idxS = idxS+1;
% 
% NsM = size(trajM.pos,1);
% NsS = size(trajS.pos,1);
% 
% % --- validity mask: coherence + finite indices within bounds ---
% valid = isfinite(Ifilt) & isfinite(coh) & (coh >= opts.cohThresh);
% inbM  = isfinite(idxM) & idxM>=1 & idxM<=NsM;
% inbS  = isfinite(idxS) & idxS>=1 & idxS<=NsS;
% valid = valid & inbM & inbS;
% % 
% % % --- gather master/slave positions per pixel (vectorized, NaN-safe) ---
% % % pull coordinate vectors for convenient indexing
% % Mx = trajM.pos(:,1); My = trajM.pos(:,2); Mz = trajM.pos(:,3);
% % Sx = trajS.pos(:,1); Sy = trajS.pos(:,2); Sz = trajS.pos(:,3);
% % 
% % % Use placeholder index 1 where invalid to avoid indexing errors, then re-NaN
% % tmpM = idxM; tmpM(~valid) = 1;
% % tmpS = idxS; tmpS(~valid) = 1;
% % 
% % rmx = Mx(tmpM); rmy = My(tmpM); rmz = Mz(tmpM);   % MxN each
% % rsx = Sx(tmpS); rsy = Sy(tmpS); rsz = Sz(tmpS);
% % 
% % rmx(~valid)=NaN; rmy(~valid)=NaN; rmz(~valid)=NaN;
% % rsx(~valid)=NaN; rsy(~valid)=NaN; rsz(~valid)=NaN;
% 
% [tQM, modeM] = resolveParam(idxM, tM_idx, sM, NsM);
% [tQS, modeS] = resolveParam(idxS, tS_idx, sS, NsS);
% 
% % --- sample master/slave positions at per-pixel parameters (interp1) ---
% rmx = NaN(size(idxM)); rmy = rmx; rmz = rmx;
% rsx = NaN(size(idxS)); rsy = rsx; rsz = rsx;
% 
% maskM = isfinite(tQM);
% if strcmp(modeM,'index')
%     rmx(maskM) = interp1(tM_idx, PM(:,1), tQM(maskM), 'linear', NaN);
%     rmy(maskM) = interp1(tM_idx, PM(:,2), tQM(maskM), 'linear', NaN);
%     rmz(maskM) = interp1(tM_idx, PM(:,3), tQM(maskM), 'linear', NaN);
% else % 'meters'
%     rmx(maskM) = interp1(sM, PM(:,1), tQM(maskM), 'linear', NaN);
%     rmy(maskM) = interp1(sM, PM(:,2), tQM(maskM), 'linear', NaN);
%     rmz(maskM) = interp1(sM, PM(:,3), tQM(maskM), 'linear', NaN);
% end
% 
% maskS = isfinite(tQS);
% if strcmp(modeS,'index')
%     rsx(maskS) = interp1(tS_idx, PS(:,1), tQS(maskS), 'linear', NaN);
%     rsy(maskS) = interp1(tS_idx, PS(:,2), tQS(maskS), 'linear', NaN);
%     rsz(maskS) = interp1(tS_idx, PS(:,3), tQS(maskS), 'linear', NaN);
% else % 'meters'
%     rsx(maskS) = interp1(sS, PS(:,1), tQS(maskS), 'linear', NaN);
%     rsy(maskS) = interp1(sS, PS(:,2), tQS(maskS), 'linear', NaN);
%     rsz(maskS) = interp1(sS, PS(:,3), tQS(maskS), 'linear', NaN);
% end
% 
% % tighten validity: NaNs from interp are invalid
% valid = valid & isfinite(rmx) & isfinite(rmy) & isfinite(rmz) ...
%               & isfinite(rsx) & isfinite(rsy) & isfinite(rsz);
% 
% % --- LOS unit vector from master pos to DEM point (per pixel) ---
% dx = X - rmx;  dy = Y - rmy;  dz = Z - rmz;
% R  = sqrt(dx.^2 + dy.^2 + dz.^2);  R(R==0) = NaN;
% sx = dx ./ R; sy = dy ./ R; sz = dz ./ R;
% 
% % --- per-pixel baseline vector (slave - master) ---
% Bx = rsx - rmx;  By = rsy - rmy;  Bz = rsz - rmz;
% 
% % --- predicted flat-earth phase:  -k * (B · ŝ) ---
% phi_FE = -k * (Bx.*sx + By.*sy + Bz.*sz);
% phi_FE(~valid) = NaN;
% 
% % optional gentle low-pass on the PREDICTION only (keeps it long-wave)
% if opts.sigmaFE > 0
%     phi_FE = nanGaussLP(phi_FE, opts.sigmaFE);
% end
% 
% % remove global gauge for neatness
% phi_FE = phi_FE - median(phi_FE(valid), 'omitnan');
% 
% % --- subtract in complex domain (wrap-safe, preserves amplitude) ---
% I_corr = Ifilt .* exp(-1i * phi_FE);
% 
% % pack outputs
% model.phi_pred = phi_FE;
% model.valid    = valid;
% model.note     = 'Per-pixel geometry: phi = -(4π/λ) * (B_slave-master · ŝ_master→DEM)';
% end
% 
% % -------- helper: NaN-aware separable Gaussian LP --------
% function Aout = nanGaussLP(A, sigma)
% if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
% L = max(3, ceil(5*sigma)); x = (-L:L);
% g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
% mask = isfinite(A); A(~mask)=0;
% num = conv2(conv2(A, g, 'same'), g', 'same');
% den = conv2(conv2(double(mask), g, 'same'), g', 'same');
% Aout = num ./ max(den, eps); Aout(den==0)=NaN;
% end
% 
% % --- decide what idxM/idxS mean (heuristics) ---
% function [tQuery, mode] = resolveParam(idx, t_idx, s, Ns)
%     v = idx(isfinite(idx));
%     mode = 'index'; tQuery = idx;          % default
%     if isempty(v), tQuery(:) = NaN; return; end
% 
%     % zero-based integer?
%     if all(v>=0 & v<=Ns) && median(abs(v-round(v)),'omitnan') < 0.1 && any(v==0)
%         tQuery = idx + 1; mode='index';
%         tQuery(tQuery<1 | tQuery>Ns) = NaN; return;
%     end
%     % 1-based (possibly fractional)?
%     if all(v>0 & v<=Ns)
%         mode='index'; return;
%     end
%     % normalized 0..1 ?
%     if min(v) >= -0.05 && max(v) <= 1.05
%         tQuery = idx .* s(end);  mode='meters'; return;
%     end
%     % meters along track?
%     if min(v) >= -s(end)*0.05 && max(v) <= s(end)*1.05
%         tQuery = idx; mode='meters'; return;
%     end
%     % fallback: affine map to meters
%     a = (s(end)-s(1)) / (max(v)-min(v) + eps);
%     b = s(1) - a*min(v);
%     tQuery = a*idx + b; mode='meters';
% end
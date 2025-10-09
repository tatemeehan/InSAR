function [I_corr, model] = remove_phzramp_deltaB_column(Ifilt, coh, geom, trajM, lambda, opts)
% Physics-rooted residual ramp removal via small along-track baseline error δB(s).
% Ifilt : MxN complex interferogram (filtered)
% coh   : MxN coherence
% geom  : struct with X,Y,Z (MxN) and closestIndex_master (MxN)
% trajM.pos : Ns x 3 master positions (meters, lever arms applied), same frame as X,Y,Z
% lambda : wavelength (m)
% opts.cohThresh : default 0.30
% opts.sigmaLOS  : LP on ŝ to suppress noise (px), default 20
% opts.sigmaB    : 1D LP on δB(j) along columns, default 120
% opts.verbose   : true/false (print diagnostics)

if nargin<6, opts = struct(); end
if ~isfield(opts,'cohThresh'), opts.cohThresh = 0.30; end
if ~isfield(opts,'sigmaLOS'),  opts.sigmaLOS  = 20;   end
if ~isfield(opts,'sigmaB'),    opts.sigmaB    = 120;  end
if ~isfield(opts,'verbose'),   opts.verbose   = true; end

k = 4*pi/lambda;
[M,N] = size(Ifilt);

% --- 1) Build per-pixel LOS ŝ from master positions at fractional indices ---
X = geom.X; Y = geom.Y; Z = geom.Z;
idxM = double(geom.closestIndex_master);        % MxN, fractional OK
PM = trajM.pos;  NsM = size(PM,1);  tM = (1:NsM)';

% Interpolate master position at each pixel (chain-safe)
rmx = interp1(tM, PM(:,1), idxM, 'linear', NaN);
rmy = interp1(tM, PM(:,2), idxM, 'linear', NaN);
rmz = interp1(tM, PM(:,3), idxM, 'linear', NaN);

dx = X - rmx;  dy = Y - rmy;  dz = Z - rmz;
R  = sqrt(dx.^2 + dy.^2 + dz.^2);   R(R==0) = NaN;
sx = dx./R;  sy = dy./R;  sz = dz./R;

% Optional gentle LOS smoothing to kill pixel noise (NaN-robust)
sx = nanGaussLP(sx, opts.sigmaLOS);
sy = nanGaussLP(sy, opts.sigmaLOS);
sz = nanGaussLP(sz, opts.sigmaLOS);

valid = isfinite(Ifilt) & isfinite(coh) & (coh >= opts.cohThresh) & ...
        isfinite(sx) & isfinite(sy) & isfinite(sz);

% --- 2) Measured horizontal (between columns) wrap-safe phase gradient ---
[ gx, ~, Wx, ~ ] = phaseGradAndWeights(Ifilt, coh);  % Mx(N-1), Mx(N-1)

% --- 3) Predictors: Δŝ between adjacent columns (Mx(N-1)) ---
sx1 = sx(:,1:N-1); sx2 = sx(:,2:N);
sy1 = sy(:,1:N-1); sy2 = sy(:,2:N);
sz1 = sz(:,1:N-1); sz2 = sz(:,2:N);
ms  = isfinite(sx1)&isfinite(sx2)&isfinite(sy1)&isfinite(sy2)&isfinite(sz1)&isfinite(sz2);

dsx = NaN(M, N-1); dsy = dsx; dsz = dsx;
dsx(ms) = sx2(ms) - sx1(ms);
dsy(ms) = sy2(ms) - sy1(ms);
dsz(ms) = sz2(ms) - sz1(ms);

% --- 4) Solve per-column δB_edge(j) from gx ≈ -k*(Δŝ · δB_edge(j)) ---
dB_edge = nan(N-1,3);
for j = 1:N-1
    v = ms(:,j) & isfinite(gx(:,j)) & isfinite(Wx(:,j));
    if ~any(v), continue; end

    % Weighted robust WLS: y = X*c  with X = [-k*dsx, -k*dsy, -k*dsz]
    Xj = [ -k*dsx(v,j), -k*dsy(v,j), -k*dsz(v,j) ];
    yj = gx(v,j);
    wj = sqrt(max(Wx(v,j),0));     % pair weight
    dB_edge(j,:) = robustWLS3(Xj, yj, wj);
end

% --- 5) Put δB at column centers, smooth along columns ---
dB_col = nan(N,3);
if N>=2
    dB_col(1,:)     = dB_edge(1,:);
    dB_col(2:N-1,:) = 0.5*(dB_edge(1:N-2,:) + dB_edge(2:N-1,:));
    dB_col(N,:)     = dB_edge(end,:);
end
% Smooth δB along columns (NaN-aware)
for q=1:3, dB_col(:,q) = smooth1D_nan(dB_col(:,q), opts.sigmaB); end

% --- 6) Phase to subtract: φ_δB(x,y) = ŝ(x,y) · δB_col(j), then multiply by (-k) ---
phi_dB = nan(M,N);
for j = 1:N
    vj = isfinite(sx(:,j)) & isfinite(sy(:,j)) & isfinite(sz(:,j)) & all(isfinite(dB_col(j,:)));
    if any(vj)
        phi_dB(vj,j) = -k * ( sx(vj,j)*dB_col(j,1) + sy(vj,j)*dB_col(j,2) + sz(vj,j)*dB_col(j,3) );
    end
end
phi_dB = nanGaussLP(phi_dB, 10);                 % tiny LP for cleanliness
phi_dB = phi_dB - median(phi_dB(valid),'omitnan');  % remove gauge

% --- 7) Subtract in complex domain (wrap-safe) ---
I_corr = Ifilt .* exp(-1i * phi_dB);

% --- 8) Diagnostics (optional) ---
if opts.verbose
    [gx2,~,Wx2,~] = phaseGradAndWeights(I_corr, coh);
    m1 = isfinite(gx) & isfinite(Wx) & isfinite(dsx) & isfinite(dsy) & isfinite(dsz);
    m2 = isfinite(gx2)& isfinite(Wx2);
    rms1 = sqrt(nansum((sqrt(Wx(m1)).*gx(m1)).^2) / max(nnz(m1),1));
    rms2 = sqrt(nansum((sqrt(Wx2(m2)).*gx2(m2)).^2) / max(nnz(m2),1));
    fprintf('δB-column: RMS|gx| before %.4g, after %.4g rad\n', rms1, rms2);
end

% Output model
model.deltaB_col = dB_col;     % N x 3
model.phi_dB     = phi_dB;     % M x N (what was subtracted)
end

% ---------- helpers (chain-safe) ----------
function [gx, gy, Wx, Wy] = phaseGradAndWeights(I, coh)
    [M,N] = size(I);
    U = I ./ max(abs(I), eps);
    gx = NaN(M,N-1); gy = NaN(M-1,N);
    Wx = NaN(M,N-1); Wy = NaN(M-1,N);

    U1x = U(:,1:N-1);  U2x = U(:,2:N);
    C1x = coh(:,1:N-1); C2x = coh(:,2:N);
    mx  = isfinite(U1x)&isfinite(U2x)&isfinite(C1x)&isfinite(C2x);
    tmp = U2x .* conj(U1x);
    wpx = sqrt(max(C1x,0) .* max(C2x,0));
    gx(mx) = angle(tmp(mx));
    Wx(mx) = wpx(mx);

    U1y = U(1:M-1,:);  U2y = U(2:M,:);
    C1y = coh(1:M-1,:); C2y = coh(2:M,:);
    my  = isfinite(U1y)&isfinite(U2y)&isfinite(C1y)&isfinite(C2y);
    tmp = U2y .* conj(U1y);
    wpy = sqrt(max(C1y,0) .* max(C2y,0));
    gy(my) = angle(tmp(my));
    Wy(my) = wpy(my);
end

function c = robustWLS3(X, y, w)
    % Huber-IRLS on weighted 3-param LS
    w = w(:); y = y(:); X = X;
    Xw = bsxfun(@times, X, w); yw = w.*y;
    c = Xw \ yw; if any(~isfinite(c)), c = zeros(3,1); end
    delta = 1.345; maxit = 25;
    for it=1:maxit
        r = y - X*c; s = 1.4826*mad(r,1); if ~(isfinite(s)&&s>0), break; end
        u = r/s; psi = min(1, delta ./ max(abs(u), eps));
        ww = w .* psi;
        Xw = bsxfun(@times, X, ww); yw = ww.*y;
        cnew = Xw \ yw;
        if any(~isfinite(cnew)), break; end
        if norm(X*(cnew-c),2) <= 1e-6*max(1,norm(y,2)), c=cnew; break; end
        c = cnew;
    end
end

function Aout = nanGaussLP(A, sigma)
    if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
    L = max(3, ceil(5*sigma)); x = (-L:L);
    g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
    mask = isfinite(A); A(~mask)=0;
    num = conv2(conv2(A,g,'same'), g','same');
    den = conv2(conv2(double(mask),g,'same'), g','same');
    Aout = num ./ max(den, eps); Aout(den==0)=NaN;
end

function b = smooth1D_nan(a, sigma)
    if sigma<=0, b=a; return; end
    L = max(3, ceil(5*sigma)); x = (-L:L);
    g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
    m = isfinite(a); a2=a; a2(~m)=0; m2=double(m);
    num = conv(a2,g,'same'); den = conv(m2,g,'same');
    b = num ./ max(den, eps); b(den==0)=NaN;
end

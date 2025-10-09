function [Icorr, model] = remove_phzramp_geom2col(Ifilt, coh, dR, dTheta, opts)
% Column-wise geometry ramp: grad(phi) ~ betaR*grad(dR) + betaT*grad(dTheta).
% Ifilt, coh, dR, dTheta are MxN (rows=range, cols=azimuth).
% opts fields (all optional):
%   cohThresh (0.35), minSamples (150), sigmaAz (80), ridge (1e-3),
%   nearRangeSigma ([])  % in pixels; emphasize near range if set
%   addConst (false)     % allow a constant range gradient per column

if nargin<5, opts=struct(); end
if ~isfield(opts,'cohThresh'),     opts.cohThresh     = 0.35; end
if ~isfield(opts,'minSamples'),    opts.minSamples    = 150;  end
if ~isfield(opts,'sigmaAz'),       opts.sigmaAz       = 80;   end
if ~isfield(opts,'sigmaPhiAz'),  opts.sigmaPhiAz  = 40;   end % smooth φ_pred along azimuth
if ~isfield(opts,'ridge'),         opts.ridge         = 1e-3; end
if ~isfield(opts,'lambda'),      opts.lambda      = 1e-2; end % ridge strength
if ~isfield(opts,'nearRangeSigma'),opts.nearRangeSigma= [];   end
if ~isfield(opts,'addConst'),      opts.addConst      = false;end

[M,N] = size(Ifilt);

% ---- wrap-safe range gradient of phase (rad/px) ----
U   = Ifilt ./ max(abs(Ifilt), eps);
gF  = angle( U(2:end,:).*conj(U(1:end-1,:)) );      % (M-1) x N

% ---- predictor range gradients ----
gR   = diff(dR,     1, 1);                          % (M-1) x N
gTh  = diff(dTheta, 1, 1);                          % (M-1) x N

% ---- valid pairs and weights ----
Wpair = min(coh(2:end,:), coh(1:end-1,:));          % (M-1)xN
V     = isfinite(gF) & isfinite(gR) & isfinite(gTh) & (Wpair>=opts.cohThresh);
W     = Wpair; W(~V) = 0;

% Optional near-range emphasis
if ~isempty(opts.nearRangeSigma)
    rr = (1:(M-1))'; r0 = 1;  % “near range” at row 1
    wR = exp( - (rr - r0).^2 / (2*opts.nearRangeSigma^2) );
    W = W .* wR;
end

betaR = nan(1,N); betaT = nan(1,N); beta0 = nan(1,N);

for j = 1:N
    v = V(:,j);
    if nnz(v) < opts.minSamples, continue; end

% Robust standardize predictors (per column)
[mu1,s1] = robustCenterScale(gR(v,j));
[mu2,s2] = robustCenterScale(gTh(v,j));
X = [ (gR(v,j)-mu1)/s1 , (gTh(v,j)-mu2)/s2 ];
if opts.addConst, X = [ones(nnz(v),1) X]; end

    y  = gF(v,j);
    w  = W(v,j);
    [b, ok] = irls_ridge(X, y, w, opts.ridge);
    if ~ok, continue; end

    if opts.addConst
        beta0(j) = b(1);
        b = b(2:end);
    end
    % store back in original predictor units
    betaR(j) = b(1) / s1;
    betaT(j) = b(2) / s2;
end

% ---- smooth betas along azimuth ----
betaR_s = smoothNan1D(betaR, opts.sigmaAz);
betaT_s = smoothNan1D(betaT, opts.sigmaAz);
if opts.addConst
    beta0_s = smoothNan1D(beta0, opts.sigmaAz);
else
    beta0_s = zeros(1,N);
end

% ---- predicted gradient and integration to phase ----
gP = zeros(M-1, N);
for j = 1:N
    if ~isfinite(betaR_s(j)) || ~isfinite(betaT_s(j)), continue; end
    gP(:,j) = beta0_s(j) + betaR_s(j)*gR(:,j) + betaT_s(j)*gTh(:,j);
end

% Integrate along range (cumsum), then remove column mean (reference)
phi_pred = zeros(M,N);
% phi_pred(2:end,:) = cumsum(gP,1);
% gP: (M-1)xN predicted range gradients
% validMask: (M-1)xN true where gradient is valid/confident
validMask = isfinite(gP) & (min(coh(1:end-1,:),coh(2:end,:)) >= opts.cohThresh);

phi_pred = NaN(size(Ifilt));   % MxN
[M1,N] = size(gP);
M = M1 + 1;

for j = 1:N
    v = validMask(:,j);
    if ~any(v), continue; end

    % segment starts/ends within this column
    vp = [false; v; false];
    s  = find(diff(vp)==1);
    e  = find(diff(vp)==-1) - 1;

    col = NaN(M,1);
    for k = 1:numel(s)
        idx = s(k):e(k);            % indices in gP -> between rows r and r+1
        % phase at the segment’s first row is zero; integrate within the segment
        col(idx(1):(idx(end)+1)) = [0; cumsum(gP(idx,j))];
    end
    phi_pred(:,j) = col;
end

% Remove a column mean (only over confident pixels) so we don’t inject a bias
for j = 1:N
    vv = isfinite(phi_pred(:,j)) & (coh(:,j) >= opts.cohThresh);
    if any(vv), phi_pred(:,j) = phi_pred(:,j) - mean(phi_pred(vv,j)); end
end

for j=1:N
    vv = isfinite(phi_pred(:,j)) & coh(:,j)>=opts.cohThresh;
    if any(vv), phi_pred(:,j) = phi_pred(:,j) - mean(phi_pred(vv,j)); end
end
% --- NEW: gentle azimuth-only smoothing of predicted phase ---
if isfield(opts,'sigmaPhiAz') && opts.sigmaPhiAz > 0
    for i = 1:M
        % smooth across columns for each row; NaNs are handled
        phi_pred(i,:) = smoothNan1D(phi_pred(i,:), opts.sigmaPhiAz);
    end
end
% --- OPTIONAL: global scale k so geometry sets the shape only ---
if ~isfield(opts,'fitGlobalScale'), opts.fitGlobalScale = true; end
if opts.fitGlobalScale
    gA_raw  = angle( Ifilt(:,2:end) .* conj(Ifilt(:,1:end-1)) );           % az grad (rad/px)
    gA_pred = phi_pred(:,2:end) - phi_pred(:,1:end-1);                      % predicted az grad
    Waz     = min(coh(:,1:end-1), coh(:,2:end));
    m       = isfinite(gA_raw) & isfinite(gA_pred) & (Waz >= opts.cohThresh);

    % weighted least squares: minimize ||gA_raw - k*gA_pred||^2
    num = sum( Waz(m) .* gA_raw(m) .* gA_pred(m) );
    den = sum( Waz(m) .* (gA_pred(m).^2) );
    if den > 0
        k = num / den;
        phi_pred = k * phi_pred;
    end
end

% ---- complex subtraction ----
Icorr = Ifilt .* exp(-1i * phi_pred);

% ---- pack model + quick metric ----
model.betaR = betaR; model.betaT = betaT;
model.betaR_s = betaR_s; model.betaT_s = betaT_s;
model.beta0 = beta0; model.beta0_s = beta0_s;
model.phi_pred = phi_pred; model.opts = opts;

% quick report
gradA = angle( Ifilt(:,2:end) .* conj(Ifilt(:,1:end-1)) );   % M × (N-1)
gradR = angle( Ifilt(2:end,:) .* conj(Ifilt(1:end-1,:)) );   % (M-1) × N

maskA = isfinite(gradA) & min(coh(:,1:end-1),coh(:,2:end)) >= opts.cohThresh;
maskR = isfinite(gradR) & min(coh(1:end-1,:),coh(2:end,:)) >= opts.cohThresh;

medA_raw = median(abs(gradA(maskA)),'omitnan');
medR_raw = median(abs(gradR(maskR)),'omitnan');

% after correction (Icorr)
gradA_c = angle( Icorr(:,2:end) .* conj(Icorr(:,1:end-1)) );
gradR_c = angle( Icorr(2:end,:) .* conj(Icorr(1:end-1,:)) );
medA_cor = median(abs(gradA_c(maskA)),'omitnan');
medR_cor = median(abs(gradR_c(maskR)),'omitnan');
fprintf('med|∂φ/∂Az| raw=%.4f, corr=%.4f   med|∂φ/∂R| raw=%.4f, corr=%.4f rad/px\n', ...
        medA_raw, medA_cor, medR_raw, medR_cor);
end

% ---------- helpers ----------
function [b, ok] = irls_ridge(X, y, w, lam)
    ok = true;
    w  = w(:); y = y(:); 
    Xw = X .* w; yw = y .* w;
    H  = Xw.'*Xw; rhs = Xw.'*yw;
    lam = lam * max(trace(H)/size(H,1), eps);
    b  = (H + lam*eye(size(H))) \ rhs;
    if any(~isfinite(b)), ok=false; return; end
    % Huber IRLS
    delta = 1.345; maxit = 25;
    for it=1:maxit
        r = y - X*b;
        s = 1.4826*mad(r,1); if ~(isfinite(s) && s>0), break; end
        u = r/s; psi = min(1, delta./max(abs(u),eps));
        ww = w .* psi;
        Xw = X .* ww; yw = y .* ww;
        H  = Xw.'*Xw; rhs = Xw.'*yw;
        b2 = (H + lam*eye(size(H))) \ rhs;
        if any(~isfinite(b2)), break; end
        if norm(X*(b2-b),2) <= 1e-6*max(1,norm(y,2)), b = b2; break; end
        b = b2;
    end
end

function s = smoothNan1D(a, sigma)
    if sigma<=0, s=a; return; end
    L=max(3,ceil(5*sigma)); x=(-L:L);
    g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
    num=conv(a2,g,'same'); den=conv(m2,g,'same');
    s=num./max(den,eps); s(den==0)=NaN;
end

function [mu,s] = robustCenterScale(x)
    mu = median(x);
    s  = 1.4826*mad(x,1);
    if ~isfinite(s) || s < eps
        s = std(x);
    end
    if ~isfinite(s) || s < eps
        s = 1;
    end
end
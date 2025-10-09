function [Icorr, model] = remove_phzramp_geomScaleCol(Ifilt, coh, dR, opts)
% Column-wise geometry scaling using range gradients.
% Ifilt, coh, dR: MxN.  opts: cohThresh, minSamples, sigmaAz, ridge.

if nargin<4, opts=struct(); end
if ~isfield(opts,'cohThresh'),  opts.cohThresh  = 0.35; end
if ~isfield(opts,'minSamples'), opts.minSamples = 150;  end
if ~isfield(opts,'sigmaAz'),    opts.sigmaAz    = 80;   end
if ~isfield(opts,'ridge'),      opts.ridge      = 1e-3; end

[M,N] = size(Ifilt);
% wrap-safe range gradient of phase (rad/px)
U  = Ifilt ./ max(abs(Ifilt),eps);
gF = angle( U(2:end,:).*conj(U(1:end-1,:)) );        % (M-1)xN
% range gradient of ΔR (m/px) on same pairs of rows
gdR = diff(dR,1,1);                                   % (M-1)xN

% valid mask for gradient pairs
V = isfinite(gF) & isfinite(gdR);
W = min(coh(2:end,:),coh(1:end-1,:));                 % pair weight
V = V & (W>=opts.cohThresh);
w = W; w(~V)=0;

gamma = nan(1,N);     % per-column scale
for j=1:N
    v = V(:,j);
    if nnz(v)<opts.minSamples, continue; end
    x = gdR(v,j);                     % predictor (m/px)
    y = gF(v,j);                      % target (rad/px)
    ww = w(v,j);
    % weighted ridge scalar regression: gamma = (x'Wx)/(x'Wx + λ)* (x'Wy)/(x'Wx)
    Xw = x.*ww; Yw = y.*ww;
    xx = sum(Xw.*x);  xy = sum(Xw.*y);
    lam = opts.ridge * max(xx,eps);
    gamma(j) = xy / (xx + lam);       % robust-ish; small ridge
end

% smooth gamma along azimuth to follow slow B⊥ drift
gamma_s = smoothNan1D(gamma, opts.sigmaAz);

% predicted ramp and complex correction
phi_pred = zeros(M,N);
for j=1:N
    if ~isfinite(gamma_s(j)), continue; end
    pr = gamma_s(j) * dR(:,j);        % rad = (rad/m)*m  (units via gamma)
    % remove column mean so we don't change overall reference
    pr = pr - mean(pr(isfinite(pr) & coh(:,j)>=opts.cohThresh), 'omitnan');
    phi_pred(:,j) = pr;
end
Icorr = Ifilt .* exp(-1i*phi_pred);

model.gamma = gamma; model.gamma_s = gamma_s; model.phi_pred = phi_pred;
end

function s = smoothNan1D(a, sigma)
    if sigma<=0, s=a; return; end
    L=max(3,ceil(5*sigma)); x=(-L:L);
    g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
    num=conv(a2,g,'same'); den=conv(m2,g,'same');
    s=num./max(den,eps); s(den==0)=NaN;
end

function [I_corr, model] = remove_phzramp_rangeRidge(Ifilt, coh, R, opts)
% Coherence-weighted per-column φ ≈ a_j + b_j*R(x,j), ridge-stable.
% Ifilt, coh, R are MxN; R = slant range per pixel.
% opts: cohThresh(0.35), minSamples(100), sigmaAB(150), ridge(1e-3), useFE(false), phiFE([])

if nargin<4, opts = struct(); end
if ~isfield(opts,'cohThresh'),  opts.cohThresh = 0.35; end
if ~isfield(opts,'minSamples'), opts.minSamples = 100; end
if ~isfield(opts,'sigmaAB'),    opts.sigmaAB   = 150;  end
if ~isfield(opts,'ridge'),      opts.ridge     = 1e-3; end
if ~isfield(opts,'useFE'),      opts.useFE     = false; end
if ~isfield(opts,'phiFE'),      opts.phiFE     = []; end

[M,N] = size(Ifilt);
phW = angle(Ifilt);

% --- unwrap per column (only where coherent) ---
ph = NaN(M,N);
mask = isfinite(coh) & coh >= opts.cohThresh & isfinite(phW);
for j=1:N
    v = mask(:,j);
    if nnz(v) >= opts.minSamples
        ph(v,j) = unwrap(phW(v,j));
    end
end

% --- LS per column with ridge & coherence weights ---
a = nan(1,N); b = nan(1,N);
for j=1:N
    v = isfinite(ph(:,j)) & isfinite(R(:,j));
    if nnz(v) < opts.minSamples, continue; end

    y = ph(v,j);
    x = R(v,j);

    % (optional) remove column medians to make intercept stable
    y = y - median(y,'omitnan');
    x = x - median(x,'omitnan');

    % design with intercept
    X = [ones(nnz(v),1), x];

    % coherence weights
    w = coh(v,j);
    w = w / max(median(w,'omitnan'), eps);
    Xw = X .* w; yw = y .* w;

    % ridge (scaled by X'X trace to be unitless)
    H = Xw.'*Xw;
    lam = opts.ridge * trace(H)/size(H,1);
    beta = (H + lam*eye(2)) \ (Xw.'*yw);
    if all(isfinite(beta))
        a(j) = beta(1);
        b(j) = beta(2);
    end
end

% --- smooth along columns (NaN-aware) ---
a = smooth1D_nan(a, opts.sigmaAB);
b = smooth1D_nan(b, opts.sigmaAB);

% --- predict & subtract; optionally add FE term with tiny scale γ(j) ---
phi_pred = repmat(a, M, 1) + R .* repmat(b, M, 1);

if opts.useFE && ~isempty(opts.phiFE)
    % fit a tiny per-column γ on residual after removing a+ bR
    phiFE = opts.phiFE; phiFE(~isfinite(phiFE)) = NaN;
    gam = nan(1,N);
    for j=1:N
        v = isfinite(ph(:,j)) & isfinite(phiFE(:,j));
        if nnz(v) < opts.minSamples, continue; end
        y = (ph(v,j) - (a(j) + b(j)*R(v,j))) - median(ph(v,j) - (a(j) + b(j)*R(v,j)),'omitnan');
        x = phiFE(v,j) - median(phiFE(v,j),'omitnan');
        w = coh(v,j); w = w / max(median(w,'omitnan'), eps);
        Xw = x .* w;  yw = y .* w;
        den = Xw.'*x; num = Xw.'*y;
        if den>0, gam(j) = num/den; end
    end
    gam = smooth1D_nan(gam, 200);
    phi_pred = phi_pred + phiFE .* repmat(gam, M, 1);
end

% subtract in complex domain (wrap-safe)
I_corr = Ifilt .* exp(-1i*phi_pred);

% pack
model.a = a; model.b = b; model.phi_pred = phi_pred;
if exist('gam','var'), model.gamma = gam; end
end

% --- helpers ---
function s = smooth1D_nan(a, sigma)
    if sigma<=0, s=a; return; end
    L=max(3,ceil(5*sigma)); x=(-L:L);
    g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
    num=conv(a2,g,'same'); den=conv(m2,g,'same');
    s = num ./ max(den,eps); s(den==0)=NaN;
end

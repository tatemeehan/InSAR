function [Icorr, model] = remove_phzramp_sinusoid(Ifilt, coh, lambda_px, opts)
% Coherence-weighted per-column sinusoid removal at f0 = 1/lambda_px.
% Ifilt, coh: MxN; lambda_px in pixels (≈ 9.5 from your spectrum).
% opts: cohThresh(0.35), minSamples(120), sigmaAz(60), ridge(1e-3),
%       includeTrend(false) — fit DC/slope but don't subtract them.

if nargin<4, opts=struct(); end
if ~isfield(opts,'cohThresh'),    opts.cohThresh    = 0.35; end
if ~isfield(opts,'minSamples'),   opts.minSamples   = 120;  end
if ~isfield(opts,'sigmaAz'),      opts.sigmaAz      = 60;   end
if ~isfield(opts,'ridge'),        opts.ridge        = 1e-3; end
if ~isfield(opts,'includeTrend'), opts.includeTrend = false; end

[M,N] = size(Ifilt);
f0 = 1/max(lambda_px, eps);                 % cycles/pixel
y  = (1:M)';                                % range index

phW = angle(Ifilt);
mask = isfinite(coh) & coh>=opts.cohThresh & isfinite(phW);

A = nan(1,N); B = nan(1,N);     % cosine/sine amplitudes
C = nan(1,N); D = nan(1,N);     % optional DC/slope (fitted, not removed)

for j=1:N
    v = mask(:,j);
    if nnz(v) < opts.minSamples, continue; end

    % unwrap only valid rows; center indices to stabilize intercepts
    yy = y(v);  y0 = yy - median(yy);
    phi = unwrap(phW(v,j));
    phi = phi - median(phi,'omitnan');

    % design matrix: [cos, sin, (optional 1), (optional y0)]
    X = [cos(2*pi*f0*y0), sin(2*pi*f0*y0)];
    if opts.includeTrend, X = [X, ones(nnz(v),1), y0]; end

    % coherence weights + ridge for numerical stability
    w  = coh(v,j); w = w / max(median(w,'omitnan'), eps);
    Xw = X .* w; yw = phi .* w;
    H = Xw.'*Xw;
    lam = opts.ridge * trace(H)/size(H,1);
    beta = (H + lam*eye(size(H))) \ (Xw.'*yw);
    if ~all(isfinite(beta)), continue; end

    A(j) = beta(1);  B(j) = beta(2);
    if opts.includeTrend, C(j)=beta(3); D(j)=beta(4); end
end

% smooth amplitudes along azimuth to suppress striping
A = smooth1D_nan(A, opts.sigmaAz);
B = smooth1D_nan(B, opts.sigmaAz);

% build predicted ripple and subtract in complex domain
Icorr   = Ifilt;
phi_rip = zeros(M,N,'like',phW);
for j=1:N
    v = isfinite(phW(:,j));
    yy = y(v);  y0 = yy - median(yy);
    if ~isfinite(A(j)) || ~isfinite(B(j)), continue; end
    pr = A(j)*cos(2*pi*f0*y0) + B(j)*sin(2*pi*f0*y0);
    phi_rip(v,j) = pr;
    Icorr(v,j)   = Ifilt(v,j) .* exp(-1i*pr);
end

model.f0 = f0; model.lambda_px = lambda_px;
model.A = A; model.B = B; model.C = C; model.D = D;
model.phi_rip = phi_rip;
end

function s = smooth1D_nan(a, sigma)
    if sigma<=0, s=a; return; end
    L=max(3,ceil(5*sigma)); x=(-L:L);
    g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
    num=conv(a2,g,'same'); den=conv(m2,g,'same');
    s=num./max(den,eps); s(den==0)=NaN;
end

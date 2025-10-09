function [Icorr, model] = remove_phzramp_azRamp(Ifilt, coh, opts)
% Wrap-safe per-column ramp (in range) whose slope varies with azimuth.
% Ifilt, coh: MxN (rows=range, cols=azimuth)
% opts: deg (1 or 2), cohThresh(0.35), minSamples(120), sigmaAz(60), ridge(1e-3)

if nargin<3, opts=struct(); end
if ~isfield(opts,'deg'),        opts.deg        = 1;   end  % try 1, then 2
if ~isfield(opts,'cohThresh'),  opts.cohThresh  = 0.35;end
if ~isfield(opts,'minSamples'), opts.minSamples = 120; end
if ~isfield(opts,'sigmaAz'),    opts.sigmaAz    = 80;  end % smooth along az
if ~isfield(opts,'ridge'),      opts.ridge      = 1e-3;end

[M,N] = size(Ifilt);
phW   = angle(Ifilt);
mask  = isfinite(coh) & coh>=opts.cohThresh & isfinite(phW);
r     = (1:M)';                             % "range" index (px or meters)
Icorr = Ifilt;

% Coefficients per column (azimuth)
b0 = nan(1,N);                 % intercept (we do NOT subtract it)
b1 = nan(1,N);                 % slope vs range
b2 = nan(1,N);                 % optional quadratic

for j=1:N
    v = mask(:,j);
    if nnz(v) < opts.minSamples, continue; end

    rr  = r(v);   r0 = median(rr);   x = rr - r0;
    phi = unwrap(phW(v,j));
    phi = phi - median(phi,'omitnan');
    w   = coh(v,j); w = w / max(median(w,'omitnan'), eps);

    % Design matrix (centered in range)
    if opts.deg==1
        X = [ones(nnz(v),1), x];
    else
        X = [ones(nnz(v),1), x, x.^2];
    end

    % Weighted ridge regression (robust to ill-conditioning)
    Xw = X .* w;  yw = phi .* w;
    H  = Xw.'*Xw;
    lam = opts.ridge * trace(H)/size(H,1);
    beta = (H + lam*eye(size(H))) \ (Xw.'*yw);
    if ~all(isfinite(beta)), continue; end

    b0(j)=beta(1);  b1(j)=beta(2);
    if opts.deg>1, b2(j)=beta(3); end
end

% Smooth coefficients along azimuth to avoid striping
b1s = smoothNan1D(b1, opts.sigmaAz);
if opts.deg>1, b2s = smoothNan1D(b2, opts.sigmaAz); else, b2s = []; end

% Build predicted ramp and subtract (complex-domain, wrap-safe)
phi_pred = zeros(M,N,'like',phW);
for j=1:N
    if ~isfinite(b1s(j)), continue; end
    rr = r; r0 = median(r);
    x  = rr - r0;
    pr = b1s(j)*x;
    if opts.deg>1 && isfinite(b2s(j)), pr = pr + b2s(j)*x.^2; end
    phi_pred(:,j) = pr;                 % no intercept subtraction
    Icorr(:,j)    = Ifilt(:,j) .* exp(-1i*pr);
end

model.b0=b0; model.b1=b1; model.b2=b2;
model.b1s=b1s; model.b2s=b2s; model.phi_pred=phi_pred;
model.opts=opts;
end

function s = smoothNan1D(a, sigma)
    if sigma<=0, s=a; return; end
    L=max(3,ceil(5*sigma)); x=(-L:L);
    g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    m=isfinite(a); a2=a; a2(~m)=0; m2=double(m);
    num=conv(a2,g,'same'); den=conv(m2,g,'same');
    s=num./max(den,eps); s(den==0)=NaN;
end

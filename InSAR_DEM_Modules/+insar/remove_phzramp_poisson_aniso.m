function [I_corr, model] = remove_phzramp_poisson_aniso(I_filt, coh, ...
        sigmaX, sigmaY, lambda, rslant, beta, p, cohThresh)

% I_filt  (MxN) complex interferogram
% coh     (MxN) coherence
% sigmaX  smoothing (px) across RANGE (columns)
% sigmaY  smoothing (px) across AZIMUTH (rows)  [make this >> sigmaX]
% lambda  ridge for the normal equations (e.g., 1e-3)
% rslant  (MxN) slant range (optional; if empty, uses column index)
% beta,p  taper strength & shape toward near range (default beta=2, p=2)
% cohThresh mask threshold (default 0.30)

if nargin<3 || isempty(sigmaX),   sigmaX = 50;  end
if nargin<4 || isempty(sigmaY),   sigmaY = 150; end
if nargin<5 || isempty(lambda),   lambda = 1e-3; end
if nargin<6, rslant = []; end
if nargin<7 || isempty(beta),     beta   = 2;   end
if nargin<8 || isempty(p),        p      = 2;   end
if nargin<9 || isempty(cohThresh),cohThresh = 0.10; end

[M,N] = size(I_filt);
valid = isfinite(I_filt) & isfinite(coh) & (coh>=cohThresh);

% ---------------- unit phasor & edge masks ----------------
U  = I_filt ./ max(abs(I_filt), eps);
ex = valid(:,1:N-1) & valid(:,2:N);   % horizontal (range) edges
ey = valid(1:M-1,:) & valid(2:M,:);   % vertical  (azimuth) edges

% complex edge phasors (no chained indexing)
Gx = NaN(M,N-1); Gy = NaN(M-1,N);
tmpx = U(:,2:N).*conj(U(:,1:N-1));  Gx(ex) = tmpx(ex);
tmpy = U(2:M,:).*conj(U(1:M-1,:));  Gy(ey) = tmpy(ey);

% coherence edge weights
Wx = NaN(M,N-1); Wy = NaN(M-1,N);
tmpwx = sqrt(coh(:,2:N).*coh(:,1:N-1));  Wx(ex)=tmpwx(ex);
tmpwy = sqrt(coh(2:M,:).*coh(1:M-1,:));  Wy(ey)=tmpwy(ey);

% ---------------- near-range taper (τ >= 1) ----------------
if ~isempty(rslant) && any(isfinite(rslant(:)))
    r = rslant; r(~valid) = NaN;
    rmin = prctile(r(valid), 2); rmax = prctile(r(valid), 98);
    s = (rmax - r) ./ max(rmax - rmin, eps); % 0 far → 1 near
    % horizontal edges live between columns c and c+1: average s
    Sx = 0.5*(s(:,1:N-1) + s(:,2:N));
    % vertical edges between rows i and i+1: average s
    Sy = 0.5*(s(1:M-1,:) + s(2:M,:));
else
    % fallback: use column index as a proxy for range (col 1 = near)
    [~,col] = ndgrid(1:M, 1:N);
    s_col = (max(col(valid)) - col) ./ max(max(col(valid))-min(col(valid)),1);
    Sx = 0.5*(s_col(:,1:N-1) + s_col(:,2:N));
    Sy = 0.5*(s_col(1:M-1,:) + s_col(2:M,:));
end
tauX = 1 + beta * (max(Sx,0)).^p;   % strengthen near range
tauY = 1 + beta * (max(Sy,0)).^p;

Wx = Wx .* tauX;  Wy = Wy .* tauY;

% ---------------- anisotropic complex smoothing ----------------
Gx_s = circSmoothAniso(Gx, Wx, sigmaY, sigmaX);  % rows=sigmaY, cols=sigmaX
Gy_s = circSmoothAniso(Gy, Wy, sigmaY, sigmaX);

% target gradients (radians)
vX = isfinite(Gx_s)&isfinite(Wx)&(Wx>0);
vY = isfinite(Gy_s)&isfinite(Wy)&(Wy>0);
gx_t = NaN(M,N-1); gy_t = NaN(M-1,N);
gx_t(vX) = angle(Gx_s(vX));
gy_t(vY) = angle(Gy_s(vY));

% ---------------- mean-gradient plane anchor ----------------
ax = 0; by = 0;
if any(vX(:)), ax = sum(Wx(vX).*gx_t(vX)) / max(sum(Wx(vX)), eps); end
if any(vY(:)), by = sum(Wy(vY).*gy_t(vY)) / max(sum(Wy(vY)), eps); end
[xg,yg] = meshgrid(1:N,1:M);
x0 = median(xg(valid)); y0 = median(yg(valid));
phi_plane = (xg-x0)*ax + (yg-y0)*by;

% ---------------- assemble normal equations ----------------
idx = zeros(M,N); idx(valid) = 1:nnz(valid);  n = nnz(valid);
I=[]; J=[]; V=[]; rhs = zeros(n,1);

% horizontal edges
[rh,ch] = find(ex);
p = idx(sub2ind([M,N], rh,   ch));
q = idx(sub2ind([M,N], rh,   ch+1));
w = Wx(ex); d = gx_t(ex);
I = [I; p; q; p; q]; J = [J; p; q; q; p]; V = [V; w; w; -w; -w];
rhs = rhs + accumarray([p; q], [-w.*d; +w.*d], [n,1], @sum, 0);

% vertical edges
[rv,cv] = find(ey);
p = idx(sub2ind([M,N], rv,   cv));
q = idx(sub2ind([M,N], rv+1, cv));
w = Wy(ey); d = gy_t(ey);
I = [I; p; q; p; q]; J = [J; p; q; q; p]; V = [V; w; w; -w; -w];
rhs = rhs + accumarray([p; q], [-w.*d; +w.*d], [n,1], @sum, 0);

Nmat = sparse(I,J,V,n,n);
Nmat = Nmat + lambda*speye(n) + 1e-6*speye(n);

% ---------------- solve & combine with plane ----------------
try
    S.type='ict'; S.droptol=1e-3;
    L = ichol(Nmat,S);
    [phi,flag] = pcg(Nmat, rhs, 1e-4, 300, L, L');
    if flag~=0, phi = Nmat\rhs; end
catch
    phi = Nmat\rhs;
end
phi_img = NaN(M,N); phi_img(valid)=phi;
mphi = median(phi(isfinite(phi)),'omitnan'); 
phi_img = phi_img - mphi;

phi_img = phi_img + phi_plane;
phi_img = nanGaussLP_aniso(phi_img, 0.5*sigmaY, 0.5*sigmaX); % gentle

I_corr = I_filt .* exp(-1i*phi_img);

% bookkeeping
model.phi_pred = phi_img;
model.ax = ax; model.by = by;
model.sigmaX=sigmaX; model.sigmaY=sigmaY; model.lambda=lambda;
model.beta=beta; model.p=p;
end
function Gs = circSmoothAniso(G, W, sigY, sigX)
    % smooth Re/Im with separable Gaussians; normalize by smoothed weights
    R=zeros(size(G)); I=zeros(size(G)); M=zeros(size(G));
    v = isfinite(G)&isfinite(W)&(W>0);
    R(v)=real(G(v)).*W(v); I(v)=imag(G(v)).*W(v); M(v)=W(v);
    Rf = nanGaussLP_aniso(R,sigY,sigX);
    If = nanGaussLP_aniso(I,sigY,sigX);
    Mf = nanGaussLP_aniso(M,sigY,sigX);
    Rm = Rf ./ max(Mf,eps); Im = If ./ max(Mf,eps);
    Gs = Rm + 1i*Im;
    mag = abs(Gs); mag(mag==0)=1; Gs = Gs ./ mag;
end

function Aout = nanGaussLP_aniso(A, sigY, sigX)
    if (sigY<=0 && sigX<=0) || all(~isfinite(A(:))), Aout=A; return; end
    Ly = max(3, ceil(5*sigY)); y = (-Ly:Ly);
    Lx = max(3, ceil(5*sigX)); x = (-Lx:Lx);
    gy = exp(-(y.^2)/(2*sigY^2)); gy = gy/sum(gy);
    gx = exp(-(x.^2)/(2*sigX^2)); gx = gx/sum(gx);
    mask = isfinite(A); A(~mask)=0;
    num = conv2(conv2(A, gy(:), 'same'), gx, 'same');
    den = conv2(conv2(double(mask), gy(:), 'same'), gx, 'same');
    Aout = num ./ max(den,eps); Aout(den==0)=NaN;
end

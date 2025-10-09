function [I_out, phi_total, logtab] = remove_phzramp_multiscale(Ifilt, coh, sigmas, lambda, cohThresh, tol, maxit)
% Multi-scale Poisson ramp removal with early stop.
% sigmas: vector, e.g. [120 90 70 60]
if nargin<3 || isempty(sigmas), sigmas = [120 90 70 60 50 40 30 25 15 5]; end
if nargin<4 || isempty(lambda), lambda = 1e-3; end
if nargin<5 || isempty(cohThresh), cohThresh = 0.30; end
if nargin<6 || isempty(tol), tol = 5e-3; end
if nargin<7 || isempty(maxit), maxit = numel(sigmas); end

I = Ifilt; phi_total = 0;
[M,N] = size(I);
valid = isfinite(I) & isfinite(coh) & (coh >= cohThresh);

logtab = struct('it',{},'sigma',{},'medgx',{},'medgy',{},'meangx',{},'meangy',{},'med_dphi',{});
prev_med_dphi = inf;

for it = 1:min(maxit, numel(sigmas))
    sg = sigmas(it);
    [I_next, mdl] = insar.remove_phzramp_poisson(I, coh, sg, lambda);

    % improvement metric: complex difference
    dphi = angle(I_next .* conj(I));         % (-pi,pi]
    med_dphi = median(abs(dphi(isfinite(dphi))), 'omitnan');

    % target gradients at this scale (same as inside the method)
    U = I ./ max(abs(I), eps);
    ex = valid(:,1:N-1) & valid(:,2:N);  ey = valid(1:M-1,:) & valid(2:M,:);
    Gx = NaN(M,N-1); Gy = NaN(M-1,N);
    tmpx = U(:,2:N).*conj(U(:,1:N-1)); Gx(ex)=tmpx(ex);
    tmpy = U(2:M,:).*conj(U(1:M-1,:)); Gy(ey)=tmpy(ey);
    Wx = NaN(M,N-1); Wy = NaN(M-1,N);
    tmpwx = sqrt(coh(:,2:N).*coh(:,1:N-1)); Wx(ex)=tmpwx(ex);
    tmpwy = sqrt(coh(2:M,:).*coh(1:M-1,:)); Wy(ey)=tmpwy(ey);
    Gx_s = circSmooth(Gx, Wx, sg);  Gy_s = circSmooth(Gy, Wy, sg);
    vX = isfinite(Gx_s)&isfinite(Wx)&(Wx>0);  vY = isfinite(Gy_s)&isfinite(Wy)&(Wy>0);
    gx_t = NaN(M,N-1); gy_t = NaN(M-1,N); gx_t(vX)=angle(Gx_s(vX)); gy_t(vY)=angle(Gy_s(vY));

    meangx = sum(Wx(vX).*gx_t(vX)) / max(sum(Wx(vX)), eps);
    meangy = sum(Wy(vY).*gy_t(vY)) / max(sum(Wy(vY)), eps);
    medgx  = median(abs(gx_t(vX)), 'omitnan');
    medgy  = median(abs(gy_t(vY)), 'omitnan');

    logtab(it).it = it; logtab(it).sigma = sg;
    logtab(it).medgx = medgx; logtab(it).medgy = medgy;
    logtab(it).meangx = meangx; logtab(it).meangy = meangy;
    logtab(it).med_dphi = med_dphi;

    phi_total = phi_total + mdl.phi_pred;
    I = I_next;

    % early stop: absolute target + diminishing returns
    if abs(meangx) < 2e-5 && abs(meangy) < 2e-5 && ...
       (prev_med_dphi - med_dphi) < tol
        break;
    end
    prev_med_dphi = med_dphi;
end

I_out = I;
end
function Gs = circSmooth(G, W, sigma)
    R = zeros(size(G)); I = zeros(size(G)); M = zeros(size(G));
    v = isfinite(G) & isfinite(W) & (W>0);
    R(v) = real(G(v)) .* W(v);
    I(v) = imag(G(v)) .* W(v);
    M(v) = W(v);
    Rf = nanGaussLP(R, sigma); If = nanGaussLP(I, sigma); Mf = nanGaussLP(M, sigma);
    Rm = Rf ./ max(Mf, eps); Im = If ./ max(Mf, eps);
    Gs = Rm + 1i*Im;
    mag = abs(Gs); mag(mag==0) = 1;
    Gs = Gs ./ mag;
end

function Aout = nanGaussLP(A, sigma)
    if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
    L=max(3,ceil(5*sigma)); x=-L:L; g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    mask=isfinite(A); A(~mask)=0;
    num=conv2(conv2(A,g,'same'),g','same');
    den=conv2(conv2(double(mask),g,'same'),g','same');
    Aout=num./max(den,eps); Aout(den==0)=NaN;
end
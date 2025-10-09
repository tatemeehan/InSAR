function [I_corr, model] = remove_phzramp_poisson(I_filt, coh, sigmaGrad, lambda)
% Coherence-weighted smooth phase screen (no external predictors).
% I_filt  : complex filtered interferogram (MxN)
% coh     : coherence (MxN)
% sigmaGrad : Gaussian sigma in pixels to smooth complex gradients (50–150 typical)
% lambda  : small Laplacian ridge (e.g., 1e-3)

    if nargin<3 || isempty(sigmaGrad), sigmaGrad = 25; end
    if nargin<4 || isempty(lambda),    lambda    = 1e-3; end

    [M,N] = size(I_filt);
    valid = isfinite(I_filt) & isfinite(coh) & (coh>=0.3);

    % ---- unit phasor ----
    U = I_filt ./ max(abs(I_filt), eps);

    % ---- edge masks ----
    ex = valid(:,1:N-1) & valid(:,2:N);     % horizontal edges
    ey = valid(1:M-1,:) & valid(2:M,:);     % vertical edges

    % ---- complex edge phasors and weights (no chained indexing) ----
    Gx = NaN(M, N-1);  Gy = NaN(M-1, N);
    tmp = U(:,2:N) .* conj(U(:,1:N-1));     Gx(ex) = tmp(ex);
    tmp = U(2:M,:)  .* conj(U(1:M-1,:));    Gy(ey) = tmp(ey);

    Wx = NaN(M, N-1); Wy = NaN(M-1, N);
    tmp = sqrt(coh(:,2:N)  .* coh(:,1:N-1)); Wx(ex) = tmp(ex);
    tmp = sqrt(coh(2:M,:)  .* coh(1:M-1,:)); Wy(ey) = tmp(ey);

    % ---- smooth complex gradients by vector averaging ----
    Gx_s = circSmooth(Gx, Wx, sigmaGrad);
    Gy_s = circSmooth(Gy, Wy, sigmaGrad);

    % target edge gradients (radians)
    gx_t = NaN(M,N-1); gy_t = NaN(M-1,N);
    vX = isfinite(Gx_s) & isfinite(Wx) & (Wx>0);
    vY = isfinite(Gy_s) & isfinite(Wy) & (Wy>0);
    gx_t(vX) = angle(Gx_s(vX));
    gy_t(vY) = angle(Gy_s(vY));
    fprintf('med|gx|=%.4f, med|gy|=%.4f (rad)\n', ...
    median(abs(gx_t(isfinite(gx_t))),"omitnan"), ...
    median(abs(gy_t(isfinite(gy_t))),"omitnan"));

    % --- mean-gradient plane anchor (uses same edge weights) ---
    ax = 0; by = 0;
    if any(vX(:))
        ax = sum(Wx(vX) .* gx_t(vX)) / max(sum(Wx(vX)), eps);   % rad per pixel along +x
    end
    if any(vY(:))
        by = sum(Wy(vY) .* gy_t(vY)) / max(sum(Wy(vY)), eps);   % rad per pixel along +y
    end

    % centered coordinates so plane has ~zero mean
    [xg, yg] = meshgrid(1:N, 1:M);
    x0 = median(xg(valid));   % robust center on valid set
    y0 = median(yg(valid));
    phi_plane = (xg - x0) * ax + (yg - y0) * by;

    mx = sum(Wx(vX) .* gx_t(vX)) / max(sum(Wx(vX)),eps);
    my = sum(Wy(vY) .* gy_t(vY)) / max(sum(Wy(vY)),eps);
    fprintf('mean gx=%.4g, mean gy=%.4g rad/px\n', mx, my);




    % ---- build index of valid pixels ----
    idx = zeros(M,N); idx(valid) = 1:nnz(valid);
    n = nnz(valid);

    % ---- assemble weighted Laplacian normal equations directly ----
    % Each edge contributes: w*(phi_q - phi_p - d)^2
    % -> matrix: +w at (p,p) and (q,q), -w at (p,q) and (q,p)
    % -> rhs: -w*d at p, +w*d at q
    I = []; J = []; V = [];
    rhs = zeros(n,1);

    % horizontal
    [r,c] = find(ex);
    p = idx(sub2ind([M,N], r,   c));
    q = idx(sub2ind([M,N], r,   c+1));
    w = Wx(ex); d = gx_t(ex);

    I = [I; p; q; p; q];
    J = [J; p; q; q; p];
    V = [V; w; w; -w; -w];
    rhs = rhs + accumarray([p; q], [-w.*d; +w.*d], [n,1], @sum, 0);

    % vertical
    [r,c] = find(ey);
    p = idx(sub2ind([M,N], r,   c));
    q = idx(sub2ind([M,N], r+1, c));
    w = Wy(ey); d = gy_t(ey);

    I = [I; p; q; p; q];
    J = [J; p; q; q; p];
    V = [V; w; w; -w; -w];
    rhs = rhs + accumarray([p; q], [-w.*d; +w.*d], [n,1], @sum, 0);

    % sparse system
    Nmat = sparse(I, J, V, n, n);

    % small ridge & gauge fix (zero-mean)
    Nmat = Nmat + lambda*speye(n) + 1e-6*speye(n);

    % ---- solve (PCG with incomplete Cholesky) ----
    try
        setup.type = 'ict';  setup.droptol = 1e-3;
        L = ichol(Nmat, setup);
        [phi,flag] = pcg(Nmat, rhs, 1e-4, 300, L, L');
        if flag ~= 0, phi = Nmat \ rhs; end
    catch
        % if ichol not available, fall back
        phi = Nmat \ rhs;
    end

    % put on grid, subtract mean to set gauge
    phi_img = NaN(M,N); phi_img(valid) = phi;
    mphi = median(phi(isfinite(phi)), 'omitnan'); 
    phi_img = phi_img - mphi;

    phi_img = phi_img + phi_plane;
    % gentle LP (optional)
    phi_img = nanGaussLP(phi_img, 0.5*sigmaGrad);

    % remove ramp in complex domain
    I_corr = I_filt .* exp(-1i * phi_img);

    % output
    model.phi_pred = phi_img;
    model.valid    = valid;
    model.sigmaGrad= sigmaGrad;
    model.lambda   = lambda;
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

% function [I_corr, model] = remove_phzramp_poisson(I_filt, coh, sigmaGrad, lambda)
% % Coherence-weighted smooth phase screen from data (no predictors).
% % sigmaGrad: Gaussian (px) to smooth the complex gradients (50–120 typical)
% % lambda   : small Laplacian ridge for stability (e.g., 1e-3)
% 
% if nargin<3 || isempty(sigmaGrad), sigmaGrad = 80; end
% if nargin<4 || isempty(lambda),    lambda    = 1e-3; end
% 
% [M,N] = size(I_filt);
% valid = isfinite(I_filt) & isfinite(coh) & (coh>0);
% 
% % Unit phasor
% U = I_filt ./ max(abs(I_filt), eps);
% 
% % Complex “gradient phasors”
% Gx = NaN(M,N-1); gyW = NaN(M-1,N);  % (we keep sizes of edge grids)
% Wx = NaN(M,N-1); Wy = NaN(M-1,N);
% 
% ix = valid(:,1:N-1) & valid(:,2:N);
% iy = valid(1:M-1,:) & valid(2:M,:);
% 
% Gx(ix) = U(:,2:N)(ix) .* conj(U(:,1:N-1)(ix));
% Gy(iy) = U(2:M,:)(iy) .* conj(U(1:M-1,:)(iy));
% 
% % Edge weights from coherence (geometric mean)
% Wx(ix) = sqrt(coh(:,1:N-1)(ix) .* coh(:,2:N)(ix));
% Wy(iy) = sqrt(coh(1:M-1,:)(iy) .* coh(2:M,:)(iy));
% 
% % Smooth complex gradients (vector average on the unit circle)
% Gx_s = circSmooth(Gx, Wx, sigmaGrad);
% Gy_s = circSmooth(Gy, Wy, sigmaGrad);
% 
% % Target gradient (radians on edges)
% gx_t = angle(Gx_s);   gy_t = angle(Gy_s);
% 
% % Build masked sparse gradient operator A and weights w
% [idx, nPix] = buildIndex(valid);
% [A, b, w] = buildGradientSystem(gx_t, gy_t, Wx, Wy, valid, idx);
% 
% % Optional Laplacian ridge
% L = buildLaplacian(valid, idx);    % nPix x nPix
% W = spdiags(w,0,length(w),length(w));
% Nmat = A.'*W*A + lambda*L; rhs = A.'*(W*b);
% 
% % Fix gauge (zero-mean) to remove nullspace
% Nmat = Nmat + (1e-6)*speye(nPix);
% phi = zeros(nPix,1);
% phi = Nmat \ rhs;
% 
% % Put back on grid and smooth lightly
% phi_img = NaN(M,N); phi_img(valid) = phi;
% phi_img = nanGaussLP(phi_img, 0.5*sigmaGrad);
% 
% % Remove in complex domain
% I_corr = I_filt .* exp(-1i*phi_img);
% 
% model.phi_pred = phi_img;
% model.valid    = valid;
% end
% 
% function Gs = circSmooth(G, W, sigma)
%     % Vector-average smoothing of unit complex with NaNs/weights
%     R = real(G); I = imag(G);
%     R(isnan(R)) = 0; I(isnan(I)) = 0; W(isnan(W)) = 0;
%     Rf = nanGaussLP(R.*W, sigma); If = nanGaussLP(I.*W, sigma);
%     Wf = nanGaussLP(W, sigma);
%     Rm = Rf ./ max(Wf, eps); Im = If ./ max(Wf, eps);
%     Gs = Rm + 1i*Im;
%     Gs = Gs ./ max(abs(Gs), eps);   % re-normalize
% end
% 
% function [idx, nPix] = buildIndex(valid)
%     idx = zeros(size(valid)); idx(valid) = 1:nnz(valid);
%     nPix = nnz(valid);
% end
% 
% function [A, b, w] = buildGradientSystem(gx, gy, Wx, Wy, valid, idx)
%     [M,N] = size(valid);
%     rows = []; cols = []; vals = [];
%     b = []; w = [];
% 
%     % Horizontal edges
%     ex = valid(:,1:N-1) & valid(:,2:N) & isfinite(gx) & isfinite(Wx);
%     [r,c] = find(ex);
%     p = idx(sub2ind([M,N], r, c));
%     q = idx(sub2ind([M,N], r, c+1));
%     m = numel(p);
%     rows = [rows, (1:m),      (1:m)];
%     cols = [cols, p(:).',     q(:).'];
%     vals = [vals, -ones(1,m),  ones(1,m)];
%     b    = [b; gx(ex)];
%     w    = [w; Wx(ex)];
% 
%     % Vertical edges
%     ey = valid(1:M-1,:) & valid(2:M,:) & isfinite(gy) & isfinite(Wy);
%     [r,c] = find(ey);
%     p = idx(sub2ind([M,N], r,   c));
%     q = idx(sub2ind([M,N], r+1, c));
%     k = numel(p);
%     rows = [rows, (m+1:m+k),     (m+1:m+k)];
%     cols = [cols, p(:).',        q(:).'];
%     vals = [vals, -ones(1,k),     ones(1,k)];
%     b    = [b; gy(ey)];
%     w    = [w; Wy(ey)];
% 
%     A = sparse(rows, cols, vals, numel(b), nnz(valid));
%     w = max(w, 0); w(~isfinite(w)) = 0;
%     b(~isfinite(b)) = 0;
% end
% 
% function L = buildLaplacian(valid, idx)
%     [M,N] = size(valid);
%     ii = []; jj = []; vv = [];
%     % 4-neighbour Laplacian on the valid graph
%     for d = [0 1; 0 -1; 1 0; -1 0].'
%         di = d(1); dj = d(2);
%         s = valid & circshift(valid,[-di,-dj]);
%         p = idx(s);
%         q = idx(circshift(s,[-di,-dj]));
%         ii = [ii, p, p]; jj = [jj, p, q];
%         vv = [vv, 1, -1];
%     end
%     n = nnz(valid);
%     L = sparse(ii, jj, vv, n, n);
% end
% 
% function Aout = nanGaussLP(A, sigma)
%     if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
%     L=max(3,ceil(5*sigma)); x=-L:L; g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
%     mask=isfinite(A); A(~mask)=0;
%     num=conv2(conv2(A,g,'same'),g','same');
%     den=conv2(conv2(double(mask),g,'same'),g','same');
%     Aout=num./max(den,eps); Aout(den==0)=NaN;
% end

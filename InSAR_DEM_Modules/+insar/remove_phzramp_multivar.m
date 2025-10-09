function [I_corr, model] = remove_phzramp_multivar(I_filt, coh, geom, rr)
% NaN-robust ramp removal with gradient-domain regression.
% I_filt : complex filtered interferogram (MxN)
% coh    : coherence (MxN)
% geom   : predictors on same grid (any of: dR, dInc, Bperp, X, Y)
% rr     : options struct

    % ---- defaults ----
    if nargin < 4, rr = struct(); end
    if ~isfield(rr,'use'),        rr.use        = {'dR','dInc','Bperp','X','Y'}; end
    if ~isfield(rr,'cohThresh'),  rr.cohThresh  = 0.35; end
    if ~isfield(rr,'sigmaPred'),  rr.sigmaPred  = 100; end
    if ~isfield(rr,'sigmaPhi'),   rr.sigmaPhi   = 50; end
    if ~isfield(rr,'robust'),     rr.robust     = true; end
    if ~isfield(rr,'lambda'),     rr.lambda     = 1e-3;    end   % small ridge like 1e-3 optional
    if ~isfield(rr,'polyXY'),   rr.polyXY = 2; end   % 1=planar X,Y; 2=add X^2, XY, Y^2


    [M,N] = size(I_filt);

    % --- coherence-based valid mask (use it for standardization AND fitting)
    validPix = isfinite(coh) & (coh >= rr.cohThresh) & isfinite(I_filt);

    % --- Ensure X/Y exist if requested
    if any(strcmp(rr.use,'X')) && ~isfield(geom,'X'), [geom.X,~] = meshgrid(1:N,1:M); end
    if any(strcmp(rr.use,'Y')) && ~isfield(geom,'Y'), [~,geom.Y] = meshgrid(1:N,1:M); end

    % --- Build NaN-robust low-passed predictors (ramp-scale)
    names = {};
    P = struct();
    for k=1:numel(rr.use)
        nm = rr.use{k};
        if isfield(geom,nm) && ~isempty(geom.(nm))
            Pk = nanGaussLP(geom.(nm), rr.sigmaPred);
            if any(isfinite(Pk(:)))
                P.(nm) = Pk; names{end+1} = nm; %#ok<AGROW>
            end
        end
    end

    % Early exit if nothing to regress
    if isempty(names)
        I_corr = I_filt;
        model  = struct('names',{},{},'coeffs',[],'phi_pred',zeros(M,N),'mask',validPix);
        return
    end

    % --- Add polynomial XY terms (before standardization) ---
    if rr.polyXY >= 1
        % ensure X,Y present in P; if user didn't request them, create from coord grid
        if ~isfield(P,'X'), [P.X,~] = meshgrid(1:N,1:M); names{end+1}='X'; end
        if ~isfield(P,'Y'), [~,P.Y] = meshgrid(1:N,1:M); names{end+1}='Y'; end
    end
    if rr.polyXY >= 2
        P.X2 = P.X.^2;   names{end+1}='X2';
        P.XY = P.X.*P.Y; names{end+1}='XY';
        P.Y2 = P.Y.^2;   names{end+1}='Y2';
    end

    % --- STANDARDIZE predictors *using validPix*
    [P, names, scale] = standardizePredictors(P, names, validPix);

    % --- Phase gradients (wrap-safe)
    [gx,gy] = phaseGradFinite(I_filt);

    % --- Predictor gradients
    Xx = []; Xy = [];
    for k=1:numel(names)
        [px,py] = gradsFinite(P.(names{k}));
        Xx = [Xx, px(:)]; %#ok<AGROW>
        Xy = [Xy, py(:)]; %#ok<AGROW>
    end

    % --- Masks for each gradient component
    maskX = validPix & isfinite(gx);
    maskY = validPix & isfinite(gy);

    % --- Stack observations (x; y) and design (x; y)
    y = [gx(maskX); gy(maskY)];
    X = [Xx(maskX,:); Xy(maskY,:)];

    % Guard: enough rows to solve?
    if isempty(y) || size(X,1) < size(X,2)
        I_corr = I_filt;
        model  = struct('names',{names},'coeffs',zeros(1,numel(names)), ...
                        'phi_pred',zeros(M,N),'mask',validPix,'scale',scale);
        return
    end

    % --- Weights (sqrt coherence), aligned with y
    w = [coh(maskX); coh(maskY)];
    w = sqrt(max(w,0));  w(~isfinite(w)) = 0;

    % --- Solve (robust or LS) with optional tiny ridge
    if rr.robust
        c = robustIRLS_rowsafe_ridge(X, y, w, rr.lambda);
    else
        c = ridgeSolve(X, y, w, rr.lambda);
    end
    if any(~isfinite(c)), c(:) = 0; end

    % --- Predicted ramp (combine standardized predictors) and smooth lightly
    phi_pred = zeros(M,N);
    for k=1:numel(names)
        phi_pred = phi_pred + c(k) * P.(names{k});
    end
    phi_pred = nanGaussLP(phi_pred, rr.sigmaPhi);

    % --- Subtract in complex domain
    I_corr = I_filt .* exp(-1i * phi_pred);

    model.names   = names;
    model.coeffs  = c(:).';
    model.phi_pred= phi_pred;
    model.mask    = validPix;
    model.scale   = scale;
end

% ---------- helpers ----------

function [Pstd, names, stats] = standardizePredictors(P, names, mask)
Pstd = struct(); stats = struct();
for k = 1:numel(names)
    nm = names{k}; A = P.(nm);
    v = isfinite(A) & mask;
    if ~any(v(:)), continue; end
    mu = median(A(v)); s = 1.4826*mad(A(v),1);
    if ~isfinite(s) || s < eps, s = std(A(v)); end
    if ~isfinite(s) || s < eps, s = 1; end
    Astd = (A - mu) / s;
    Pstd.(nm) = Astd;
    stats.(nm).mu = mu; stats.(nm).sigma = s;
end
end

function [gx,gy] = phaseGradFinite(I)
    [M,N] = size(I);
    U = I ./ max(abs(I), eps);
    gx = NaN(M,N); gy = NaN(M,N);
    fin = isfinite(U);
    finX = fin(:,1:N-1) & fin(:,2:N);
    Vx   = U(:,2:N) .* conj(U(:,1:N-1));
    Gx = NaN(M,N-1); Gx(finX) = angle(Vx(finX)); gx(:,1:N-1) = Gx;
    finY = fin(1:M-1,:) & fin(2:M,:);
    Vy   = U(2:M,:) .* conj(U(1:M-1,:));
    Gy = NaN(M-1,N); Gy(finY) = angle(Vy(finY)); gy(1:M-1,:) = Gy;
end

function [gx,gy] = gradsFinite(A)
    [M,N] = size(A);
    gx = NaN(M,N); gy = NaN(M,N); fin = isfinite(A);
    if N>=3
        finXi = fin(:,1:N-2) & fin(:,2:N-1) & fin(:,3:N);
        Dx = 0.5*(A(:,3:N)-A(:,1:N-2)); T=NaN(M,N-2); T(finXi)=Dx(finXi); gx(:,2:N-1)=T;
    end
    if N>=2
        maskL = fin(:,1)&fin(:,2); dl=NaN(M,1); tmp=A(:,2)-A(:,1); dl(maskL)=tmp(maskL); gx(:,1)=dl;
        maskR = fin(:,N-1)&fin(:,N); dr=NaN(M,1); tmp=A(:,N)-A(:,N-1); dr(maskR)=tmp(maskR); gx(:,N)=dr;
    end
    if M>=3
        finYi = fin(1:M-2,:)&fin(2:M-1,:)&fin(3:M,:);
        Dy = 0.5*(A(3:M,:)-A(1:M-2,:)); T=NaN(M-2,N); T(finYi)=Dy(finYi); gy(2:M-1,:)=T;
    end
    if M>=2
        maskT=fin(1,:)&fin(2,:); dt=NaN(1,N); tmp=A(2,:)-A(1,:); dt(maskT)=tmp(maskT); gy(1,:)=dt;
        maskB=fin(M-1,:)&fin(M,:); db=NaN(1,N); tmp=A(M,:)-A(M-1,:); db(maskB)=tmp(maskB); gy(M,:)=db;
    end
end

function Aout = nanGaussLP(A,sigma)
    if sigma<=0 || all(~isfinite(A(:))), Aout=A; return; end
    L=max(3,ceil(5*sigma)); x=(-L:L); g=exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    mask=isfinite(A); A(~mask)=0;
    num=conv2(conv2(A,g,'same'),g','same');
    den=conv2(conv2(double(mask),g,'same'),g','same');
    Aout=num./max(den,eps); Aout(den==0)=NaN;
end

function c = ridgeSolve(X,y,w,lambda)
    w = w(:); Xw = bsxfun(@times,X,w); yw = w.*y;
    if nargin<4 || isempty(lambda), lambda=0; end
    XtX = Xw.'*Xw; if lambda>0, XtX = XtX + lambda*eye(size(XtX)); end
    c = XtX \ (Xw.'*yw);
end

function c = robustIRLS_rowsafe_ridge(X,y,w,lambda)
    if nargin<4, lambda=0; end
    c = ridgeSolve(X,y,w,lambda);
    delta=1.345; maxit=30;
    for it=1:maxit
        r=y - X*c; s=1.4826*mad(r,1); if ~(isfinite(s)&&s>0), break; end
        u=r/s; psi=min(1,  delta./max(abs(u),eps));
        ww = w.*psi;
        c_new = ridgeSolve(X,y,ww,lambda);
        if any(~isfinite(c_new)), break; end
        if norm(X*(c_new-c),2) <= 1e-6*max(1,norm(y,2)), c=c_new; break; end
        c=c_new;
    end
end

% function [I_corr, model] = remove_phzramp_multivar(I_filt, coh, geom, rr)
% % NaN-robust ramp removal with gradient-domain regression.
% % I_filt : complex filtered interferogram (MxN)
% % coh    : coherence (MxN)
% % geom   : struct of predictors on same grid (any of: dR, dInc, Bperp, X, Y)
% % rr     : options
% 
%     % ---- defaults ----
%     if nargin < 4, rr = struct(); end
%     if ~isfield(rr,'use'),        rr.use        = {'dR','dInc','Bperp'}; end
%     if ~isfield(rr,'cohThresh'),  rr.cohThresh  = 0.35; end
%     if ~isfield(rr,'sigmaPred'),  rr.sigmaPred  = 60; end
%     if ~isfield(rr,'sigmaPhi'),   rr.sigmaPhi   = 40; end
%     if ~isfield(rr,'robust'),     rr.robust     = true; end
% 
%     [M,N] = size(I_filt);
%     phi = angle(I_filt);
% 
%     % Ensure X/Y exist if requested
%     if any(strcmp(rr.use,'X')) && ~isfield(geom,'X'), [geom.X,~] = meshgrid(1:N,1:M); end
%     if any(strcmp(rr.use,'Y')) && ~isfield(geom,'Y'), [~,geom.Y] = meshgrid(1:N,1:M); end
% 
%     % --- Build NaN-robust low-passed predictors ---
%     names = {};
%     P = struct();
%     for k=1:numel(rr.use)
%         nm = rr.use{k};
%         if isfield(geom,nm) && ~isempty(geom.(nm))
%             Pk = nanGaussLP(geom.(nm), rr.sigmaPred);   % NaN-robust LP
%             if ~all(~isfinite(Pk(:)))                   % skip if all NaN
%                 P.(nm) = Pk;
%                 names{end+1} = nm; %#ok<AGROW>
%             end
%         end
%     end
% 
%     [Pz, names, scale] = standardizePredictors(P, names, validPix);
%     P = Pz; model.scale = scale;
% 
%     % If no usable predictors, pass-through
%     if isempty(names)
%         I_corr = I_filt;
%         model = struct('names',{},{},'coeffs',[],'phi_pred',zeros(M,N),'mask',true(M,N));
%         return
%     end
% 
%     % --- Phase gradients (wrap-safe) ---
%     [gx,gy] = phaseGradFinite(I_filt);  % finite diffs only where both neighbors valid
% 
%     % --- Predictor gradients ---
%     Xx = []; Xy = [];
%     for k=1:numel(names)
%         [px,py] = gradsFinite(P.(names{k}));
%         Xx = [Xx, px(:)]; %#ok<AGROW>
%         Xy = [Xy, py(:)]; %#ok<AGROW>
%     end
% 
%     % --- Mask of valid observations ---
%     validPix = isfinite(coh) & (coh >= rr.cohThresh);
%     maskX = validPix & isfinite(gx);
%     maskY = validPix & isfinite(gy);
% 
%     y = [gx(maskX); gy(maskY)];
%     X = [Xx(maskX,:); Xy(maskY,:)];
% 
%     % Guard: enough observations?
%     if isempty(y) || size(X,1) < size(X,2)
%         % Not enough data to fit; return input unchanged
%         I_corr = I_filt;
%         model = struct('names',{names},'coeffs',zeros(1,numel(names)),'phi_pred',zeros(M,N),'mask',validPix);
%         return
%     end
% 
%     % Weights (sqrt coherence), aligned with y
%     w = [coh(maskX); coh(maskY)];
%     w = sqrt(max(w,0));
%     w(~isfinite(w)) = 0;
% 
%     % --- Fit coefficients (robust or LS) ---
%     if rr.robust
%         c = robustIRLS_rowsafe(X, y, w);
%     else
%         [Xw,yw] = applyWeights(X,y,w);
%         c = Xw \ yw;
%     end
%     if any(~isfinite(c)), c(:) = 0; end
% 
%     % --- Predicted ramp & smoothing (NaN-robust) ---
%     phi_pred = zeros(M,N);
%     for k=1:numel(names)
%         phi_pred = phi_pred + c(k) * P.(names{k});
%     end
%     phi_pred = nanGaussLP(phi_pred, rr.sigmaPhi);
% 
%     % --- Subtract in complex domain ---
%     I_corr = I_filt .* exp(-1i * phi_pred);
% 
%     model.names = names;
%     model.coeffs = c(:).';
%     model.phi_pred = phi_pred;
%     model.mask = validPix;
% end
% 
% % ---------- helpers ----------
% function Aout = nanGaussLP(A, sigma)
%     % NaN-aware Gaussian low-pass: conv(A*mask)/conv(mask)
%     if sigma <= 0 || all(~isfinite(A(:)))
%         Aout = A; return
%     end
%     L = max(3, ceil(5*sigma)); x = (-L:L);
%     g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
%     mask = isfinite(A);
%     A(~mask) = 0;
%     num = conv2(conv2(A, g, 'same'), g', 'same');
%     den = conv2(conv2(double(mask), g, 'same'), g', 'same');
%     Aout = num ./ max(den, eps);
%     Aout(den == 0) = NaN;
% end
% 
% function [gx,gy] = phaseGradFinite(I)
% % Wrap-safe phase gradients from complex interferogram I (MxN).
% % Uses unit phasors and guards against NaNs; no chained indexing.
% 
%     [M,N] = size(I);
%     U = I ./ max(abs(I), eps);   % unit phasors; avoid /0
% 
%     gx = NaN(M,N);
%     gy = NaN(M,N);
% 
%     fin = isfinite(U);
% 
%     % ---- x-direction (size M x (N-1)) ----
%     finX = fin(:,1:N-1) & fin(:,2:N);
%     Vx   = U(:,2:N) .* conj(U(:,1:N-1));   % M x (N-1)
%     Gx   = NaN(M, N-1);
%     Gx(finX) = angle(Vx(finX));
%     gx(:,1:N-1) = Gx;
% 
%     % ---- y-direction (size (M-1) x N) ----
%     finY = fin(1:M-1,:) & fin(2:M,:);
%     Vy   = U(2:M,:) .* conj(U(1:M-1,:));   % (M-1) x N   <-- fixed
%     Gy   = NaN(M-1, N);
%     Gy(finY) = angle(Vy(finY));
%     gy(1:M-1,:) = Gy;
% end
% 
% 
% 
% function [gx,gy] = gradsFinite(A)
% % Central differences where neighbors are finite; guarded edges; no chained indexing
%     [M,N] = size(A);
%     gx = NaN(M,N);
%     gy = NaN(M,N);
%     fin = isfinite(A);
% 
%     % ----- X gradients -----
%     % interior columns 2..N-1 (need j-1, j+1)
%     if N >= 3
%         finXi = fin(:,1:N-2) & fin(:,3:N) & fin(:,2:N-1);
%         Dx    = 0.5 * (A(:,3:N) - A(:,1:N-2));   % size M x (N-2)
%         T     = NaN(M, N-2);
%         T(finXi) = Dx(finXi);
%         gx(:,2:N-1) = T;
%     end
%     % left edge (col 1 needs 1 & 2)
%     if N >= 2
%         maskL = fin(:,1) & fin(:,2);
%         dl    = NaN(M,1);
%         tmpL  = A(:,2) - A(:,1);
%         dl(maskL) = tmpL(maskL);
%         gx(:,1) = dl;
% 
%         % right edge (col N needs N-1 & N)
%         maskR = fin(:,N-1) & fin(:,N);
%         dr    = NaN(M,1);
%         tmpR  = A(:,N) - A(:,N-1);
%         dr(maskR) = tmpR(maskR);
%         gx(:,N) = dr;
%     end
% 
%     % ----- Y gradients -----
%     % interior rows 2..M-1 (need i-1, i+1)
%     if M >= 3
%         finYi = fin(1:M-2,:) & fin(3:M,:) & fin(2:M-1,:);
%         Dy    = 0.5 * (A(3:M,:) - A(1:M-2,:));   % size (M-2) x N
%         T2    = NaN(M-2, N);
%         T2(finYi) = Dy(finYi);
%         gy(2:M-1,:) = T2;
%     end
%     % top edge (row 1 needs 1 & 2)
%     if M >= 2
%         maskT = fin(1,:) & fin(2,:);
%         dt    = NaN(1,N);
%         tmpT  = A(2,:) - A(1,:);
%         dt(maskT) = tmpT(maskT);
%         gy(1,:) = dt;
% 
%         % bottom edge (row M needs M-1 & M)
%         maskB = fin(M-1,:) & fin(M,:);
%         db    = NaN(1,N);
%         tmpB  = A(M,:) - A(M-1,:);
%         db(maskB) = tmpB(maskB);
%         gy(M,:) = db;
%     end
% end
% 
% function [Xw,yw] = applyWeights(X,y,w)
%     w = w(:);
%     if ~isempty(X)
%         Xw = bsxfun(@times, X, w);
%     else
%         Xw = X;
%     end
%     yw = w .* y;
% end
% 
% function c = robustIRLS_rowsafe(X,y,w)
%     [Xw,yw] = applyWeights(X,y,w);
%     if isempty(yw), c = zeros(size(X,2),1); return; end
%     c = Xw \ yw;
% 
%     delta = 1.345; maxit = 30;
%     for it=1:maxit
%         r = y - X*c;
%         s = 1.4826*mad(r,1); if ~(isfinite(s) && s > 0), break; end
%         u = r / s;
%         psi = min(1, delta ./ max(abs(u), eps));
%         ww = w .* psi;
%         [Xww,yww] = applyWeights(X,y,ww);
%         c_new = Xww \ yww;
%         if any(~isfinite(c_new)), break; end
%         if norm(X*(c_new - c), 2) <= 1e-6 * max(1, norm(y,2)), c = c_new; break; end
%         c = c_new;
%     end
% end
% 
% function [Pstd, names, stats] = standardizePredictors(P, names, mask)
% % P is a struct of predictors already low-passed and co-registered
% % mask selects valid pixels (e.g., coh >= thresh)
% 
% Pstd = struct(); stats = struct();
% for k = 1:numel(names)
%     nm = names{k};
%     A  = P.(nm);
%     v  = isfinite(A) & mask;
%     if ~any(v(:)), continue; end
% 
%     % robust center/scale
%     mu  = median(A(v));
%     s   = 1.4826 * mad(A(v), 1);    % robust std; fallback:
%     if ~isfinite(s) || s < eps, s = std(A(v)); end
%     if ~isfinite(s) || s < eps, s = 1; end
% 
%     Astd = (A - mu) / s;
%     Pstd.(nm) = Astd;
% 
%     stats.(nm).mu = mu;
%     stats.(nm).sigma = s;
% end
% end
% 
% 
% % function [I_corr, model] = remove_phzramp_multivar(I_filt, coh, geom, rr)
% % % REMOVE_PHZRAMP_MULTIVAR  Remove low-frequency phase ramp before unwrapping.
% % %
% % % Inputs:
% % %   I_filt : filtered complex interferogram (what you unwrap), MxN
% % %   coh    : coherence map (same size), used for weights/mask
% % %   geom   : struct of predictors on same grid; any subset of:
% % %            .dR, .dInc, .Bperp, .X, .Y   (others ok too)
% % %   rr     : options struct (all optional)
% % %            .use        = {'dR','dInc','Bperp','X','Y'}
% % %            .cohThresh  = 0.35
% % %            .sigmaPred  = 60   % px, LPF predictors (scale-select the ramp)
% % %            .sigmaPhi   = 40   % px, LPF predicted phase before subtracting
% % %            .robust     = true % IRLS/Huber
% % %
% % % Outputs:
% % %   I_corr : complex interferogram with ramp removed
% % %   model  : .names, .coeffs, .phi_pred, .mask
% % 
% %     % ---- defaults ----
% %     if nargin < 4, rr = struct(); end
% %     if ~isfield(rr,'use'),        rr.use        = {'dR','dInc','Bperp','X','Y'}; end
% %     if ~isfield(rr,'cohThresh'),  rr.cohThresh  = 0.35; end
% %     if ~isfield(rr,'sigmaPred'),  rr.sigmaPred  = 60; end
% %     if ~isfield(rr,'sigmaPhi'),   rr.sigmaPhi   = 40; end
% %     if ~isfield(rr,'robust'),     rr.robust     = true; end
% % 
% %     [M,N] = size(I_filt);
% %     hasIPT = exist('imgaussfilt','file')==2;
% %     LP = @(A,s) hasIPT*imgaussfilt(A,s) + ~hasIPT*gaussSep(A,s);
% % 
% %     phi = angle(I_filt);
% % 
% %     % --- ensure X,Y predictors exist if requested ---
% %     if any(strcmp(rr.use,'X')) && ~isfield(geom,'X'), [geom.X,~] = meshgrid(1:N,1:M); end
% %     if any(strcmp(rr.use,'Y')) && ~isfield(geom,'Y'), [~,geom.Y] = meshgrid(1:N,1:M); end
% % 
% %     % --- low-pass predictors to ramp scale ---
% %     names = {};
% %     P = struct();
% %     for k=1:numel(rr.use)
% %         nm = rr.use{k};
% %         if isfield(geom,nm) && ~isempty(geom.(nm))
% %             P.(nm) = LP(geom.(nm), rr.sigmaPred);
% %             names{end+1} = nm; %#ok<AGROW>
% %         end
% %     end
% %     if isempty(names)
% %         % If nothing to regress, return input
% %         I_corr = I_filt; 
% %         model = struct('names',{},'coeffs',[],'phi_pred',zeros(M,N),'mask',true(M,N));
% %         return
% %     end
% % 
% %     % --- phase gradients via complex-product trick (wrap-safe) ---
% %     % gx ≈ angle( I .* conj(I shifted -x) ), gy likewise.
% %     [gx,gy] = phaseGrad(I_filt);
% % 
% %     % --- predictor gradients ---
% %     Xx = []; Xy = [];
% %     for k=1:numel(names)
% %         [px,py] = grads(P.(names{k}));
% %         Xx = [Xx, px(:)]; %#ok<AGROW>
% %         Xy = [Xy, py(:)]; %#ok<AGROW>
% %     end
% % 
% %     % --- vectorize & mask by coherence ---
% %     mask = isfinite(gx) & isfinite(gy) & isfinite(coh) & (coh >= rr.cohThresh);
% %     y = [gx(mask); gy(mask)];
% %     X = [Xx(mask,:); Xy(mask,:)];   % stacked x then y equations
% %     % Weight by coherence (square-root for stability)
% %     w = sqrt([coh(mask); coh(mask)]);
% % 
% %     % --- fit coefficients (shared across x/y by stacking) ---
% %     if rr.robust
% %         c = robustIRLS(X, y, w);
% %     else
% %         c = (X.*w) \ (y.*w);
% %     end
% % 
% %     % --- rebuild predicted ramp φ_pred = Σ c_k * P_k, smooth a bit ---
% %     phi_pred = zeros(M,N);
% %     for k=1:numel(names)
% %         phi_pred = phi_pred + c(k) * P.(names{k});
% %     end
% %     phi_pred = LP(phi_pred, rr.sigmaPhi);
% % 
% %     % --- subtract ramp in complex domain ---
% %     I_corr = I_filt .* exp(-1i * phi_pred);
% % 
% %     model.names = names;
% %     model.coeffs = c(:).';
% %     model.phi_pred = phi_pred;
% %     model.mask = mask;
% % 
% % end
% % 
% % % ===== helpers =====
% % function [gx,gy] = phaseGrad(I)
% %     I = I ./ max(abs(I), eps);
% %     gx = nan(size(I)); gy = nan(size(I));
% %     Ix = I(:,2:end)   .* conj(I(:,1:end-1));
% %     Iy = I(2:end,:)   .* conj(I(1:end-1,:));
% %     gx(:,1:end-1) = angle(Ix);
% %     gy(1:end-1,:) = angle(Iy);
% % end
% % function [gx,gy] = grads(A)
% %     gx = zeros(size(A)); gy = zeros(size(A));
% %     gx(:,2:end-1) = 0.5*(A(:,3:end)-A(:,1:end-2));
% %     gx(:,1)       = A(:,2)-A(:,1);
% %     gx(:,end)     = A(:,end)-A(:,end-1);
% %     gy(2:end-1,:) = 0.5*(A(3:end,:)-A(1:end-2,:));
% %     gy(1,:)       = A(2,:)-A(1,:);
% %     gy(end,:)     = A(end,:)-A(end-1,:);
% % end
% % function out = gaussSep(A,s)
% %     L = max(3,ceil(5*s)); x=(-L:L);
% %     g = exp(-(x.^2)/(2*s^2)); g=g/sum(g);
% %     out = conv2(conv2(A,g,'same'), g','same');
% % end
% % function c = robustIRLS(X,y,w)
% %     c = (X.*w) \ (y.*w);
% %     delta = 1.345; maxit=30;
% %     for it=1:maxit
% %         r = y - X*c; s = 1.4826*mad(r,1); if s<eps, break; end
% %         u = r/s; psi = min(1, delta./max(abs(u),eps));
% %         ww = w .* psi;
% %         c = (X.*ww) \ (y.*ww);
% %         if norm((X.*w)*c - (y.*w))/max(1,norm(y.*w)) < 1e-6, break; end
% %     end
% % end

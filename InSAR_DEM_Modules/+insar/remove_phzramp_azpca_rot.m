function [Icorr, model] = remove_phzramp_azpca_rot(Ifilt, coh, thetaDeg, opts)
% Rotate → azPCA → rotate-back ramp → subtract, with anti-clipping padding.
% Requires remove_phzramp_azpca.m on path.

if nargin < 4, opts = struct(); end
if ~isfield(opts,'cohThresh'),  opts.cohThresh  = 0.25;  end
if ~isfield(opts,'sigmaAzTilt'),opts.sigmaAzTilt= 60;   end  % smooth slopes across az

[M,N] = size(Ifilt);

if abs(thetaDeg) < 1e-6
    % no rotation needed
    [Icorr, model] = remove_phzramp_azpca(Ifilt, coh, opts);
    model.thetaDeg = 0; model.phi_pred_back = model.phi_pred;
    return
end

% 1) Compute padded canvas that safely contains the rotated image
[Mpad,Npad] = paddedSizeForRotation(M, N, thetaDeg);

% 2) Center-pad inputs
[If_pad, idxCrop] = centerPad(Ifilt, Mpad, Npad);
coh_pad           = centerPad(coh,   Mpad, Npad);

% optional protectMask support
if isfield(opts,'protectMask') && ~isempty(opts.protectMask)
    prot_pad = centerPad(opts.protectMask, Mpad, Npad) > 0;
    optsRot  = opts; optsRot.protectMask = prot_pad;
else
    optsRot  = opts;
end

% 3) Rotate padded arrays into azimuth-aligned frame (no clipping now)
[If_rot, maskI]  = rotLikePad(If_pad, -thetaDeg);
[coh_rot, maskC] = rotLikePad(coh_pad, -thetaDeg);
maskRot = maskI & maskC;
If_rot(~maskRot)  = NaN;
coh_rot(~maskRot) = NaN;

% % 4) Run azPCA in rotated frame
% opts = struct('preRotated', true, 'nearRangeSigma', 100, 'returnRotated', true);
% topt = struct('improveMin', 0.005, 'lambdaRG', 0.8, 'slopeMADCap', 2.5, 'sigmaAzTilt', 80,'cohThresh',0.3);
% roiMask_rot = isfinite(coh_rot);
% [Icorr_roi, touch] = touchup_range_tilt_roi(If_rot, coh_rot, 0, roiMask_rot, opts, topt);
% [Icorr_rot, modelRot] = insar.remove_phzramp_azpca(Icorr_roi, coh_rot, optsRot);
% 
% [Icorr_rot, modelRot] = insar.remove_phzramp_azpca(If_rot, coh_rot, optsRot);
% opts = struct('preRotated', true, 'nearRangeSigma', 250, 'returnRotated', true);
% topt = struct('improveMin', 0.005, 'lambdaRG', 0.8, 'slopeMADCap', 2.5, 'sigmaAzTilt', 80,'cohThresh',0.3);
% roiMask_rot = isfinite(coh_rot);
% [Icorr_roi, touch] = touchup_range_tilt_roi(Icorr_rot, coh_rot, 0, roiMask_rot, opts, topt);
% % [Icorr2, touch] = touchup_range_tilt_balanced(Icorr_rot, coh_rot,thetaDeg, opts);
% % [Icorr3, touch2] = touchup_range_tilt(Icorr_rot, coh_rot,opts);
% 
% 
% 
% % 5) Rotate predicted ramp back to padded original frame, then crop to M×N
% [phi_back_pad, mask_back] = rotLikePad(modelRot.phi_pred, +thetaDeg);
% phi_back_pad(~mask_back) = 0;
% [phi_back_pad_touch, mask_back_touch] = rotLikePad(touch.phi_lin, +thetaDeg);
% phi_back_pad_touch(~mask_back) = 0;
% 
% % crop to original extent
% phi_back = phi_back_pad(idxCrop.rows, idxCrop.cols);
% phi_back_touch = phi_back_pad_touch(idxCrop.rows, idxCrop.cols);
% 
% phi_pred = phi_back+phi_back_touch;
% eta = fminbnd(@(e) gradCostAzRot(Ifilt, coh, phi_back*e, thetaDeg, getOpt(opts,'cohThresh',0.25)), 0, 1.5);
% 
% % eta = fminbnd(@(e) gradCostAzRot(Ifilt, coh, phi_back*e, thetaDeg, getOpt(opts,'cohThresh',0.25)), 0, 1.5);
% eta1 = fminbnd(@(e) gradCostAzRot(Ifilt, coh, phi_back*e, thetaDeg, getOpt(opts,'cohThresh',0.25)), 0, 1.5);
% eta2 = fminbnd(@(e) gradCostAzRot(Ifilt, coh, phi_back_touch*e, thetaDeg, getOpt(opts,'cohThresh',0.25)), 0, 1.5);
% 
% % apply optimal scale
% % Icorr = Ifilt .* exp(-1i * (eta * phi_pred));
% Icorr = Ifilt .* exp(-1i * (eta * phi_back + touch.eta.*phi_back_touch));
% --- 0) Inputs ---
% If_rot, coh_rot are already rotated into true-az frame (thetaDeg applied upstream)
% roiMask_rot is the valid ROI in the rotated frame
optsTilt = struct('preRotated', true, 'nearRangeSigma', 100, 'returnRotated', true);
topt     = struct('improveMin', 0.005, 'lambdaRG', 0.8, 'slopeMADCap', 2.5, ...
                  'sigmaAzTilt', 80, 'cohThresh', 0.3);
roiMask_rot = isfinite(coh_rot);

% --- 1) Range-tilt stage (returns corrected residual and tilt ramp/eta) ---
[I1_rot, touch] = touchup_range_tilt_roi(If_rot, coh_rot, 0, roiMask_rot, optsTilt, topt);
% touch.phi_lin  : tilt ramp shape (in rotated frame)
% touch.eta      : fitted scale for tilt (already applied inside I1_rot)

% --- 2) PCA stage on the residual (no global eta here) ---
optsPCA = struct('K',3,'sigmaAzPrep',15,'sigmaR',25,'cohThresh',0.3);
[I2_rot, modelRot] = insar.remove_phzramp_azpca(I1_rot, coh_rot, optsPCA);
% modelRot.phi_pred : PCA ramp shape that was applied once (η=1 implicitly)

% OPTIONAL refinement: re-scale PCA with a small global η2 on the residual
% (use a weighted az-gradient LS on I1_rot, not on Ifilt)
eta2 = fminbnd(@(e) gradCostAzRot(I1_rot, coh_rot, modelRot.phi_pred*e, 0, getOpt(opts,'cohThresh',0.25)), 0, 1.5);
% [eta2, infoEta] = fit_eta_scalar(I1_rot, coh_rot, modelRot.phi_pred, ...
                   % struct('cohThresh',0.25, 'etaMax',2.5, 'bandSigmaLo',8)); % light HP to ignore DC

% If remove_phzramp_azpca already applied phi_pred once, adjust by (eta2-1):
I2_rot = I2_rot .* exp(-1i * (eta2 - 1) * modelRot.phi_pred);

% --- 3) Rotate each predicted ramp back and crop to original frame ---
[phi_back_pad_pca,   mask_back_pca]   = rotLikePad(modelRot.phi_pred, +thetaDeg);
[phi_back_pad_tilt,  mask_back_tilt]  = rotLikePad(touch.phi_lin,     +thetaDeg);
% IMPORTANT: use the correct mask for each
phi_back_pad_pca(~mask_back_pca)   = 0;
phi_back_pad_tilt(~mask_back_tilt) = 0;
phi_back_pca  = phi_back_pad_pca( idxCrop.rows, idxCrop.cols );
phi_back_tilt = phi_back_pad_tilt( idxCrop.rows, idxCrop.cols );

% Combined predicted phase (for logging/visuals only)
phi_pred = eta2 * phi_back_pca + touch.eta * phi_back_tilt;

% --- 4) Final corrected interferogram in original frame (apply from original Ifilt) ---
Icorr = Ifilt .* exp(-1i * (eta2 * phi_back_pca + touch.eta * phi_back_tilt));


% 6) Apply correction on ORIGINAL grid (single resample of ramp only)
% Icorr = Ifilt .* exp(-1i * phi_back);

% 7) Quick metric (original grid)
% gradA  = angle( Ifilt(:,2:end).*conj(Ifilt(:,1:end-1)) );
% gradAc = angle( Icorr(:,2:end).*conj(Icorr(:,1:end-1)) );
% mAz    = isfinite(gradA) & min(coh(:,1:end-1),coh(:,2:end)) >= getOpt(opts,'cohThresh',0.3);
% medAz0 = median(abs(gradA(mAz)),'omitnan');
% medAz1 = median(abs(gradAc(mAz)),'omitnan');
% fprintf('azPCA-ROT(θ=%.2f°): med|∂φ/∂Az| %.4f → %.4f rad/px (Keff=%d)\n', ...
%         thetaDeg, medAz0, medAz1, size(modelRot.modes,2));

% 8) Pack model
model = modelRot;
model.thetaDeg      = thetaDeg;
model.phi_pred = phi_pred;
model.phi_pred_tilt = phi_back_tilt;
model.phi_pred_pca = phi_back_pca;
model.tilt = touch;
model.eta_tilt = touch.eta; 
model.eta_pca = eta2;
end

% ---------- helpers ----------
function [eta, info] = fit_eta_scalar(Ifilt, coh, phiPred, opts)
% Weighted LS fit for a single scale η using azimuth gradients.
% Ifilt, phiPred, coh are MxN; opts fields (all optional):
%   .cohThresh = 0.25
%   .thetaDeg  = 0          % rotate to true-az frame if needed
%   .roiMask   = []         % MxN logical; restricts where we fit
%   .bandSigmaLo = 0        % az lowpass (px) for baseline; 0 disables
%   .bandSigmaHi = 0        % az highpass (px); if >0, removes DC/tilt
%   .etaMax = 1.0           % cap
% Returns eta and simple diagnostics.

    if nargin<4, opts = struct(); end
    cohThresh = getf(opts,'cohThresh',0.25);
    thetaDeg  = getf(opts,'thetaDeg',0);
    R         = getf(opts,'roiMask',[]);
    sLo       = getf(opts,'bandSigmaLo',0);
    sHi       = getf(opts,'bandSigmaHi',0);
    etaMax    = getf(opts,'etaMax',1.0);
    etaMin    = getf(opts,'etaMin',-etaMax);


    I = Ifilt; P = phiPred; C = coh;

    % rotate if requested
    if abs(thetaDeg) > 1e-6
        I = imrotate(I,  -thetaDeg, 'bilinear','crop');
        P = imrotate(P,  -thetaDeg, 'bilinear','crop');
        C = imrotate(C,  -thetaDeg, 'bilinear','crop');
        if ~isempty(R)
            R = imrotate(double(R), -thetaDeg, 'nearest','crop')>0.5;
        end
    end

    % azimuth gradients (wrap-safe for I; simple diff for phi)
    U   = I ./ max(abs(I),eps);
    gA  = angle(U(:,2:end) .* conj(U(:,1:end-1)));     % M x (N-1)
    gAp = diff(P,1,2);                                  % M x (N-1)

    % band-limit (optional): gaussian high-pass / low-pass along azimuth
    if sLo>0, gA  = gA  - gauss1d_nan(gA, sLo);  gAp = gAp - gauss1d_nan(gAp, sLo); end
    if sHi>0, gA  = gauss1d_nan(gA, sHi);        gAp = gauss1d_nan(gAp, sHi);       end

    % valid mask and weights
    Cpair = min(C(:,1:end-1), C(:,2:end));
    M = isfinite(gA) & isfinite(gAp) & (Cpair >= cohThresh);
    if ~isempty(R), M = M & R(:,1:end-1) & R(:,2:end); end

    w = sqrt(max(Cpair,0)); w(~M) = 0;

    y = gA(M);
    x = gAp(M);

    % guard
    if isempty(y) || all(abs(x) < eps)
        eta = 0; info = struct('usedPairs',0,'numel',numel(M),'ok',false); return;
    end

    % weighted LS closed-form: eta = <x,y>_w / <x,x>_w
    num = sum( (w(M) .* x) .* y );
    den = sum( (w(M) .* x) .* x );
    eta = num / max(den, eps);
    eta = min(max(eta, etaMin), etaMax);
    % cap to [0, etaMax]
    % if ~isfinite(eta), eta = 0; end
    % eta = max(0, min(eta, etaMax));
    

    % simple report
    info = struct;
    info.usedPairs = nnz(M);
    info.numel     = numel(M);
    info.etaRaw    = num/den;
    info.eta       = eta;
end

function A = gauss1d_nan(A, sigma)
% NaN-aware azimuth (dim 2) gaussian smoother
    if sigma<=0, return; end
    [M,N] = size(A);
    L = max(3, ceil(5*sigma));
    x = (-L:L); g = exp(-(x.^2)/(2*sigma^2)); g=g/sum(g);
    for i=1:M
        r = A(i,:); m = isfinite(r);
        r2 = r; r2(~m)=0; m2 = double(m);
        num = conv(r2, g, 'same'); den = conv(m2, g, 'same');
        rsm = num ./ max(den, eps); rsm(den==0) = NaN;
        A(i,:) = rsm;
    end
end
function [eta1, eta2, info] = fit_eta_two(Ifilt, coh, phi1, phi2, opts)
% Joint fit of η1*phi1 + η2*phi2 in azimuth-gradient domain, with
% weighted Gram-Schmidt to de-correlate phi2 from phi1.
    if nargin<5, opts=struct(); end
    [e1,info1] = fit_eta_scalar(Ifilt, coh, phi1, opts); % initial η1
    % Build gradients as in fit_eta_scalar
    thetaDeg  = getf(opts,'thetaDeg',0);
    cohThresh = getf(opts,'cohThresh',0.25);

    I = Ifilt; C = coh; P1 = phi1; P2 = phi2;
    if abs(thetaDeg)>1e-6
        I = imrotate(I,  -thetaDeg, 'bilinear','crop');
        C = imrotate(C,  -thetaDeg, 'bilinear','crop');
        P1= imrotate(P1, -thetaDeg, 'bilinear','crop');
        P2= imrotate(P2, -thetaDeg, 'bilinear','crop');
    end
    U = I ./ max(abs(I),eps);
    gA  = angle(U(:,2:end).*conj(U(:,1:end-1)));
    gP1 = diff(P1,1,2);
    gP2 = diff(P2,1,2);

    Cpair = min(C(:,1:end-1),C(:,2:end));
    M = isfinite(gA)&isfinite(gP1)&isfinite(gP2)&(Cpair>=cohThresh);
    w = sqrt(max(Cpair,0)); w(~M)=0;

    y = gA(M); x1 = gP1(M); x2 = gP2(M);

    % De-correlate x2 from x1 (weighted)
    a = sum((w(M).*x1).*x2) / max(sum((w(M).*x1).*x1), eps);
    x2o = x2 - a*x1;

    % Weighted ridge LS for [eta1 eta2]
    lam = getf(opts,'ridge',1e-3);
    Xw = [w(M).*x1, w(M).*x2o];
    yw = w(M).*y;
    H = Xw.'*Xw + lam*eye(2);
    b = H \ (Xw.'*yw);
    eta1 = max(0, b(1));
    eta2 = max(0, b(2));
    info = struct('proj',a,'eta1',eta1,'eta2',eta2,'ok',all(isfinite(b)));
end


function v = getf(s,f,d); if ~isstruct(s)||~isfield(s,f)||isempty(s.(f)), v=d; else, v=s.(f); end, end

function c = gradCostAzRot(Ifilt, coh, phi, thetaDeg, cohThresh)
    If = Ifilt .* exp(-1i*phi);
    % rotate to az-aligned frame (no clipping; same helper from the rot fn)
    [If_rot, maskI]  = rotLikePad(If, -thetaDeg);
    [coh_rot, maskC] = rotLikePad(coh, -thetaDeg);
    m = maskI & maskC & (coh_rot >= cohThresh);
    % azimuth gradient in rotated frame
    g = angle( If_rot(:,2:end).*conj(If_rot(:,1:end-1)) );
    m2 = m(:,1:end-1) & m(:,2:end);
    c = median(abs(g(m2)),'omitnan');   % objective
end

function [Mpad,Npad] = paddedSizeForRotation(M,N,thetaDeg)
% Size of bounding box of MxN after rotation by theta
th = abs(thetaDeg) * pi/180;
Mrot = ceil(abs(M*cos(th)) + abs(N*sin(th)));
Nrot = ceil(abs(M*sin(th)) + abs(N*cos(th)));
% guard both directions & add margin
Mpad = max([M, Mrot]) + 4;
Npad = max([N, Nrot]) + 4;
end

function [B, idx] = centerPad(A, Mpad, Npad)
% Center pad A (MxN) to (Mpad x Npad) with NaNs (or 0 for logical)
[M,N] = size(A);
top  = floor((Mpad - M)/2);
bot  = Mpad - M - top;
left = floor((Npad - N)/2);
right= Npad - N - left;

if islogical(A)
    padVal = false;
else
    padVal = NaN;
end
B = padarray(A, [top left], padVal, 'pre');
B = padarray(B, [bot right], padVal, 'post');

idx.rows = (top+1):(top+M);
idx.cols = (left+1):(left+N);
end

function [B, maskOut] = rotLikePad(A, thetaDeg)
% Robust rotation with NaN support; input size preserved ('crop')
% Works for real or complex A.
    if isreal(A)
        mask = isfinite(A);
        A2 = A; A2(~mask) = 0;
        B  = imrotate(A2, thetaDeg, 'bilinear', 'crop');
    else
        mask = isfinite(real(A)) & isfinite(imag(A));
        Ar = real(A); Ai = imag(A);
        Ar(~mask) = 0; Ai(~mask) = 0;
        Br = imrotate(Ar, thetaDeg, 'bilinear', 'crop');
        Bi = imrotate(Ai, thetaDeg, 'bilinear', 'crop');
        B  = complex(Br, Bi);
    end
    maskOut = imrotate(mask, thetaDeg, 'nearest', 'crop') > 0.5;
    B(~maskOut) = NaN;
end

function v = getOpt(s, name, default)
    if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name)), v = default; else, v = s.(name); end
end

function [Iout, info] = touchup_range_tilt_roi(Iin, cohIn, thetaDeg, roiMask, opts, touchupOpts)
% ROI wrapper for your existing touchup_range_tilt_balanced.m
%
% [Iout, info] = touchup_range_tilt_roi(Iin, cohIn, thetaDeg, roiMask, opts, touchupOpts)
%
% Inputs
%   Iin, cohIn : MxN complex interferogram and coherence (same frame)
%   thetaDeg   : rotation that puts scene into true-az frame.
%                If inputs are already rotated to true-az, set opts.preRotated=true
%   roiMask    : MxN logical mask in the SAME frame as (Iin, cohIn).
%                Pixels outside ROI are down-weighted (coherence->0).
%                Pass [] to disable ROI gating.
%   opts       : struct (all optional)
%       .preRotated      (false)  % true if Iin/cohIn already in true-az frame
%       .nearRangeSigma  ([])     % rows std-dev (px). If set, boosts near range
%       .returnRotated   ([])     % default: return in SAME frame as Iin
%   touchupOpts : passed straight into touchup_range_tilt_balanced (e.g. improveMin, etaMax, etc.)
%
% Outputs
%   Iout : corrected interferogram in the SAME frame as Iin
%   info : passthrough from base touch-up + wrapper metadata

if nargin < 5 || isempty(opts),        opts = struct(); end
if nargin < 6 || isempty(touchupOpts), touchupOpts = struct(); end
preRotated    = getf(opts,'preRotated',false);
nearSigma     = getf(opts,'nearRangeSigma',[]);
returnRotated = getf(opts,'returnRotated',[]);

% By default, return in SAME frame as input
if isempty(returnRotated), returnRotated = true; end

I   = Iin;
coh = cohIn;
R   = roiMask;

% 1) Rotate to true-az frame if needed (and rotate ROI mask too)
rotApplied = false;
if ~preRotated && abs(thetaDeg) > 1e-6
    [I,  M_I ] = imrotate_with_mask(Iin,  -thetaDeg);
    [coh,M_coh] = imrotate_with_mask(cohIn, -thetaDeg);
    if ~isempty(roiMask)
        R = imrotate(double(roiMask), -thetaDeg, 'nearest', 'crop') > 0.5;
    end
    % retain NaNs outside the valid rotated footprint
    I(~M_I)   = NaN;
    coh(~M_coh) = NaN;
    rotApplied = true;
else
    % ensure R aligns to current frame if provided
    if ~isempty(roiMask), R = logical(roiMask); end
end

% 2) Apply ROI gating by modulating coherence
if ~isempty(R)
    % Outside ROI, set coherence ~ 0 so base touch-up ignores it
    coh = coh .* double(R);
end

% 3) Optional near-range emphasis (rows Gaussian) via coherence weights
if ~isempty(nearSigma) && isfinite(nearSigma) && nearSigma > 0
    [M,N] = size(coh);
    r0 = nan(1,N);

    % Choose anchor per column: top-most ROI pixel (preferred), else top-most valid coherence
    for j = 1:N
        if ~isempty(R)
            idx = find(R(:,j) & isfinite(coh(:,j)), 1, 'first');
        else
            idx = find(isfinite(coh(:,j)), 1, 'first');
        end
        if ~isempty(idx), r0(j) = idx; end
    end

    % Build 2-D Gaussian weights centered at r0(j)
    [rr, cc] = ndgrid(1:M, 1:N);
    Wnr = zeros(M,N);
    goodCols = isfinite(r0);
    if any(goodCols)
        r0full = r0(cc);                               % broadcast r0 over rows
        Wnr = exp( - (rr - r0full).^2 / (2*nearSigma^2) );
        Wnr(:, ~goodCols) = 0;                         % no anchor → no boost
    end

    % Apply as multiplicative weight; keep NaNs outside ROI/footprint as-is
    coh = coh .* Wnr;

    % Optional: if this knocks too many pixels below your cohThresh,
    % consider lowering touchupOpts.cohThresh a bit (e.g., 0.25 → 0.20).
end
% if ~isempty(nearSigma)
%     [M,~] = size(I);
%     rr = (0:M-1)';
%     wR = exp(-(rr-0).^2 / (2*nearSigma^2));        % peak at top (near range)
%     coh = coh .* repmat(wR,1,size(coh,2));
% end

% 4) Call your existing balanced touch-up in this frame (no further rotation inside)
if ~exist('touchup_range_tilt_balanced','file')
    error('touchup_range_tilt_balanced.m not found on path.');
end
[IrotCorr, infoBase] = touchup_range_tilt_balanced(I, coh, 0, touchupOpts);

% 5) Return in the same frame as input
if rotApplied && ~preRotated && ~returnRotated
    % keep rotated frame (user asked to keep the rotated)
    Iout = IrotCorr;
else
    % rotate back if we rotated in step 1 and the caller wants original frame
    if rotApplied && returnRotated
        Iout = imrotate(IrotCorr, thetaDeg, 'bilinear', 'crop');
        % re-mask outside original footprint to NaN
        [~, M_back] = imrotate_with_mask(ones(size(Iin)), thetaDeg); %#ok<ASGLU>
        Iout(~M_back) = NaN;
    else
        Iout = IrotCorr;
    end
end

% 6) Collect info
info = infoBase;
info.wrapper.preRotated    = preRotated;
info.wrapper.nearSigma     = nearSigma;
info.wrapper.rotApplied    = rotApplied;
info.wrapper.returnRotated = returnRotated;
info.wrapper.roiCoverage   = nnz(R)/numel(R);
end

% ----------------- helpers -----------------
function [Ir, Mr] = imrotate_with_mask(I, theta)
% Legacy-friendly NaN-aware rotation using 'bilinear'+'crop'
% Mr marks the valid (non-padded) region after rotation.
sz = size(I);
baseMask = true(sz);
Ir = imrotate(I, theta, 'bilinear', 'crop');
Mr = imrotate(baseMask, theta, 'nearest',  'crop');
Ir(~Mr) = NaN;
end

% function v = getf(s, f, d)
% if ~isstruct(s) || ~isfield(s,f) || isempty(s.(f)), v = d; else, v = s.(f); end
% end

function [azP, rgMed] = metrics_rot_roi(I, coh, thetaDeg, cohThresh, roiMask, azBand)
Ir=I; cohr=coh;
if abs(thetaDeg) > 1e-6
    Ir   = imrotate(I,   -thetaDeg, 'bilinear', 'crop');
    cohr = imrotate(coh, -thetaDeg, 'bilinear', 'crop');
end
Ur  = Ir ./ max(abs(Ir), eps);
gAz = angle( Ur(:,2:end).*conj(Ur(:,1:end-1)) );
mAz = isfinite(gAz) & min(cohr(:,1:end-1),cohr(:,2:end))>=cohThresh;
if ~isempty(roiMask)
    mr = imrotate(roiMask, -thetaDeg, 'nearest', 'crop');
    mAz = mAz & mr(:,1:end-1) & mr(:,2:end);
end
if ~isempty(azBand)
    % band-limit along azimuth (columns)
    gAz = bandlimit_az(gAz, azBand);
end
azP = prctile(abs(gAz(mAz)), 90);       % p90 in ROI (and band)

U   = I ./ max(abs(I), eps);
gR  = angle( U(2:end,:).*conj(U(1:end-1,:)) );
mR  = isfinite(gR) & min(coh(1:end-1,:),coh(2:end,:))>=cohThresh;
rgMed = median(abs(gR(mR)), 'omitnan');
end

function J = cost_balanced_roi(I, coh, phi, thetaDeg, cohThresh, lambdaRG, roiMask, azBand)
[az0, rg0] = metrics_rot_roi(I,    coh, thetaDeg, cohThresh, roiMask, azBand);
[az1, rg1] = metrics_rot_roi(I.*exp(-1i*phi), coh, thetaDeg, cohThresh, roiMask, azBand);
J = az1 + lambdaRG * max(0, rg1 - rg0);
end

function G = bandlimit_az(G, band)
% Simple FFT bandpass along azimuth (columns) for each row of G
[M, N1] = size(G);
N = N1; 
F = fftshift(fft(G, [], 2), 2);
% frequency axis in cycles/pixel
f = linspace(-0.5, 0.5, N);
mask = (abs(f) >= band(1)) & (abs(f) <= band(2));
F = F .* mask; 
G = real(ifft(ifftshift(F,2), [], 2));
end

function [Icorr2, info] = touchup_range_tilt_balanced(Icorr, coh, thetaDeg, opts)
% Small, safe per-column range-tilt removal with automatic scaling ζ.
% Inputs:
%   Icorr     : complex interferogram after your az-PCA step (MxN)
%   coh       : coherence (MxN)
%   thetaDeg  : rotation angle used for az-metric (same one you liked)
%   opts      : struct (all optional)
%       .cohThresh   (0.45)   – mask for confident diffs
%       .sigmaAzTilt (60)     – az smoothing of column slopes
%       .slopeMADCap (3.0)    – cap |slope| ≤ cap*MAD of valid slopes
%       .lambdaRG    (0.5)    – penalty weight for ↑ range gradient
%       .improveMin  (0.03)   – required relative az-metric drop (3%)
%       .etaMax      (1.0)    – max ζ tried in line search

if nargin<4, opts=struct(); end
cohThresh   = getf(opts,'cohThresh', 0.45);
sigmaAzTilt = getf(opts,'sigmaAzTilt', 60);
slopeMADCap = getf(opts,'slopeMADCap', 2.5);
lambdaRG    = getf(opts,'lambdaRG', 0.5);
improveMin  = getf(opts,'improveMin', 0.005);
etaMax      = getf(opts,'etaMax', 1.0);

[M,N] = size(Icorr);
U  = Icorr ./ max(abs(Icorr), eps);
gR = angle( U(2:end,:).*conj(U(1:end-1,:)) );           % (M-1)×N residual range-grad
W  = min(coh(1:end-1,:), coh(2:end,:));
V  = isfinite(gR) & (W >= cohThresh);

% robust per-column slope (median of valid diffs)
slope = nan(1,N);
for j=1:N
    v = V(:,j);
    if nnz(v) >= 30
        slope(j) = median(gR(v,j), 'omitnan');
    end
end
validSlope = ~isnan(slope);
% cap extreme slopes by MAD to avoid over-tilting
slMAD = 1.4826*mad(slope(isfinite(slope)), 1);
cap   = slopeMADCap * max(slMAD, eps);
slope(validSlope) = max(min(slope(validSlope), cap), -cap);

% smooth slopes along azimuth
slope = nanconv1D(slope, sigmaAzTilt);

% build unit-ζ tilt phase surface φ_lin(ζ) = ζ * r * slope(j)
r = (0:M-1).';
phi_lin_1 = r .* slope;                              % M×N, NaNs propagate
% de-mean per column over confident pixels so we don’t inject bias
for j=1:N
    vv = isfinite(phi_lin_1(:,j)) & (coh(:,j) >= cohThresh);
    if any(vv), phi_lin_1(:,j) = phi_lin_1(:,j) - mean(phi_lin_1(vv,j)); end
end

% --- metrics BEFORE ---
[az0, rg0] = metrics_rot(Icorr, coh, thetaDeg, cohThresh);

% line-search ζ in [0, etaMax] to balance az improvement vs range penalty
cost = @(eta) cost_balanced(Icorr, coh, phi_lin_1*eta, thetaDeg, cohThresh, lambdaRG);
eta  = fminbnd(cost, 0, etaMax);

Icorr2 = Icorr .* exp(-1i * (phi_lin_1 * eta));

% --- metrics AFTER ---
[az1, rg1] = metrics_rot(Icorr2, coh, thetaDeg, cohThresh);

relImprove = (az0 - az1) / max(az0, eps);   % rotated-frame az metric

% optional guard: if not enough improvement, skip tilt
if ~(relImprove >= improveMin || az1 < az0)
    Icorr2 = Icorr; eta = 0; az1 = az0; rg1 = rg0; phi_lin_1 = zeros(size(Icorr));
end

% report
fprintf('tilt ζ=%.3f | p90(|∂φ/∂Az_rot|): %.4f→%.4f | med|∂φ/∂R|: %.4f→%.4f\n',...
        eta, az0, az1, rg0, rg1);

info.phi_lin  = phi_lin_1;
info.eta      = eta;
info.slope    = slope;
info.azBefore = az0; info.azAfter = az1;
info.rgBefore = rg0; info.rgAfter = rg1;
info.opts     = opts;
end

% ---- helpers ----
% function v = getf(s, f, d)
% if ~isstruct(s) || ~isfield(s,f) || isempty(s.(f)), v=d; else, v=s.(f); end
% end

function a = nanconv1D(a, sigma)
if sigma<=0, return; end
L = max(3, ceil(5*sigma)); x = (-L:L);
g = exp(-(x.^2)/(2*sigma^2)); g = g/sum(g);
m = isfinite(a); a2=a; a2(~m)=0; m2=double(m);
num = conv(a2, g, 'same'); den = conv(m2, g, 'same');
a = num ./ max(den, eps); a(den==0)=NaN;
ix1 = find(a2,1,'first'); ix2 = find(a2,1,'last');
a([1:ix1-1,ix2+1:numel(a)]) = NaN;
end

function [azP90, rgMed] = metrics_rot(I, coh, thetaDeg, cohThresh)
if nargin < 5, preRotated = true; end
Ir = I;  cohr = coh;
if ~preRotated && abs(thetaDeg) > 1e-6
    % If your imrotate lacks FillValues/loose, see mask trick below
    Ir   = imrotate(I,   -thetaDeg, 'bilinear', 'crop');
    cohr = imrotate(coh, -thetaDeg, 'bilinear', 'crop');
end
Ur   = Ir ./ max(abs(Ir), eps);
gAz  = angle(Ur(:,2:end).*conj(Ur(:,1:end-1)));
mAz  = isfinite(gAz) & min(cohr(:,1:end-1),cohr(:,2:end)) >= cohThresh;
azP90= prctile(abs(gAz(mAz)), 90);

U    = I ./ max(abs(I), eps);
gR   = angle(U(2:end,:).*conj(U(1:end-1,:)));
mR   = isfinite(gR) & min(coh(1:end-1,:),coh(2:end,:)) >= cohThresh;
rgMed= median(abs(gR(mR)), 'omitnan');
end

function J = cost_balanced(I, coh, phi, thetaDeg, cohThresh, lambdaRG)
% J = p90(|∂φ/∂Az_rot| after) + λ * max(0, med|∂φ/∂R|_after − med|∂φ/∂R|_before)
[az0, rg0] = metrics_rot(I, coh, thetaDeg, cohThresh);
I2 = I .* exp(-1i*phi);
[az1, rg1] = metrics_rot(I2, coh, thetaDeg, cohThresh);
J = az1 + lambdaRG * max(0, rg1 - rg0);
end

function [Icorr2, touch] = touchup_range_tilt(Icorr, coh, opts)
% tiny per-column range tilt removal on the residual (after main ramp)
if ~isfield(opts,'cohThresh'),  opts.cohThresh  = 0.4;  end
if ~isfield(opts,'sigmaAzTilt'),opts.sigmaAzTilt= 50;   end  % smooth slopes across az

[M,N] = size(Icorr);
U  = Icorr ./ max(abs(Icorr), eps);                 % unit phasors
gR = angle( U(2:end,:).*conj(U(1:end-1,:)) );       % (M-1)×N residual range-grad
W  = min(coh(1:end-1,:), coh(2:end,:));
V  = isfinite(gR) & (W >= opts.cohThresh);

% median range slope per column (robust)
slope = nan(1,N);
for j = 1:N
    v = V(:,j);
    if nnz(v) >= 30
        slope(j) = median(gR(v,j),'omitnan');
    end
end

% smooth the slopes along azimuth
L = max(100, ceil(5.*opts.sigmaAzTilt)); x = (-L:L);
g = exp(-(x.^2)/(2*opts.sigmaAzTilt^2)); g = g/sum(g);
m = isfinite(slope); s2 = slope; s2(~m) = 0; m2 = double(m);
slope_s = conv(s2, g, 'same') ./ max(conv(m2, g, 'same'), eps);
slope_s(~m) = NaN;

% build a column-wise linear phase φ_lin(r,j) = slope_s(j) * r (r starts at 0)
phi_lin = zeros(M,N);
for j = 1:N
    if isfinite(slope_s(j))
        phi_lin(:,j) = (0:M-1).' * slope_s(j);
    else
        phi_lin(:,j) = NaN;
    end
end
% de-mean per column over confident pixels so we don’t inject bias
for j = 1:N
    vv = isfinite(phi_lin(:,j)) & (coh(:,j) >= opts.cohThresh);
    if any(vv), phi_lin(:,j) = phi_lin(:,j) - mean(phi_lin(vv,j)); end
end

% apply
Icorr2 = Icorr .* exp(-1i * phi_lin);

% small report (same mask)
% gR0 = gR;  % before
% U2  = Icorr2 ./ max(abs(Icorr2),eps);
% gR1 = angle( U2(2:end,:).*conj(U2(1:end-1,:)) );
% medR0 = median(abs(gR0(V)),'omitnan');
% medR1 = median(abs(gR1(V)),'omitnan');
% fprintf('touch-up range tilt: med|∂φ/∂R| %.4f → %.4f rad/px\n', medR0, medR1);

touch.slope = slope_s;
touch.phi_lin = phi_lin;
end


% function [Icorr, model] = remove_phzramp_azpca_rot(Ifilt, coh, thetaDeg, opts)
% % Rotate → azPCA → rotate-back ramp → subtract (in original grid).
% % thetaDeg: positive = CCW rotation of image axes; we rotate by -thetaDeg
% %           so that azimuth becomes (approximately) column-aligned.
% %
% % Requires: remove_phzramp_azpca.m in path.
% 
% if nargin < 4, opts = struct(); end
% [M,N] = size(Ifilt);
% 
% % 1) Rotate interferogram & coherence into azimuth-aligned frame
% % [If_rot, maskI]  = rotLike(Ifilt, -thetaDeg);       % complex-safe, NaN outside support
% % [coh_rot, maskC] = rotLike(coh,   -thetaDeg);       % real, NaN outside support
% [If_rot]  = rotfill(Ifilt, -thetaDeg);       % complex-safe, NaN outside support
% [coh_rot] = rotfill(coh,   -thetaDeg);       % real, NaN outside support
% % maskRot = maskI & maskC;
% maskRot = isfinite(If_rot) & isfinite(coh_rot);
% If_rot(~maskRot)  = NaN;
% coh_rot(~maskRot) = NaN;
% 
% % 2) Run azPCA in rotated frame (your fixed version)
% [Icorr_rot, modelRot] = insar.remove_phzramp_azpca(If_rot, coh_rot, opts);
% 
% % 3) Rotate predicted ramp back, apply on original grid
% % [phi_back, mask_back] = rotLike(modelRot.phi_pred, +thetaDeg);
% [phi_back] = rotfill(modelRot.phi_pred, +thetaDeg);
% maskback = 1;
% phi_back(~mask_back) = 0;                           % no correction outside support
% 
% Icorr = Ifilt .* exp(-1i * phi_back);
% 
% % 4) Quick metric (same as before, on original grid)
% gradA  = angle( Ifilt(:,2:end).*conj(Ifilt(:,1:end-1)) );
% gradAc = angle( Icorr(:,2:end).*conj(Icorr(:,1:end-1)) );
% mAz    = isfinite(gradA) & min(coh(:,1:end-1),coh(:,2:end)) >= getOpt(opts,'cohThresh',0.3);
% medAz0 = median(abs(gradA(mAz)),'omitnan');
% medAz1 = median(abs(gradAc(mAz)),'omitnan');
% fprintf('azPCA-ROT(θ=%.2f°): med|∂φ/∂Az| %.4f → %.4f rad/px (Keff=%d)\n', ...
%         thetaDeg, medAz0, medAz1, size(modelRot.modes,2));
% 
% % Pack model
% model = modelRot;
% model.thetaDeg = thetaDeg;
% model.phi_pred_back = phi_back;
% end
% 
% % ---------- helpers ----------
% function [B, maskOut] = rotLike(A, thetaDeg)
% % Robust rotation with NaN support; works on real or complex.
% % Uses 'bilinear' for data, 'nearest' for mask, 'crop' bbox.
%     if isreal(A)
%         mask = isfinite(A);
%         A2 = A; A2(~mask) = 0;
%         B  = imrotate(A2, thetaDeg, 'bilinear', 'crop');
%     else
%         mask = isfinite(real(A)) & isfinite(imag(A));
%         Ar = real(A); Ai = imag(A);
%         Ar(~mask) = 0; Ai(~mask) = 0;
%         Br = imrotate(Ar, thetaDeg, 'bilinear', 'crop');
%         Bi = imrotate(Ai, thetaDeg, 'bilinear', 'crop');
%         B  = complex(Br, Bi);
%     end
%     maskOut = imrotate(mask, thetaDeg, 'nearest', 'crop') > 0.5;
%     B(~maskOut) = NaN;
% end
% 
% function Arot = rotfill(A, angDeg, method, bbox, fillVal)
% % Rotate array A by angDeg (CCW), keep size ('crop'), and fill newly
% % exposed pixels with fillVal (NaN). Works for real or complex.
% if nargin<3||isempty(method), method='bilinear'; end
% if nargin<4||isempty(bbox),   bbox='crop';      end
% if nargin<5, fillVal = NaN; end
% 
% % rotate data
% if isreal(A)
%     Arot = imrotate(A, angDeg, method, bbox);
% else
%     Arot = imrotate(real(A), angDeg, method, bbox) ...
%          + 1i*imrotate(imag(A), angDeg, method, bbox);
% end
% 
% % build a mask of newly created pixels and set them to fillVal
% M = imrotate(ones(size(A),'single'), angDeg, method, bbox);
% edge = (M < 0.99);     % robust threshold
% if isreal(Arot)
%     Arot(edge) = fillVal;
% else
%     Arot(edge) = complex(fillVal, fillVal);
% end
% end
% 
% function v = getOpt(s, name, default)
%     if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name)), v = default; else, v = s.(name); end
% end

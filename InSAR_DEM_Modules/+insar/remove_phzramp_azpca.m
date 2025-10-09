function [Icorr, model] = remove_phzramp_azpca(Ifilt, coh, opts)
% Remove low-rank azimuth undulations using coherence-weighted PCA/SVD.
% Ifilt : MxN complex interferogram (filtered)
% coh   : MxN coherence
% opts (all optional):
%   K=3, cohThresh=0.3, lambda=1e-2, sigmaR=60, sigmaAzPrep=15, maxRows=Inf, protectMask=[]
%   rowCoverageFrac=0.10, rowMinCols=max(10,K+3), rowMinSamples=max(10,2*K)
%
% Output:
%   Icorr         : complex interferogram after ramp removal
%   model.modes   : N x Keff azimuth modes used
%   model.alpha   : M x Keff per-row coefficients (smoothed)
%   model.phi_pred: M x N predicted ramp (zero where not fitted)
%   model.opts    : options used

if nargin < 3, opts = struct(); end
K               = getopt(opts,'K',3);
cohThresh       = getopt(opts,'cohThresh',0.25);
lambda          = getopt(opts,'lambda',1e-2);
sigmaR          = getopt(opts,'sigmaR',25);
sigmaAzPrep     = getopt(opts,'sigmaAzPrep',15);
maxRows         = getopt(opts,'maxRows',Inf);
protectMask     = getopt(opts,'protectMask',[]);
rowCoverageFrac = getopt(opts,'rowCoverageFrac',0.10);
rowMinCols      = getopt(opts,'rowMinCols', max(10, K+3));
rowMinSamples   = getopt(opts,'rowMinSamples', max(10, 2*K));
modeAzSigma     = getopt(opts,'modeAzSigma',25);

[M,N] = size(Ifilt);
phi = angle(Ifilt);

% valid mask (+ optional “do not fit” protection)
valid = isfinite(phi) & isfinite(coh) & (coh >= cohThresh);
if ~isempty(protectMask), valid = valid & ~protectMask; end

% quick coverage debug
nPerRow = sum(valid,2);
% fprintf('azPCA: valid cols/row median=%g (min=%d, max=%d) of N=%d, cohThresh=%.2f\n', ...
%         median(nPerRow), min(nPerRow), max(nPerRow), N, cohThresh);

if nnz(valid) == 0
    Icorr = Ifilt;
    model = struct('modes',[],'alpha',[],'phi_pred',zeros(M,N),'opts',opts);
    return
end

% light azimuth pre-smoothing (NaN-aware)
phi_s = phi;
if sigmaAzPrep > 0
    L = max(3, ceil(5*sigmaAzPrep)); x = (-L:L);
    gk = exp(-(x.^2)/(2*sigmaAzPrep^2)); gk = gk/sum(gk);
    for i = 1:M
        r = phi(i,:); m = valid(i,:);
        r2 = r; r2(~m) = 0; m2 = double(m);
        num = conv(r2, gk, 'same'); den = conv(m2, gk, 'same');
        rsm = num ./ max(den, eps); rsm(den==0) = NaN;
        phi_s(i,:) = rsm;
    end
end

% weighted, row-demeaned matrix A for SVD
W = sqrt(max(coh,0)); W(~valid) = 0;
rm = nan(M,1);
for i=1:M
    m = valid(i,:);
    if any(m), rm(i) = mean(phi_s(i,m)); end
end
rm(~isfinite(rm)) = 0;
A0 = bsxfun(@minus, phi_s, rm); % row demean
A0(~valid) = 0;                 % zero out invalid samples before weighting
A  = W .* A0;

% choose rows for SVD (relaxed)
colsPerRow = sum(valid,2);
rows = find(colsPerRow >= max(rowMinCols, ceil(rowCoverageFrac*N)));
if isempty(rows)
    rows = find(sum(isfinite(phi_s),2) >= rowMinCols);
end
A_svd = A(rows,:);

% % get azimuth modes (fallback to low-freq DCT)
% if ~isempty(A_svd) && any(A_svd(:) ~= 0)
%     [~,S,V] = svd(A_svd,'econ');          % V: N×r
%     r = sum(diag(S) > 0);
%     Keff = min([K, size(V,2), max(1,r)]);
%     modes = V(:,1:Keff);
%     modes = modes - mean(modes,1,'omitnan');    % zero-mean per mode
% else
%     Keff  = max(1, min(K, floor(N/4)));
%     modes = dct_modes(N, Keff);
% end
[~,S,V] = svd(A_svd,'econ');          % singular values s = diag(S)
s = diag(S);
if isempty(s), modes = dct_modes(N, max(1,min(K, floor(N/4)))); 
else
    e = cumsum(s.^2)/max(sum(s.^2),eps);          % cumulative energy
    Keff = min([K, find(e >= 0.95, 1, 'first')]); % ~95% energy
    Keff = max(1, Keff);
    modes = V(:,1:Keff);
    modes = modes - mean(modes,1,'omitnan');      % zero-mean per mode
end

if exist("modeAzSigma",'var') && modeAzSigma > 0
    sig = modeAzSigma;                % e.g., 20–40 px
    L = max(3, ceil(5*sig)); x = -L:L;
    g = exp(-(x.^2)/(2*sig^2)); g = g/sum(g);
    for k = 1:size(modes,2)
        mk = modes(:,k).';                     % row vector (1×N)
        modes(:,k) = conv(mk, g, 'same').';    % smooth along azimuth
    end
end

% per-row weighted ridge fit to modes
alpha = nan(M, Keff);
for i = 1:M
    m = valid(i,:);
    if nnz(m) < rowMinSamples, continue; end
    y  = phi_s(i,m).';
    X  = modes(m,1:Keff);
    w  = W(i,m).';
    Xw = X .* w;
    H  = Xw.'*Xw + lambda*eye(Keff);
    rhs= Xw.'*(w.*y);
    a  = H \ rhs;
    if all(isfinite(a)), alpha(i,1:Keff) = a(:).'; end
end

% smooth coefficients along range
alpha_s = alpha;
if sigmaR > 0
    L = max(3, ceil(5*sigmaR)); xr = (-L:L);
    gr = exp(-(xr.^2)/(2*sigmaR^2)); gr = gr/sum(gr);
    for k = 1:Keff
        ak = alpha(:,k); mk = isfinite(ak);
        a2 = ak; a2(~mk) = 0; m2 = double(mk);
        num = conv(a2, gr, 'same'); den = conv(m2, gr, 'same');
        akS = num ./ max(den, eps); akS(den==0) = NaN;
        alpha_s(:,k) = akS;
    end
end

% reconstruct predicted ramp (start as zeros so no-NaN output)
phi_pred = zeros(M,N,'like',phi);
for i = 1:M
    if any(isfinite(alpha_s(i,:)))
        phi_pred(i,:) = modes * alpha_s(i,:).';
    end
end
% demean per row over valid pixels to avoid DC injection
for i = 1:M
    m = valid(i,:);
    if any(m), phi_pred(i,m) = phi_pred(i,m) - mean(phi_pred(i,m)); end
end

% phi_pred = imgaussfilt(phi_pred,sigmaR./4);
% apply in complex domain
Icorr = Ifilt .* exp(-1i * phi_pred);

% quick metric
% gradA  = angle( Ifilt(:,2:end).*conj(Ifilt(:,1:end-1)) );
% gradAc = angle( Icorr(:,2:end).*conj(Icorr(:,1:end-1)) );
% mAz    = isfinite(gradA) & min(coh(:,1:end-1),coh(:,2:end))>=cohThresh;
% medAz0 = median(abs(gradA(mAz)),'omitnan');
% medAz1 = median(abs(gradAc(mAz)),'omitnan');
% fprintf('azPCA: med|∂φ/∂Az| %.4f → %.4f rad/px (Keff=%d)\n', medAz0, medAz1, size(modes,2));

% pack model
model.modes    = modes;
model.alpha    = alpha_s;
model.phi_pred = phi_pred;
model.opts     = opts;
end

% ---------- helpers ----------
function v = getopt(s, name, default)
if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name))
    v = default;
else
    v = s.(name);
end
end

function B = dct_modes(N,K)
% Low-frequency cosine basis (zero-mean columns)
t = (0:N-1)'; B = zeros(N,K);
for k = 1:K
    B(:,k) = cos(pi*(k-1)*(t + 0.5)/N);
end
B = B - mean(B,1,'omitnan');
end

% function [Icorr, model] = remove_phzramp_azpca(Ifilt, coh, opts)
% % Coherence-weighted PCA (SVD) removal of low-rank azimuth undulations.
% % Ifilt : MxN complex interferogram
% % coh   : MxN coherence
% % opts fields (all optional):
% %   K=3, cohThresh=0.3, lambda=1e-2, sigmaR=60, sigmaAzPrep=15,
% %   maxRows=Inf, protectMask=[]
% 
% if nargin < 3, opts = struct(); end
% K           = getopt(opts,'K',3);
% cohThresh   = getopt(opts,'cohThresh',0.3);
% lambda      = getopt(opts,'lambda',1e-2);
% sigmaR      = getopt(opts,'sigmaR',60);
% sigmaAzPrep = getopt(opts,'sigmaAzPrep',15);
% maxRows     = getopt(opts,'maxRows',Inf);
% protectMask = getopt(opts,'protectMask',[]);
% rowCoverageFrac = getopt(opts,'rowCoverageFrac',0.15);   % was 0.5; much looser
% rowMinCols      = getopt(opts,'rowMinCols',  max(20, K+5));
% 
% 
% [M,N] = size(Ifilt);
% phi = angle(Ifilt);
% 
% % valid mask (+ optional “do not fit” protection)
% valid = isfinite(phi) & isfinite(coh) & (coh >= cohThresh);
% if ~isempty(protectMask), valid = valid & ~protectMask; end
% 
% if nnz(valid) == 0
%     Icorr = Ifilt;
%     model = struct('modes',[],'alpha',[],'phi_pred',zeros(M,N),'opts',opts);
%     return
% end
% 
% % light azimuth pre-smoothing to calm speckle (NaN-aware conv)
% phi_s = phi;
% if sigmaAzPrep > 0
%     L = max(3, ceil(5*sigmaAzPrep)); x = (-L:L);
%     gk = exp(-(x.^2)/(2*sigmaAzPrep^2)); gk = gk/sum(gk);
%     for i = 1:M
%         r = phi(i,:); m = valid(i,:);
%         r2 = r; r2(~m) = 0; m2 = double(m);
%         num = conv(r2, gk, 'same'); den = conv(m2, gk, 'same');
%         rsm = num ./ max(den, eps); rsm(den==0) = NaN;
%         phi_s(i,:) = rsm;
%     end
% end
% 
% % build weighted matrix for SVD (demean per row over valid columns)
% W = sqrt(max(coh,0)); W(~valid) = 0;
% rm = nan(M,1);
% for i=1:M
%     m = valid(i,:);
%     if any(m), rm(i) = mean(phi_s(i,m)); end
% end
% A = W .* (phi_s - rm);                     % NaNs auto-zeroed by W
% 
% % % choose rows for SVD
% % rows = find(sum(valid,2) > 0.5*N);
% % if numel(rows) > maxRows
% %     rows = rows(round(linspace(1,numel(rows),maxRows)));
% % end
% % A_svd = A(rows,:);
% % 
% % % SVD → azimuth modes
% % [~,~,V] = svd(A_svd,'econ');               % V: N×N
% % K = min([K, size(V,2)]);                    % guard
% % modes = V(:,1:K);
% % modes = modes - mean(modes,1,'omitnan');    % zero-mean per mode
% % --- choose rows for SVD (relaxed) ---
% colsPerRow = sum(valid,2);
% rows = find(colsPerRow >= max(rowMinCols, ceil(rowCoverageFrac*N)));
% 
% % if still empty, ignore coherence for row selection as a fallback
% if isempty(rows)
%     rows = find(sum(isfinite(phi_s),2) >= rowMinCols);
% end
% 
% A_svd = A(rows,:);
% 
% % --- get azimuth modes (with fallback if SVD has no rank) ---
% modes = [];
% if ~isempty(A_svd) && any(A_svd(:) ~= 0)
%     [~,S,V] = svd(A_svd, 'econ');          % V: N×r
%     r = sum(diag(S) > 0);
%     Keff = min([K, size(V,2), max(1,r)]);
%     modes = V(:,1:Keff);
%     modes = modes - mean(modes,1,'omitnan');
% else
%     % Fallback: low-frequency DCT modes if SVD yields nothing
%     Keff = max(1, min(K, floor(N/4)));
%     modes = dct_modes(N, Keff);            % helper below
% end
% Keff = size(modes,2);                      % effective K here
% 
% 
% % % per-row weighted ridge fit to modes
% % alpha = nan(M,K);
% % for i = 1:M
% %     m = valid(i,:);
% %     if nnz(m) < K+5, continue; end
% %     y  = phi_s(i,m).';
% %     X  = modes(m,:);
% %     w  = W(i,m).';
% %     Xw = X .* w;
% %     H  = Xw.'*Xw + lambda*eye(K);
% %     rhs= Xw.'*(w.*y);
% %     a  = H \ rhs;
% %     if all(isfinite(a)), alpha(i,:) = a(:).'; end
% % end
% 
% alpha = nan(M,Keff);
% for i = 1:M
%     m = valid(i,:);
%     % require enough columns for this row; relax if needed
%     rowMinSamples = getopt(opts,'rowMinSamples', max(30, 3*Keff));
%     if nnz(m) < rowMinSamples, continue; end
% 
%     y  = phi_s(i,m).';
%     X  = modes(m,1:Keff);
%     w  = W(i,m).';
%     Xw = X .* w;
%     H  = Xw.'*Xw + lambda*eye(Keff);
%     rhs= Xw.'*(w.*y);
%     a  = H \ rhs;
%     if all(isfinite(a)), alpha(i,1:Keff) = a(:).'; end
% end
% 
% 
% % % smooth coefficients along range (rows)
% % alpha_s = alpha;
% % if sigmaR > 0
% %     L = max(3, ceil(5*sigmaR)); xr = (-L:L);
% %     gr = exp(-(xr.^2)/(2*sigmaR^2)); gr = gr/sum(gr);
% %     for k = 1:K
% %         ak = alpha(:,k); mk = isfinite(ak);
% %         a2 = ak; a2(~mk) = 0; m2 = double(mk);
% %         num = conv(a2, gr, 'same'); den = conv(m2, gr, 'same');
% %         akS = num ./ max(den, eps); akS(den==0) = NaN;
% %         alpha_s(:,k) = akS;
% %     end
% % end
% % smooth coefficients along range (rows)
% alpha_s = alpha;
% if sigmaR > 0
%     L = max(3, ceil(5*sigmaR)); xr = (-L:L);
%     gr = exp(-(xr.^2)/(2*sigmaR^2)); gr = gr/sum(gr);
%     for k = 1:Keff                  % <<< was 1:K
%         ak = alpha(:,k); mk = isfinite(ak);
%         a2 = ak; a2(~mk) = 0; m2 = double(mk);
%         num = conv(a2, gr, 'same'); den = conv(m2, gr, 'same');
%         akS = num ./ max(den, eps); akS(den==0) = NaN;
%         alpha_s(:,k) = akS;
%     end
% end
% 
% 
% % reconstruct predicted ramp (demean per row over valid)
% phi_pred = zeros(M,N,'like',phi);
% for i = 1:M
%     if any(isfinite(alpha_s(i,:)))
%         phi_pred(i,:) = modes * alpha_s(i,:).';
%     else
%         phi_pred(i,:) = NaN;
%     end
% end
% for i = 1:M
%     m = valid(i,:);
%     if any(m), phi_pred(i,m) = phi_pred(i,m) - mean(phi_pred(i,m)); end
% end
% 
% % apply in complex domain
% Icorr = Ifilt .* exp(-1i * phi_pred);
% 
% % quick metric
% gradA  = angle( Ifilt(:,2:end).*conj(Ifilt(:,1:end-1)) );
% gradAc = angle( Icorr(:,2:end).*conj(Icorr(:,1:end-1)) );
% mAz    = isfinite(gradA) & min(coh(:,1:end-1),coh(:,2:end))>=cohThresh;
% medAz0 = median(abs(gradA(mAz)),'omitnan');
% medAz1 = median(abs(gradAc(mAz)),'omitnan');
% fprintf('azPCA: med|∂φ/∂Az| %.4f → %.4f rad/px (Keff=%d)\n', ...
%         medAz0, medAz1, Keff);
% model.modes = modes;
% model.alpha = alpha_s;
% model.phi_pred = phi_pred;
% model.opts = opts;
% end
% 
% function v = getopt(s, name, default)
% % safe: never touches a missing field
% if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name))
%     v = default;
% else
%     v = s.(name);
% end
% end
% 
% function B = dct_modes(N,K)
%     t = (0:N-1)'; B = zeros(N,K);
%     for k = 1:K
%         % k=1 is DC; we remove mean later, so these are low-freq cosines
%         B(:,k) = cos(pi*(k-1)*(t + 0.5)/N);
%     end
%     B = B - mean(B,1);  % zero-mean per mode
% end

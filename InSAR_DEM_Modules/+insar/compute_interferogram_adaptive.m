function [phz, cor, ccor] = compute_interferogram_adaptive( ...
    slc1, slc2, qualityThresh, winSize, sigma, alpha)

if nargin < 3 || isempty(qualityThresh), qualityThresh = 0.75; end
if nargin < 4 || isempty(winSize),       winSize       = 5;    end
if nargin < 5 || isempty(sigma),         sigma         = winSize/2; end
if nargin < 6 || isempty(alpha),         alpha         = 2;    end

% Base coherence to build weights (same kernel as your coherence)
[gamma, cgamma, meta] = insar.compute_coherence(slc1, slc2, 'gaussian', winSize, sigma);
[gammar, cgammar, metar] = insar.compute_coherence_robust(slc1, slc2, 'gaussian', winSize, sigma);

% Assume: ccoh_base, ccoh_robust are complex coherences (same window)
gamma0 = abs(cgamma);
gammaR = abs(cgammar);
phi0   = angle(cgamma);
phiR   = angle(cgammar);

dphi = angle(exp(1i*(phiR - phi0)));   % (-pi,pi]
mask = gammaR > 0.3;                   % or 0.4–0.5, depending on your scene

% 1) Phase shift stats
mean_dphi   = nanmean(dphi(mask));
p95_dphi    = prctile(double(dphi(mask)),95);
fprintf('mean Δφ=%.4f rad (%.2f°), 95th=%.3f rad (%.1f°)\n', ...
        mean_dphi, mean_dphi*180/pi, p95_dphi, p95_dphi*180/pi);

% 2) Coherence gain stats
dg = gammaR - gamma0;
fprintf('median Δ|γ|=%.3f, p95 Δ|γ|=%.3f\n', median(dg(mask),'omitnan'), prctile(double(dg(mask)),95));


G = meta.kernel;

% Raw interferogram and unit phasors
Iraw = slc1 .* conj(slc2);
U    = Iraw ./ max(abs(Iraw), eps);

% Coherence weights W = gamma^alpha
W = max(gamma,0).^alpha;

% % One complex, coherence-weighted pass
% num   = imfilter(U.*W, G, 'symmetric');   % complex
% den   = imfilter(W,   G, 'symmetric');    % real
% den(den==0) = eps;
% ccor = num ./ den;
% % Circular Averaged Mulit-Looked Phase
% phz = angle(ccor);

% Multi-looked Amplitude Weighted Complex Coherence
P1  = abs(slc1).^2;  P2 = abs(slc2).^2;
num_ml = imfilter(Iraw.*W, G, 'symmetric');
d1     = imfilter(P1  .*W, G, 'symmetric');
d2     = imfilter(P2  .*W, G, 'symmetric');
ccor = num_ml ./ sqrt(d1.*d2 + eps);  % |gamma_ml| ≤ 1 by construction
phz = angle(ccor);
cor = abs(ccor);                 % <-- coherence MATCHED to this phase
phz(cor < qualityThresh) = NaN;
end

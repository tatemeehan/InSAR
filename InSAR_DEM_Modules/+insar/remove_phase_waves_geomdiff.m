function out = remove_phase_waves_geomdiff(phi, coh, X, Y, H, geom, opts)
% Geometry-aware residual removal using BOTH trajectories (no LOS cube).
% Uses linearization of phi_topo = kz * h with respect to small errors in
% Bperp scale, incidence, and range.
%
% Inputs:
%   phi  : interferometric phase [rad], size MxN
%   coh  : coherence [0..1], size MxN
%   X,Y  : coordinates [m], same size
%   H    : DEM [m], same size
%   geom : struct with fields:
%          .Bperp, .slant, .incidence, .slant2, .incidence2, .lambda, .lookMask
%   opts : struct with fields (all optional):
%          .cohThresh (default 0.6)
%          .useRobust (default true)
%          .usePlane  (default true)  % include ax*x + ay*y
%
% Output:
%   out.phi_corr  : corrected phase
%   out.fit       : fitted residual
%   out.beta      : coefficients [a0 ax ay alpha beta gamma]
%   out.deltaB_med: implied ΔB⊥ [m] from alpha term (global quick QA)
%   out.stats     : regression stats (robustfit) if used

if nargin < 7, opts = struct; end
if ~isfield(opts,'cohThresh'), opts.cohThresh = 0.6; end
if ~isfield(opts,'useRobust'), opts.useRobust = true; end
if ~isfield(opts,'usePlane'),  opts.usePlane  = true; end

% Shorthand (master-side maps)
B  = geom.Bperp;
R1 = geom.slant;
th1= deg2rad(geom.incidence);
R2 = geom.slant2;
th2= deg2rad(geom.incidence2);

lambda = geom.lambda;
mask = geom.lookMask & isfinite(phi) & isfinite(H) & ...
       isfinite(B) & isfinite(R1) & isfinite(th1) & ...
       isfinite(R2) & isfinite(th2) & (coh >= opts.cohThresh);

% Core kz (already uses both r1/r2 via Bperp and (R,theta) of master)
kz = (4*pi/lambda) .* (B ./ max(R1 .* sin(th1), 1e-6));

% Linearized sensitivity features (per pixel):
dth = th2 - th1;                 % incidence mismatch
dR  = R2  - R1;                  % range mismatch

% ∂(kz)/∂θ times h  →  (4π/λ) * h * B * ( -cosθ / (R*sin^2θ) )
F_theta = (4*pi/lambda) .* H .* B .* ( -cos(th1) ./ max(R1 .* (sin(th1).^2), 1e-6) ) .* dth;

% ∂(kz)/∂R times h  →  (4π/λ) * h * B * ( -1 / (R^2 * sinθ) )
F_range = (4*pi/lambda) .* H .* B .* ( -1 ./ max((R1.^2) .* sin(th1), 1e-6) ) .* dR;

% B⊥ scale / kz*h term (captures small multiplicative error in B⊥ or overall kz)
F_kzh = kz .* H;

% Design matrix
if opts.usePlane
    P = [ones(nnz(mask),1), X(mask), Y(mask), F_kzh(mask), F_theta(mask), F_range(mask)];
    labels = {'a0','ax','ay','alpha(kz*h)','beta(theta)','gamma(range)'};
else
    P = [ones(nnz(mask),1), F_kzh(mask), F_theta(mask), F_range(mask)];
    labels = {'a0','alpha(kz*h)','beta(theta)','gamma(range)'};
end

z = phi(mask);

% Fit
if opts.useRobust
    [b,stats] = robustfit(P(:,2:end), z);    % intercept handled internally
    b = [b(1); b(2:end)];
else
    b = P \ z; stats = struct();
end

% Build fitted surface everywhere
if opts.usePlane
    fit = b(1) + b(2).*X + b(3).*Y + b(4).*F_kzh + b(5).*F_theta + b(6).*F_range;
else
    fit = b(1) + b(2).*F_kzh + b(3).*F_theta + b(4).*F_range;
end

phi_corr = phi - fit;

% Quick QA: translate alpha to implied ΔB⊥ at a representative (R,θ)
Rmed  = median(R1(mask),'omitnan');
smed  = median(sin(th1(mask)),'omitnan');
alpha = b( find(strcmp(labels,'alpha(kz*h)'),1,'first') );
if ~isempty(alpha)
    % phi ≈ a0 + alpha*(kz*h) → alpha ~ fractional scale on kz
    % ΔB ≈ alpha * (R sinθ) * (λ/(4π)) / h  → use a representative h scale
    Hmed = median(H(mask),'omitnan'); Hmed = max(Hmed, 1); % avoid divide by ~0
    deltaB_med = alpha * (Rmed*smed) * (lambda/(4*pi)) / Hmed;
else
    deltaB_med = NaN;
end

out = struct('phi_corr',phi_corr,'fit',fit,'beta',b,'labels',{labels},...
             'deltaB_med',deltaB_med,'stats',stats,'mask_used',mask);
end

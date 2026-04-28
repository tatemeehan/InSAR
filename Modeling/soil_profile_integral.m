function [Ig, diag] = soil_profile_integral(f_Hz, kx, prof, Lg_default)
%SOIL_PROFILE_INTEGRAL
% Computes depth integral Ig for the soil return in profile mode:
%   Ig = ∫_0^{zmax} w(z) * exp(-2∫alpha dz) * exp(-i 2∫beta dz) dz
% with w(z)=exp(-z/Lg).
%
% prof required fields:
%   .zmax, .dz
%   .type: 'exp' or 'linear'
%   .VWC0, .VWCinf, .z0
%   .SWC, .sigma_bulk
%   .Tc (scalar) or .Tc_profile struct
%   .soilParams, .wparams, .tauBW, .clayParams
% optional:
%   .Lg (override benchmark)
%   .sigma_profile struct (optional, currently default constant if provided)
%   .freeze.enable (logical) + freeze params for wrapper

% --- Depth grid ---
if ~isfield(prof,'zmax'), prof.zmax = 0.5; end
if ~isfield(prof,'dz'),   prof.dz   = 0.01; end
z = (0:prof.dz:prof.zmax).';
dz = prof.dz;

% --- Lg used in weighting ---
Lg = Lg_default;
if isfield(prof,'Lg') && ~isempty(prof.Lg)
    Lg = prof.Lg;
end
if ~(isfinite(Lg) && Lg > 0)
    error('Profile mode requires finite positive Lg.');
end

% --- T(z) ---
Tc_z = build_T_profile(z, prof);

% --- VWC(z) ---
VWC_z = build_VWC_profile(z, prof);

% Clamp VWC to [0, SWC]
SWC = prof.SWC;
VWC_z = max(min(VWC_z, SWC), 0);

% --- Conductivity profile (optional hook) ---
sigma_bulk = prof.sigma_bulk; % S/m
if isfield(prof,'sigma_profile') && ~isempty(prof.sigma_profile)
    sigma_z = build_sigma_profile(z, prof.sigma_profile, sigma_bulk, VWC_z, Tc_z);
else
    sigma_z = sigma_bulk * ones(size(z));
end

% --- eps(z) ---
use_freeze = isfield(prof,'freeze') && isfield(prof.freeze,'enable') && prof.freeze.enable;
sigma_used_z = nan(size(z));
if ~use_freeze
    % Vectorized evaluation (unfrozen)
    eps_z = soilPermittivity_advanced( ...
        VWC_z, SWC, sigma_z, prof.soilParams, Tc_z, f_Hz, ...
        prof.wparams, prof.tauBW, prof.clayParams);
    sigma_used_z = sigma_z;
    theta_w_liq = VWC_z;
    theta_w_ice = zeros(size(VWC_z));
    fL_z        = ones(size(VWC_z));

else
    % Choose frozen model
    mode = 'physical';
    if isfield(prof.freeze,'mode') && ~isempty(prof.freeze.mode)
        mode = lower(prof.freeze.mode);
    end

    eps_z = complex(zeros(size(z)));
    theta_w_liq = zeros(size(z));
    theta_w_ice = zeros(size(z));
    fL_z        = zeros(size(z));

    for j = 1:numel(z)
        if strcmp(mode,'physical')
            [eps_tmp, comps] = soilPermittivity_frozen_physical( ...
                VWC_z(j), SWC, sigma_z(j), prof.soilParams, Tc_z(j), f_Hz, ...
                prof.wparams, prof.tauBW, prof.clayParams, prof.freeze);
            eps_z(j)        = eps_tmp;
            theta_w_liq(j)  = comps.theta_w_liq;
            theta_w_ice(j)  = comps.theta_w_ice;
            fL_z(j)         = comps.fL;
            sigma_used_z(j) = comps.sigma_eff;
        else
            % Legacy simple wrapper (kept for comparison)
            eps_z(j) = soilPermittivity_frozen_wrapper( ...
                VWC_z(j), SWC, sigma_z(j), prof.soilParams, Tc_z(j), f_Hz, ...
                prof.wparams, prof.tauBW, prof.clayParams, prof.freeze);
            theta_w_liq(j) = VWC_z(j); % unknown split; treat as all liquid
            theta_w_ice(j) = 0;
            fL_z(j)        = 1;
            sigma_used_z(j) = sigma_z(j);
        end
    end
end

% % --- eps(z) --- V1 freeze
% use_freeze = isfield(prof,'freeze') && isfield(prof.freeze,'enable') && prof.freeze.enable;
% 
% if ~use_freeze
%     % Vectorized evaluation
%     eps_z = soilPermittivity_advanced( ...
%         VWC_z, SWC, sigma_z, prof.soilParams, Tc_z, f_Hz, ...
%         prof.wparams, prof.tauBW, prof.clayParams);
% else
%     % Frozen wrapper currently scalar -> loop ok for 30–100 layers
%     eps_z = complex(zeros(size(z)));
%     for j = 1:numel(z)
%         eps_z(j) = soilPermittivity_frozen_wrapper( ...
%             VWC_z(j), SWC, sigma_z(j), prof.soilParams, Tc_z(j), f_Hz, ...
%             prof.wparams, prof.tauBW, prof.clayParams, prof.freeze);
%     end
% end

% --- kz(z) using same convention as kz_from_eps ---
% k0 = 2*pi/lambda = 2*pi*f/c
c = 299792458;
k0 = 2*pi*f_Hz / c;

kz_z = sqrt(k0^2 .* eps_z - kx.^2);
flip = imag(kz_z) < 0;
kz_z(flip) = -kz_z(flip);

beta  = real(kz_z);
alpha = imag(kz_z);

% --- Weighting & cumulative integrals ---
w = exp(-z./Lg);

% A_cum = 2 * cumsum(alpha) * dz;   % 2∫ alpha dz
% B_cum = 2 * cumsum(beta ) * dz;   % 2∫ beta  dz

A_cum = 2 * cumtrapz(z, alpha);
B_cum = 2 * cumtrapz(z, beta);

integrand = w .* exp(-A_cum) .* exp(-1i*B_cum);
Ig = trapz(z, integrand);
% fprintf('Ig = %.6g%+.6gi\n', real(Ig), imag(Ig));
% tmp = 1/Ig;
% fprintf('1/Ig = %.6g%+.6gi\n', real(tmp), imag(tmp));

% --- Effective conductivity (energy-weighted) ---
% Two-way *field* magnitude weighting: |w|^2 * exp(-4∫alpha dz)
A4 = 4 * cumtrapz(z, alpha);          % 4∫alpha dz
W  = (abs(w).^2) .* exp(-A4);         % >=0 weight

den = trapz(z, W);
if den > 0
    sigma_eff = trapz(z, sigma_used_z .* W) / den;
else
    sigma_eff = NaN;
end

% --- Diagnostics ---
diag = struct();
diag.z = z;
diag.VWC = VWC_z;
diag.Tc  = Tc_z;
% diag.sigma_in = sigma_z;
diag.sigma = sigma_used_z;
diag.sigma_eff = sigma_eff;

diag.theta_w_liq = theta_w_liq;
diag.theta_w_ice = theta_w_ice;
diag.fL          = fL_z;

diag.eps = eps_z;
diag.kz  = kz_z;
diag.alpha = alpha;
diag.beta  = beta;

diag.w = w;
diag.A_cum = A_cum;
diag.B_cum = B_cum;
diag.integrand = integrand;

diag.Ig = Ig;
diag.Ig_cum = cumtrapz(z, integrand);

if abs(Ig) > 0
    diag.Dg_eff = 1 / Ig;
else
    diag.Dg_eff = complex(inf, inf);
end

diag.W_energy  = W;        
diag.A4_cum     = A4;       

end

% ------------------ helpers ------------------

function VWC_z = build_VWC_profile(z, prof)
VWC0 = prof.VWC0;
if ~isfield(prof,'VWCinf'), prof.VWCinf = VWC0; end
if ~isfield(prof,'z0'),     prof.z0 = 0.05; end
VWCinf = prof.VWCinf;
z0 = max(prof.z0, eps);

switch lower(prof.type)
    case 'exp'
        VWC_z = VWCinf + (VWC0 - VWCinf).*exp(-z./z0);
    case 'linear'
        VWC_z = VWC0 + (VWCinf - VWC0).*min(z./z0, 1);
    otherwise
        error('Unknown VWC profile type: %s', prof.type);
end
end

function Tc_z = build_T_profile(z, prof)
if isfield(prof,'Tc_profile') && ~isempty(prof.Tc_profile)
    Tp = prof.Tc_profile;
    switch lower(Tp.type)
        case 'linear'
            Tc_z = Tp.Tc0 + Tp.grad*z;
        case 'exp'
            zT = max(Tp.zT, eps);
            Tc_z = Tp.Tc_inf + (Tp.Tc0 - Tp.Tc_inf).*exp(-z./zT);
        otherwise
            error('Unknown Tc_profile type: %s', Tp.type);
    end
else
    Tc_z = prof.Tc * ones(size(z));
end
end

function sigma_z = build_sigma_profile(z, sigma_profile, sigma0, VWC_z, Tc_z)
%BUILD_SIGMA_PROFILE (hook)
% Default behavior: constant sigma with depth.
% Extend later: sigma(z) as function of VWC(z), T(z), salinity profiles, etc.
sigma_z = sigma0 * ones(size(z));

% Example future extension idea (commented):
% if strcmpi(sigma_profile.type,'exp')
%     sigma_inf = sigma_profile.sigma_inf;
%     zsig = sigma_profile.zsig;
%     sigma_z = sigma_inf + (sigma0 - sigma_inf).*exp(-z./max(zsig,eps));
% end
end

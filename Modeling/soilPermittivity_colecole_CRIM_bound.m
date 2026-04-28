function [eps_complex, sigma_eff] = soilPermittivity_colecole_CRIM_bound( ...
    VWC, SWC, sigma_bulk, eps_real_soil, eps_real_clay, frac_clay, Tc, f, wparams, tauBW)
% soilPermittivity_colecole_CRIM_bound
% Complex soil permittivity using:
%   - Cole-Cole water (Kaatze) with VWC-dependent bound-water tau scaling
%   - Debye-type clay dispersion
%   - CRIM sqrt mixing
%   - Optional bulk conductivity term
%
% Inputs:
%   VWC, SWC, sigma_bulk, eps_real_soil, eps_real_clay, frac_clay, Tc, f
%     as before.
%   wparams: struct for water -> Cole-Cole
%       .epssCoefs, .tauCoefs, .eps_inf, .alpha
%   tauBW: struct for bound-water scaling
%       .A, .B, .tau_scale, .max_tau
%
% Output:
%   eps_complex: complex soil permittivity

eps0  = 8.854187817e-12;
omega = 2*pi*f;

% Columnize
VWC        = VWC(:);
SWC        = SWC(:);
Tc         = Tc(:);
sigma_bulk = sigma_bulk(:);

% --- 1. Volume fractions ---
theta_air         = max(SWC - VWC, 0);
theta_solid_total = max(1 - SWC, 0);
theta_clay        = frac_clay .* theta_solid_total;
theta_solid_nonclay = theta_solid_total - theta_clay;

% --- 2. Water permittivity: Kaatze + bound-water tau scaling ---
G        = [ones(numel(Tc),1), Tc, Tc.^2, Tc.^3];
eps_s_T  = G * wparams.epssCoefs(:);
tau_free = G * wparams.tauCoefs(:);   % Kaatze tau(T)

% Bound-water scaling
A         = tauBW.A;
B         = tauBW.B;
tau_scale = tauBW.tau_scale;
max_tau   = tauBW.max_tau;

Veps = max(VWC, 1e-5);
tau_eff = tau_free .* A .* exp(B ./ Veps) .* tau_scale;
tau_eff = min(tau_eff, max_tau);

params_w = struct();
params_w.eps_inf = wparams.eps_inf;
params_w.alpha   = wparams.alpha;
params_w.sigma   = 0;    % we can keep conduction separate

eps_water = zeros(size(VWC)) + 1i*zeros(size(VWC));
for k = 1:numel(VWC)
    params_w.eps_s = eps_s_T(k);
    params_w.tau   = tau_eff(k);
    eps_water(k)   = water_permittivity_colecole(f, params_w);
end

% --- 3. Clay Debye dispersion ---
eps_inf_clay = 5;
tau_clay     = 5e-10;
alpha_clay   = 0.15;

eps_clay = eps_inf_clay + ...
    (eps_real_clay - eps_inf_clay) ./ (1 - (1i*omega*tau_clay).^(1 - alpha_clay));

% --- 4. Solids & air ---
eps_solid = eps_real_soil;
eps_air   = 1.0;

sqrt_eps_air    = sqrt(eps_air);
sqrt_eps_water  = sqrt(eps_water);
sqrt_eps_solid  = sqrt(eps_solid);
sqrt_eps_clay   = sqrt(eps_clay);

% --- 5. CRIM mixing ---
sqrt_eps_mix = ...
    theta_air           .* sqrt_eps_air   + ...
    theta_solid_nonclay .* sqrt_eps_solid + ...
    theta_clay          .* sqrt_eps_clay  + ...
    VWC                 .* sqrt_eps_water;

eps_mix = sqrt_eps_mix.^2;

% --- 6. Bulk conductivity term (optional) ---
% If sigma_bulk is your *input*, you add:
eps_cond = 1i .* sigma_bulk ./ (omega * eps0);

eps_complex = eps_mix + eps_cond;

% eps_im = imag(eps_cond);  % this will be negative for conductive loss
% Total imaginary part (negative for loss under this convention)
eps_im_total = imag(eps_complex);
eps_im_cond = imag(eps_cond);                    % < 0
eps_im_dip  = eps_im_total - eps_im_cond;        % dipolar-only imag

% sigma_eff = omega * eps0 * eps_im;
% If you want sigmas:
sigma_eff_total = omega .* eps0 .* eps_im_total;    % total
sigma_eff_cond  = sigma_bulk;                        % by construction
sigma_eff_dip   = omega .* eps0 .* eps_im_dip;      % effective dipolar sigma

end

function [eps_complex, components] = soilPermittivity_frozen_physical( ...
    VWC, SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams, freeze)
%SOILPERMITTIVITY_FROZEN_PHYSICAL
% Physical frozen extension of soilPermittivity_advanced:
% - splits pore water into liquid + ice using unfrozen water fraction fL(T)
% - CRIM sqrt mixing across air/solid/clay/liquid-water/ice
% - conductivity reduced with liquid fraction connectivity model
%
% Conventions:
%   imag(eps) > 0 is loss. We add conduction as +i*sigma/(omega*eps0).

eps0  = 8.854187817e-12;
omega = 2*pi*f;

% Ensure column
VWC = VWC(:); Tc = Tc(:); sigma_bulk = sigma_bulk(:);
if isscalar(SWC), SWC = SWC * ones(size(VWC)); else, SWC = SWC(:); end

% % Defaults for freeze
% if nargin < 10 || isempty(freeze), freeze = struct(); end
% if ~isfield(freeze,'uwf') || isempty(freeze.uwf), freeze.uwf = struct('mode','logistic'); end
% if ~isfield(freeze,'eta_sigma') || isempty(freeze.eta_sigma), freeze.eta_sigma = 2; end
% if ~isfield(freeze,'sigma_floor') || isempty(freeze.sigma_floor), freeze.sigma_floor = 0.0; end  % S/m

% Defaults for freeze
if nargin < 10 || isempty(freeze), freeze = struct(); end

% SAFETY: if freeze.enable exists and is false, behave like "unfrozen"
% (i.e., no liquid/ice split; all VWC treated as liquid)
if isfield(freeze,'enable') && ~isempty(freeze.enable) && ~freeze.enable
    freeze.uwf = struct('mode','force_liquid');  % custom mode handled below
end

if ~isfield(freeze,'uwf') || isempty(freeze.uwf)
    freeze.uwf = struct('mode','logistic');
end
if ~isfield(freeze,'eta_sigma') || isempty(freeze.eta_sigma), freeze.eta_sigma = 2; end
if ~isfield(freeze,'sigma_floor') || isempty(freeze.sigma_floor), freeze.sigma_floor = 0.0; end


% --- 1) Volume fractions for air/solids/clay (same as your advanced model) ---
theta_air         = max(SWC - VWC, 0);
theta_solid_total = max(1 - SWC, 0);

frac_clay = soilParams.frac_clay;
theta_clay          = frac_clay .* theta_solid_total;
theta_solid_nonclay = theta_solid_total - theta_clay;

% --- 2) Liquid fraction fL(T) and split water into liquid+ice ---
% [fL, uwf_info] = unfrozen_water_fraction(Tc, VWC, SWC, soilParams, freeze.uwf);
if isfield(freeze,'uwf') && isfield(freeze.uwf,'mode') && strcmpi(freeze.uwf.mode,'force_liquid')
    fL = ones(size(VWC));
    uwf_info = struct('mode','force_liquid');
else
    [fL, uwf_info] = unfrozen_water_fraction(Tc, VWC, SWC, soilParams, freeze.uwf);
end


theta_w_liq = fL .* VWC;
theta_w_ice = max(VWC - theta_w_liq, 0);

% --- 3) Permittivities of constituents ---
eps_air   = 1.0;
eps_solid = soilParams.eps_real_soil;

% Clay dispersion (your Cole–Cole clay)
eps_clay  = clay_permittivity_colecole(f, soilParams.eps_real_clay, clayParams);

% Liquid water (your dual Cole–Cole)
eps_water = water_permittivity_dual_colecole_bound(f, Tc, theta_w_liq, wparams, tauBW);

% Ice permittivity: reuse your snow ice function if available
% (If not on path, replace with a simple constant model)
if exist('ice_permittivity_maetzler06','file') == 2
    eps_ice = ice_permittivity_maetzler06(f, Tc + 273.15); % expects Kelvin in your version
else
    % fallback: weakly lossy ice
    warning('No Ice Permittivity Model Included! Using Default: 3.15 + 1i*1e-3')
    eps_ice = (3.15 + 1i*1e-3) * ones(size(Tc));
end

% --- 4) CRIM sqrt mixing ---
sqrt_mix = ...
    theta_air           .* sqrt(eps_air)   + ...
    theta_solid_nonclay .* sqrt(eps_solid) + ...
    theta_clay          .* sqrt(eps_clay)  + ...
    theta_w_liq         .* sqrt(eps_water) + ...
    theta_w_ice         .* sqrt(eps_ice);

eps_mix = sqrt_mix.^2;

% --- 5) Conductivity model (frozen-aware through theta_w_liq) ---
% NOTE: theta_w_liq already includes freezing via fL(T):
%   theta_w_liq = fL .* VWC

if ~isfield(freeze,'sigma_model') || isempty(freeze.sigma_model)
    freeze.sigma_model = struct(); % defaults in sigma_soil_model
end

% sigma_soil_model returns dS/m (HydraGO-style units)
sigma_dSm = sigma_soil_model(theta_w_liq, SWC, soilParams.frac_clay, freeze.sigma_model);

% Convert to S/m for EM conduction term
sigma_eff = 0.1 * sigma_dSm;  % dS/m -> S/m

% Optional safety floor in S/m
if isfield(freeze,'sigma_floor') && ~isempty(freeze.sigma_floor)
    sigma_eff = max(sigma_eff, freeze.sigma_floor);
end

% Add conduction with imag>0 convention
eps_cond = 1i * (sigma_eff ./ (omega * eps0));
eps_complex = eps_mix + eps_cond;


% % --- 5) Conductivity reduction tied to liquid connectivity ---
% % sigma_eff = sigma_floor + sigma_bulk * (theta_w_liq / VWC)^eta
% liq_frac = theta_w_liq ./ max(VWC, realmin);
% sigma_eff = freeze.sigma_floor + sigma_bulk .* (liq_frac .^ freeze.eta_sigma);
% 
% % Add conduction with imag>0 convention
% eps_cond = 1i * (sigma_eff ./ (omega * eps0));
% 
% eps_complex = eps_mix + eps_cond;

% Diagnostics
if nargout > 1
    components = struct();
    components.eps_mix = eps_mix;
    components.eps_cond = eps_cond;
    components.eps_water = eps_water;
    components.eps_ice = eps_ice;
    components.eps_clay = eps_clay;

    components.theta_air = theta_air;
    components.theta_solid = theta_solid_total;
    components.theta_clay = theta_clay;
    components.theta_solid_nonclay = theta_solid_nonclay;
    components.theta_w_liq = theta_w_liq;
    components.theta_w_ice = theta_w_ice;

    components.fL = fL;
    components.uwf_info = uwf_info;

    components.sigma_eff = sigma_eff;
end
end

function epsfz = soilPermittivity_frozen_wrapper(VWC, SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams, freeze)
%SOILPERMITTIVITY_FROZEN_WRAPPER
% Minimal frozen-water partition wrapper around soilPermittivity_advanced.
% Tc in degC. Produces eps with +imag loss convention.

% Defaults
if ~isfield(freeze,'T0'), freeze.T0 = 0; end
if ~isfield(freeze,'dT'), freeze.dT = 0.5; end
if ~isfield(freeze,'p_sigma'), freeze.p_sigma = 2; end

% Smooth liquid fraction (0..1)
dT = max(freeze.dT, eps);  % avoid divide-by-zero robustly
fL = 1 ./ (1 + exp(-(Tc - freeze.T0)./dT));


VWC_liq = fL .* VWC;
VWC_ice = (1 - fL) .* VWC;

% Suppress conductivity when frozen
sigma_eff = sigma_bulk .* (fL.^freeze.p_sigma);

% Liquid-only soil
eps_liq = soilPermittivity_advanced( ...
    VWC_liq, SWC, sigma_eff, soilParams, Tc, f, wparams, tauBW, clayParams);

% Ice permittivity (your snow ice model; Tc->K)
eps_ice = ice_permittivity_maetzler06(f, Tc + 273.15);

% Crude pore-water mixing between liquid and ice in sqrt domain
if VWC > 0
    w_liq = VWC_liq / VWC;
    w_ice = VWC_ice / VWC;
else
    w_liq = 0; w_ice = 0;
end

epsfz = (w_liq.*sqrt(eps_liq) + w_ice.*sqrt(eps_ice)).^2;

% Enforce +imag loss convention
epsfz = real(epsfz) + 1i*abs(imag(epsfz));
end

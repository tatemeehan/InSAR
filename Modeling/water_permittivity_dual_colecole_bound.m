function eps_water = water_permittivity_dual_colecole_bound(f, Tc, VWC, wparams, tauBW)
% WATER_PERMITTIVITY_DUAL_COLECOLE_BOUND
% Dual Cole–Cole water permittivity with a free-water mode and a
% slower bound-water mode, including mild VWC-dependent tau scaling.
%
% Inputs:
%   f       : frequency (Hz), scalar
%   Tc      : temperature (°C), N×1
%   VWC     : volumetric water content (m3/m3), N×1
%   wparams : struct with fields
%       .epssCoefs   : 4×1 poly for eps_s(T)
%       .tauCoefs    : 4×1 poly for tau_free(T)
%       .eps_inf     : high-freq permittivity
%       .alpha_free  : Cole–Cole alpha for free-water mode
%       .alpha_bound : Cole–Cole alpha for bound-water mode
%       .frac_bound  : fraction (0–1) of Δε assigned to bound mode
%   tauBW   : struct with fields
%       .free_scale  : scalar for tau_free (typically ~1)
%       .bound_scale : scalar >1 for bound-water tau (slower)
%       .B_bound     : VWC sensitivity in exp(B/θ)
%       .max_tau     : cap for bound tau
%
% Output:
%   eps_water : N×1 complex water permittivity (effective)

    eps0  = 8.854187817e-12; %#ok<NASGU>  % not used directly here
    omega = 2*pi*f;

    Tc   = Tc(:);
    VWC  = VWC(:);
    N    = numel(Tc);

    % --- Temperature-dependent eps_s(T) and tau_free(T) from Kaatze polys ---
    G        = [ones(N,1), Tc, Tc.^2, Tc.^3];
    eps_s_T  = G * wparams.epssCoefs(:);   % static permittivity
    tau_free = G * wparams.tauCoefs(:);    % relaxation time (s)

    % --- Split Δε into free and bound contributions ---
    eps_inf = wparams.eps_inf;
    delta_eps_total = eps_s_T - eps_inf;

    frac_bound = wparams.frac_bound;     % 0–1
    frac_bound = max(0, min(1, frac_bound));

    delta_eps_bound = frac_bound    .* delta_eps_total;
    delta_eps_free  = (1-frac_bound).* delta_eps_total;

    % --- Effective taus for each mode ---
    alpha_free  = wparams.alpha_free;
    alpha_bound = wparams.alpha_bound;

    % Free-water tau scaling (usually ~1)
    tau_free_eff = tau_free .* tauBW.free_scale;

    % Bound-water tau: slower, VWC-dependent
    Veps = max(VWC, 1e-5);
    tau_bound_eff = tau_free .* tauBW.bound_scale .* exp(tauBW.B_bound ./ Veps);
    tau_bound_eff = min(tau_bound_eff, tauBW.max_tau);

    % --- Cole–Cole terms (vectorized over N) ---
    jomega_tau_free  = (1i * omega .* tau_free_eff ).^(1 - alpha_free);
    jomega_tau_bound = (1i * omega .* tau_bound_eff).^(1 - alpha_bound);

    % term_free  = delta_eps_free  ./ (1 + jomega_tau_free );
    % term_bound = delta_eps_bound ./ (1 + jomega_tau_bound);
    term_free  = delta_eps_free  ./ (1 - (1i * omega * tau_free_eff).^(1 - alpha_free));
term_bound = delta_eps_bound ./ (1 - (1i * omega * tau_bound_eff).^(1 - alpha_bound));


    eps_water = eps_inf + term_free + term_bound;
end

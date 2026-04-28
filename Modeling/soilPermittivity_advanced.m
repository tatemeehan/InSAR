function [eps_complex, components] = soilPermittivity_advanced( ...
    VWC, SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams)
% SOILPERMITTIVITY_ADVANCED
% Complex soil permittivity with:
%   - Dual Cole–Cole water (free + bound) with VWC-dependent tau
%   - Tunable clay Cole–Cole dispersion
%   - CRIM sqrt mixing (air / solid / clay / water)
%   - Explicit bulk conductivity term (σ_bulk)
%
% Inputs:
%   VWC        : volumetric water content (m3/m3), N×1
%   SWC        : porosity / saturation water content, scalar or N×1
%   sigma_bulk : bulk conductivity (S/m), N×1 (e.g. HydraGO EC * 0.1)
%   soilParams : struct with fields
%       .eps_real_soil
%       .eps_real_clay
%       .frac_clay        (fraction of solid volume that is clay)
%   Tc         : soil temperature (°C), N×1
%   f          : frequency (Hz), scalar
%   wparams, tauBW : water model parameter structs (see above)
%   clayParams     : clay dispersion parameters (see above)
%
% Outputs:
%   eps_complex : N×1 complex permittivity of soil
%   components  : struct with intermediate pieces (useful for debugging)

    eps0  = 8.854187817e-12;
    omega = 2*pi*f;

    % Ensure column vectors
    VWC        = VWC(:);
    SWC        = SWC(:);
    Tc         = Tc(:);
    sigma_bulk = sigma_bulk(:);
    N          = numel(VWC);

    % --- 1. Volume fractions ---
    theta_air         = max(SWC - VWC, 0);       % air in pore space
    theta_solid_total = max(1 - SWC, 0);        % total solids
    frac_clay         = soilParams.frac_clay;

    theta_clay          = frac_clay .* theta_solid_total;
    theta_solid_nonclay = theta_solid_total - theta_clay;

    % --- 2. Water permittivity (dual Cole–Cole, VWC-dependent) ---
    eps_water = water_permittivity_dual_colecole_bound( ...
        f, Tc, VWC, wparams, tauBW);

    % --- 3. Clay permittivity (Cole–Cole) ---
    eps_clay = clay_permittivity_colecole( ...
        f, soilParams.eps_real_clay, clayParams);

    % --- 4. Solids and air permittivities ---
    eps_air   = 1.0;
    eps_solid = soilParams.eps_real_soil;

    sqrt_eps_air    = sqrt(eps_air);
    sqrt_eps_solid  = sqrt(eps_solid);
    sqrt_eps_clay   = sqrt(eps_clay);
    sqrt_eps_water  = sqrt(eps_water);

    % --- 5. CRIM mixing ---
    sqrt_eps_mix = ...
        theta_air           .* sqrt_eps_air   + ...
        theta_solid_nonclay .* sqrt_eps_solid + ...
        theta_clay          .* sqrt_eps_clay  + ...
        VWC                 .* sqrt_eps_water;

    eps_mix = sqrt_eps_mix.^2;

    % --- 6. Conductivity term (your sign convention) ---
    eps_cond = sigma_bulk ./ (-1i * omega * eps0);   % imag < 0 for sigma>0

    eps_complex = eps_mix + eps_cond;

    % --- 7. Optional components for diagnostics ---
    if nargout > 1
        eps_im_total = imag(eps_complex);
        eps_im_cond  = imag(eps_cond);
        eps_im_dip   = eps_im_total - eps_im_cond;

        sigma_eff_total = omega .* eps0 .* eps_im_total;  % effective total loss
        sigma_eff_cond  = sigma_bulk;
        sigma_eff_dip   = omega .* eps0 .* eps_im_dip;

        components = struct();
        components.eps_mix         = eps_mix;
        components.eps_water       = eps_water;
        components.eps_clay        = eps_clay;
        components.theta_air       = theta_air;
        components.theta_solid     = theta_solid_total;
        components.theta_clay      = theta_clay;
        components.theta_solid_nonclay = theta_solid_nonclay;

        components.eps_cond        = eps_cond;
        components.eps_im_total    = eps_im_total;
        components.eps_im_cond     = eps_im_cond;
        components.eps_im_dip      = eps_im_dip;

        components.sigma_eff_total = sigma_eff_total;
        components.sigma_eff_cond  = sigma_eff_cond;
        components.sigma_eff_dip   = sigma_eff_dip;
    end
end

function out = soil_eps_from_hydra_forward( ...
    VWC, SWC, EC_dSm, Tc, f, ...
    soilParams, wparams, tauBW)
% SOIL_EPS_FROM_HYDRA_FORWARD
% Convenience wrapper around soilPermittivity_colecole_CRIM_bound to
% go from HydraGO-style inputs (VWC, EC, T) to complex permittivity at f,
% split into dipolar and conductive parts.
%
% Inputs:
%   VWC      : volumetric water content (m3/m3)
%   SWC      : porosity / saturation water content (m3/m3)
%   EC_dSm   : bulk EC in dS/m (e.g. from HydraGO TC-corrected EC)
%   Tc       : soil temperature (deg C)
%   f        : frequency (Hz) [scalar]
%   soilParams (optional struct):
%       .eps_real_soil  (default 4.0)
%       .eps_real_clay  (default 25)
%       .frac_clay      (default 0.15)
%   wparams (struct, required):
%       .epssCoefs, .tauCoefs, .eps_inf, .alpha   (Kaatze-based)
%   tauBW (struct, optional):
%       .A, .B, .tau_scale, .max_tau
%
% Output:
%   out: struct with fields:
%       .eps_complex
%       .eps_real
%       .eps_im_total
%       .eps_im_cond
%       .eps_im_dip
%       .sigma_bulk          (S/m)
%       .sigma_eff_total     (S/m)
%       .sigma_eff_cond      (S/m)
%       .sigma_eff_dip       (S/m)

    % --- Defaults for soil parameters ---
    if nargin < 6 || isempty(soilParams)
        soilParams.eps_real_soil = 4.0;
        soilParams.eps_real_clay = 25.0;
        soilParams.frac_clay     = 0.15;
    else
        if ~isfield(soilParams,'eps_real_soil'), soilParams.eps_real_soil = 4.0; end
        if ~isfield(soilParams,'eps_real_clay'), soilParams.eps_real_clay = 25.0; end
        if ~isfield(soilParams,'frac_clay'),     soilParams.frac_clay     = 0.15; end
    end

    % --- Defaults for bound-water tau scaling ---
    if nargin < 8 || isempty(tauBW)
        tauBW.A         = 1.0;
        tauBW.B         = 0.01;
        tauBW.tau_scale = 1.0;
        tauBW.max_tau   = 5e-10;
    end

    % --- EC (dS/m) -> sigma_bulk (S/m) ---
    sigma_bulk = EC_dSm(:) * 0.1;   % 1 dS/m = 0.1 S/m

    % --- Call your forward model (with corrected conductivity term) ---
    [eps_complex] = soilPermittivity_colecole_CRIM_bound( ...
        VWC, SWC, sigma_bulk, ...
        soilParams.eps_real_soil, ...
        soilParams.eps_real_clay, ...
        soilParams.frac_clay, Tc, f, ...
        wparams, tauBW);

    % --- Split into components ---
    eps0 = 8.854187817e-12;
    omega = 2*pi*f;

    eps_real      = real(eps_complex);
    eps_im_total  = imag(eps_complex);                 % < 0 for loss

    % Conductive-only imag part (from sigma_bulk term)
    eps_cond = sigma_bulk ./ (-1i * omega * eps0);
    eps_im_cond = imag(eps_cond);                      % < 0

    % Dipolar imag = total - conductive
    eps_im_dip = eps_im_total - eps_im_cond;

    % Effective conductivities
    sigma_eff_total = omega .* eps0 .* eps_im_total;  % S/m
    sigma_eff_cond  = sigma_bulk;                      % S/m, by construction
    sigma_eff_dip   = omega .* eps0 .* eps_im_dip;    % S/m-equivalent dipolar loss

    % --- Pack outputs ---
    out = struct();
    out.eps_complex      = eps_complex;
    out.eps_real         = eps_real;
    out.eps_im_total     = eps_im_total;
    out.eps_im_cond      = eps_im_cond;
    out.eps_im_dip       = eps_im_dip;
    out.sigma_bulk       = sigma_bulk;
    out.sigma_eff_total  = sigma_eff_total;
    out.sigma_eff_cond   = sigma_eff_cond;
    out.sigma_eff_dip    = sigma_eff_dip;
end

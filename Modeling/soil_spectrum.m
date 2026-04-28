function results = soil_spectrum( ...
    VWC, SWC, EC_dSm, Tc, freqs, ...
    soilParams, wparams, tauBW, clayParams)
% SOIL_SPECTRUM_FROM_HYDRA
% Run the forward soil model at multiple frequencies and organize
% outputs into a struct array, one element per frequency.
%
% Inputs:
%   VWC      : volumetric water content (m3/m3), N×1
%   SWC      : porosity / saturation water content (m3/m3), scalar or N×1
%   EC_dSm   : bulk EC (dS/m), N×1 (e.g., HydraGO temp-corrected EC)
%   Tc       : soil temperature (°C), N×1
%   freqs    : vector of frequencies [Hz], 1×F
%   soilParams, wparams, tauBW as before
%
% Output:
%   results(F): struct array with fields
%       .freq
%       .eps_complex, .eps_real, .eps_im_total
%       .eps_im_cond, .eps_im_dip
%       .sigma_bulk, .sigma_eff_total
%       .sigma_eff_cond, .sigma_eff_dip

    if nargin < 6 || isempty(soilParams)
        soilParams.eps_real_soil = 4.0;
        soilParams.eps_real_clay = 25.0;
        soilParams.frac_clay     = 0.15;
    else
        if ~isfield(soilParams,'eps_real_soil'), soilParams.eps_real_soil = 4.0; end
        if ~isfield(soilParams,'eps_real_clay'), soilParams.eps_real_clay = 25.0; end
        if ~isfield(soilParams,'frac_clay'),     soilParams.frac_clay     = 0.15; end
    end
    if nargin < 7 || isempty(wparams)
        % Water parameters (Kaatze-polynomial based)
        epssCoefs = [88.119019607843190,-0.484234778121778,0.004464654282766,-2.393188854489220e-05]';
        tauCoefs  = [1.736936274509802e-11,-4.823117475060188e-13,5.488493292053644e-15,-2.366013071895416e-17]';
        wparams.epssCoefs = epssCoefs;   % your 4-coef vec
        wparams.tauCoefs  = tauCoefs;
        wparams.eps_inf   = 4.8;
        wparams.alpha     = 0.001;
    end
    if nargin < 8 || isempty(tauBW)
        tauBW.A         = 1.0;
        tauBW.B         = 0.01;
        tauBW.tau_scale = 1.0;
        tauBW.max_tau   = 5e-10;
        tauBW.free_scale  = 1.0;
        tauBW.bound_scale = 5.0;      % bound water ~5× slower than free
        tauBW.B_bound     = 0.02;     % VWC sensitivity (tune)
    end
    if nargin < 9 || isempty(clayParams)
        clayParams.eps_inf_clay = 5;
        clayParams.tau_clay     = 5e-10;
        clayParams.alpha_clay   = 0.15;
    end

    freqs  = freqs(:).';           % ensure row
    nF     = numel(freqs);
    results(nF) = struct();

    % Precompute sigma_bulk in S/m
    sigma_bulk = EC_dSm(:) * 0.1;  % dS/m -> S/m

    eps0 = 8.854187817e-12;

    for k = 1:nF
        f     = freqs(k);
        omega = 2*pi*f;

        % Call your forward model
        % eps_complex = soilPermittivity_colecole_CRIM_bound( ...
        %     VWC, SWC, sigma_bulk, ...
        %     soilParams.eps_real_soil, ...
        %     soilParams.eps_real_clay, ...
        %     soilParams.frac_clay, Tc, f, ...
        %     wparams, tauBW);

        [eps_complex, components] = soilPermittivity_advanced( ...
    VWC, SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams);

        eps_real     = real(eps_complex);
        eps_im_total = imag(eps_complex);

        % Conductive imag part from sigma_bulk term (using your sign convention)
        eps_cond   = sigma_bulk ./ (-1i * omega * eps0);
        eps_im_cond = imag(eps_cond);
        eps_im_dip  = eps_im_total - eps_im_cond;

        sigma_eff_total = omega .* eps0 .* eps_im_total;  % effective total loss (your convention)
        sigma_eff_cond  = sigma_bulk;                     % S/m (true conduction)
        sigma_eff_dip   = omega .* eps0 .* eps_im_dip;    % dipolar loss in S/m units

        % Store in struct
        results(k).freq            = f;
        results(k).eps_complex     = eps_complex;
        results(k).eps_real        = eps_real;
        results(k).eps_im_total    = eps_im_total;
        results(k).eps_im_cond     = eps_im_cond;
        results(k).eps_im_dip      = eps_im_dip;
        results(k).sigma_bulk      = sigma_bulk;
        results(k).sigma_eff_total = sigma_eff_total;
        results(k).sigma_eff_cond  = sigma_eff_cond;
        results(k).sigma_eff_dip   = sigma_eff_dip;

        % s          = struct();
        % s.freq     = f;
        % s.eps_complex     = eps_complex;
        % s.eps_real        = eps_real;
        % s.eps_im_total    = eps_im_total;
        % s.eps_im_cond     = eps_im_cond;
        % s.eps_im_dip      = eps_im_dip;
        % s.sigma_bulk      = sigma_bulk;
        % s.sigma_eff_total = sigma_eff_total;
        % s.sigma_eff_cond  = sigma_eff_cond;
        % s.sigma_eff_dip   = sigma_eff_dip;
        % 
        % results(k) = s;
    end
end

function plot_soil_dielectric_decomposition( ...
    soilParams, wparams, tauBW, clayParams, ...
    VWC, SWC, EC_dSm, Tc)

% PLOT_SOIL_DIELECTRIC_DECOMPOSITION
% Decompose complex soil permittivity vs frequency into:
%   - total soil response
%   - conduction loss
%   - dipolar (water+clay) loss
%   - material-level water & clay spectra
%
% Inputs:
%   soilParams : struct with fields
%       .eps_real_soil
%       .eps_real_clay
%       .frac_clay
%   wparams    : water permittivity struct (with Kaatze polys, alpha's, frac_bound)
%   tauBW      : bound-water tau scaling struct
%   clayParams : struct with fields
%       .eps_inf_clay
%       .tau_clay
%       .alpha_clay
%   VWC        : volumetric water content (m3/m3), scalar or vector (pick one sample)
%   SWC        : porosity
%   EC_dSm     : bulk EC (dS/m) for that sample
%   Tc         : soil temperature (°C) for that sample
%
% Note:
%   This uses soilPermittivity_advanced(VWC,SWC,sigma_bulk,soilParams,Tc,f,wparams,tauBW,clayParams)
%   returning [eps_complex, components].

    % --- Pick a single representative sample (if arrays) ---
    VWC   = VWC(1);
    SWC   = SWC(1);
    EC_dSm = EC_dSm(1);
    Tc    = Tc(1);

    sigma_bulk = EC_dSm * 0.1;   % dS/m -> S/m

    % --- Frequency grid ---
    f_min = 50e6;    % 50 MHz
    f_max = 5e9;     % 5 GHz
    Nf    = 200;
    freqs = logspace(log10(f_min), log10(f_max), Nf).';

    % Preallocate
    eps_complex_all   = complex(zeros(Nf,1));
    eps_real_all      = zeros(Nf,1);
    eps_im_total_all  = zeros(Nf,1);
    eps_im_cond_all   = zeros(Nf,1);
    eps_im_dip_all    = zeros(Nf,1);

    % Material-level water/clay (not volume-mixed, just their intrinsic spectra)
    eps_water_all     = complex(zeros(Nf,1));
    eps_clay_all      = complex(zeros(Nf,1));

    for k = 1:Nf
        f = freqs(k);

        [eps_complex, comp] = soilPermittivity_advanced( ...
            VWC, SWC, sigma_bulk, soilParams, Tc, f, ...
            wparams, tauBW, clayParams);

        eps_complex_all(k)  = eps_complex;
        eps_real_all(k)     = real(eps_complex);
        eps_im_total_all(k) = comp.eps_im_total;
        eps_im_cond_all(k)  = comp.eps_im_cond;
        eps_im_dip_all(k)   = comp.eps_im_dip;

        % Material-level
        eps_water_all(k) = comp.eps_water;
        eps_clay_all(k)  = comp.eps_clay;
    end

    % --- Plot real part ---
    figure;
    subplot(2,1,1);
    semilogx(freqs, eps_real_all, 'k-', 'LineWidth', 1.8); hold on;
    ylabel('\epsilon'' (real)');
    grid on;
    title('Soil Dielectric Decomposition');

    % Optional: show intrinsic water & clay real parts
    semilogx(freqs, real(eps_water_all), 'b--', 'LineWidth', 1.2);
    semilogx(freqs, real(eps_clay_all),  'r--', 'LineWidth', 1.2);
    legend({'Soil total', 'Water (intrinsic)', 'Clay (intrinsic)'}, ...
           'Location','best');

    % --- Plot imaginary part decomposition ---
    subplot(2,1,2);
    hold on;

    % Total soil loss (note: with your convention eps_im is negative for loss)
    semilogx(freqs, eps_im_total_all, 'k-', 'LineWidth', 1.8);

    % Conduction-only imaginary (from sigma_bulk)
    semilogx(freqs, eps_im_cond_all, 'g--', 'LineWidth', 1.2);

    % Dipolar-only imaginary (water+clay in the mixture)
    semilogx(freqs, eps_im_dip_all, 'm-.', 'LineWidth', 1.4);

    % Intrinsic water/clay imaginary parts (not mixed, just material response)
    semilogx(freqs, imag(eps_water_all), 'b:', 'LineWidth', 1.4);
    semilogx(freqs, imag(eps_clay_all),  'r:', 'LineWidth', 1.4);

    ylabel('\epsilon'''' (imag)');
    xlabel('Frequency [Hz]');
    grid on;

    legend({'Soil total', ...
            'Conduction (soil)', ...
            'Dipolar (soil mix)', ...
            'Water (intrinsic)', ...
            'Clay (intrinsic)'}, ...
           'Location','best');
end

function [eps_complex, eps_std] = invert_stuchly2pt_uncertainty(s_air, s_water, s_sample, freq, s_std)
% Invert complex permittivity using Stuchly's 2-point method with uncertainty
%
% Inputs:
%   s_air     - complex S11 in air
%   s_water   - complex S11 in water
%   s_sample  - complex S11 in sample
%   freq      - frequency vector
%   s_std     - standard deviation of S11 magnitude (assumed same for all)
%
% Outputs:
%   eps_complex - estimated complex permittivity
%   eps_std     - propagated standard deviation

    % Generate water permittivity
    % Coefficeients Derived from Fit to Kaatze(1989)
    % KaatzeData = readtable('E:\WCP\0925campaign\fieldfox\data\Kaatze_1989_data.xlsx');
    % % Water Temperature
    % Tmeas = 23.888;
    % % Extract variables
    % T = kaatzeData.T;
    % eps_s = kaatzeData.eps_s;
    % tau_ps = kaatzeData.tau;   % Relaxation time in picoseconds
    % % Convert τ to seconds (optional, depending on use case)
    % tau_s = tau_ps * 1e-12;
    % % --- Fit third-order polynomials ---
    % poly_epss = polyfit(T, eps_s, 3);
    % poly_tau  = polyfit(T, tau_s, 3);
    % % Evaluate fits
    % T_fit = linspace(min(T), max(T), 200);
    % epss_fit = polyval(poly_epss, T_fit);
    % tau_fit  = polyval(poly_tau, T_fit);

    % Generate water permittivity
    Tmeas = 23.888;
    epssCoefs = [88.119019607843190,-0.484234778121778,0.004464654282766,-2.393188854489220e-05]';
    tauCoefs  = [1.736936274509802e-11,-4.823117475060188e-13,5.488493292053644e-15,-2.366013071895416e-17]';
    G = [ones(numel(Tmeas),1),Tmeas,Tmeas.^2,Tmeas.^3];
    epss_fit = G * epssCoefs;
    tau_fit  = G * tauCoefs;
    params = struct('eps_s', epss_fit, 'eps_inf', 4.8, 'tau', tau_fit, ...
                    'alpha', 0.001, 'sigma', 0);
    eps_water = water_permittivity_colecole(freq, params);

    % Stuchly inversion
    r_air   = s_air;
    r_water = s_water;
    delta_s = r_water - r_air;
    ratio   = (s_sample - r_air) ./ delta_s;
    eps_complex = (eps_water - 1) .* ratio + 1;

    % --- Uncertainty propagation ---
    % Partial derivatives
    dE_dSsample = (eps_water - 1) ./ delta_s;
    dE_dSair    = -(eps_water - 1) .* (1 + (s_sample - r_air) ./ delta_s) ./ delta_s;
    dE_dSwater  = -(eps_water - 1) .* (s_sample - r_air) ./ delta_s.^2;

    % % Total differential: quadrature sum of uncertainties
    % eps_var = abs(dE_dSsample).^2 * s_std^2 + ...
    %           abs(dE_dSair).^2    * s_std^2 + ...
    %           abs(dE_dSwater).^2  * s_std^2;
    % eps_std = sqrt(eps_var);  % Standard deviation

    % Total differential: quadrature sum of uncertainties
    eps_var_real = real(dE_dSsample).^2 * s_std^2 + ...
        real(dE_dSair).^2    * s_std^2 + ...
        real(dE_dSwater).^2  * s_std^2;

    eps_var_imag = imag(dE_dSsample).^2 * s_std^2 + ...
        imag(dE_dSair).^2    * s_std^2 + ...
        imag(dE_dSwater).^2  * s_std^2;

    eps_std_real = sqrt(eps_var_real);
    eps_std_imag = sqrt(eps_var_imag);
    eps_std = complex(eps_std_real,eps_std_imag);

end


function [eps_complex, freq] = invert_stuchly2pt(s_air, s_water, s_sample, freq)
% invert_stuchly2pt - 2-point dielectric probe inversion using air and water
%
% Inputs:
%   s_air    - complex S11 measurement in air (vector)
%   s_water  - complex S11 measurement in DI water (~25 °C)
%   s_sample - complex S11 of unknown sample (same frequency vector)
%   freq     - frequency vector in Hz
%
% Outputs:
%   eps_complex - complex permittivity of sample (vector)
%   freq        - frequency vector in Hz

% Known permittivity of water at 25 °C (approximate)
% eps_water = 78.3 - 10j;  % adjust as needed
% eps_water = water_permittivity_tiuri80(freq, 298.15);
%   params    - struct with fields:
%               eps_s   : static permittivity
%               eps_inf : infinite frequency permittivity
%               tau     : relaxation time (seconds)
%               alpha   : Cole-Cole broadening parameter (0 for Debye)
%               sigma   : ionic conductivity (S/m)

% Water Temperature Measured at Calibration
Tmeas = 23.888;Tmeas = Tmeas(:);
% % Kaatze (1989) Data and Poly Fit
% kaatzeData = readtable('E:\WCP\0925campaign\fieldfox\data\Kaatze_1989_data.xlsx');
% % Water Temperature
% % Extract variables
% T = kaatzeData.T;
% eps_s = kaatzeData.eps_s;
% tau_ps = kaatzeData.tau;   % Relaxation time in picoseconds
% 
% % Convert τ to seconds (optional, depending on use case)
% tau_s = tau_ps * 1e-12;
% 
% % --- Fit third-order polynomials ---
% poly_epss = polyfit(T, eps_s, 3);
% poly_tau  = polyfit(T, tau_s, 3);
% % Evaluate fits
% T_fit = linspace(min(T), max(T), 200);
% epss_fit = polyval(poly_epss, T_fit);
% tau_fit  = polyval(poly_tau, T_fit);
% epss_fit = polyval(poly_epss, Tmeas);
% tau_fit  = polyval(poly_tau, Tmeas);

% Hard Coded Coefficients Developed From Kaatze (1989) 10 - 90 C
epssCoefs = [88.119019607843190,-0.484234778121778,0.004464654282766,-2.393188854489220e-05]';
tauCoefs = [1.736936274509802e-11,-4.823117475060188e-13,5.488493292053644e-15,-2.366013071895416e-17]';
% 3rd Order Polynomial
G = [ones(numel(Tmeas),1),Tmeas,Tmeas.^2,Tmeas.^3];
% Parameter Evaluation
epss_fit = G*epssCoefs;
tau_fit = G*tauCoefs;

% params.eps_s = 78.3; % 25 C
params.eps_s = epss_fit;
params.eps_inf = 4.8;
% params.tau = 8.27e-12; % 25 C
params.tau = tau_fit;
params.alpha = 0.001;
params.sigma = 0.000;
eps_water = water_permittivity_colecole(freq, params);

% Normalize calibration measurements
r_air   = s_air;
r_water = s_water;

% Apply 2-point inversion formula
eps_complex = (eps_water - 1) .* (s_sample - r_air) ./ (r_water - r_air) + 1;

end

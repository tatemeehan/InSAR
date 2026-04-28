function [eps_complex, freq] = invert_stuchly3pt(s_air, s_short, s_water, s_sample, freq)
% invert_stuchly3pt - Stuchly & Stuchly dielectric probe inversion from 3-point cal
%
% Inputs:
%   s_air    - complex S11 measurement in air (vector)
%   s_short  - complex S11 measurement on short (metal plate)
%   s_water  - complex S11 measurement in DI water (~25 °C)
%   s_sample - complex S11 of unknown sample (same frequency vector)
%   freq     - frequency vector in Hz
%
% Outputs:
%   eps_complex - complex permittivity of sample (vector)
%   freq        - frequency vector in Hz

% Known permittivity of water at 25 °C (approximate)
% eps_water = 78.3 - 10j;  % adjust if needed
params.eps_s = 78.3;
params.eps_inf = 4.8;
params.tau = 8.27e-12;
params.alpha = 0.35;
params.sigma = 0.005;
eps_water = water_permittivity_colecole(freq, params);

% Normalize calibration measurements
r_air   = s_air;
r_short = s_short;
r_water = s_water;

% Coefficients A and B from calibration
A = (r_short - r_water) ./ (r_air - r_water);
B = (eps_water - 1) ./ A;

% Invert sample
eps_complex = B .* (s_sample - r_water) ./ (r_air - s_sample) + 1;

end

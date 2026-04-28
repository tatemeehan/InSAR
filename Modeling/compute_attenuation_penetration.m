function [alpha, delta, beta] = compute_attenuation_penetration(freq, eps_complex)
% compute_attenuation_penetration - Computes attenuation coefficient and penetration depth
%
% Inputs:
%   freq         - Frequency vector [Hz]
%   eps_complex  - Complex permittivity (ε = ε' - jε'') vector
%
% Outputs:
%   alpha        - Attenuation coefficient (Np/m)
%   delta        - Penetration depth (m)
%   freq         - Frequency vector [Hz]

% Constants
eps0 = 8.854187817e-12; % [F/m] vacuum permittivity
mu0 = 4*pi*1e-7;        % [H/m] vacuum permeability
c = 1 / sqrt(eps0 * mu0); % [m/s] speed of light in vacuum

% Angular frequency
omega = 2 * pi * freq(:);  % [rad/s]

% Complex refractive index: n = sqrt(ε_r)
n_complex = sqrt(eps_complex);

% Free-space wave number
k0 = omega / c;  % [rad/m]

% Complex propagation constant
k = k0 .* n_complex;

% Attenuation coefficient (Np/m) is the imaginary part of k
alpha = imag(k);

% Penetration depth δ = 1 / α (in meters)
delta = 1 ./ alpha;

% Optional: filter very low attenuation to avoid Inf
delta(alpha < 1e-10) = Inf;

beta = real(k);

end

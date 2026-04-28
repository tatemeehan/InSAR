function eps = water_permittivity_colecole(frequency, params)
% WATER_PERMITTIVITY_COLECOLE
% Computes complex permittivity of water using Cole-Cole model.
%
% Inputs:
%   frequency - frequency in Hz (scalar or vector)
%   params    - struct with fields:
%               eps_s   : static permittivity
%               eps_inf : infinite frequency permittivity
%               tau     : relaxation time (seconds)
%               alpha   : Cole-Cole broadening parameter (0 for Debye)
%               sigma   : ionic conductivity (S/m)
%
% Output:
%   eps - complex permittivity

eps0 = 8.854187817e-12;  % vacuum permittivity [F/m]

% Unpack parameters
eps_s = params.eps_s;
eps_inf = params.eps_inf;
tau = params.tau;
alpha = params.alpha;
sigma = params.sigma;

omega = 2 * pi * frequency;

jomega_tau = (-1i * omega * tau).^(1 - alpha);

eps = eps_inf + (eps_s - eps_inf) ./ (1 + jomega_tau) + sigma ./ (-1i * omega * eps0);

end

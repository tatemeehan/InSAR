function [kz, beta, alpha, k0, kx, lambda, theta_t] = kz_from_eps(eps_r, theta_i_deg, f_Hz)
% Assumes eps_r has loss as +imag (eps'' > 0).
% Enforces imag(kz) >= 0 so alpha is positive.
% theta_t is a *complex* transmission angle diagnostic satisfying sin(theta_t)=kx/(k0*sqrt(eps_r)).

c = 299792458;
lambda = c / f_Hz;
k0 = 2*pi / lambda;

theta = deg2rad(theta_i_deg);

% Tangential wavenumber set by incidence in air (eps=1). Conserved across layers.
kx = k0 * sin(theta);

% Vertical component in the medium
kz = sqrt(k0^2 * eps_r - kx^2);

% Enforce sqrt branch so attenuation is positive
if imag(kz) < 0
    kz = -kz;
end

beta  = real(kz);
alpha = imag(kz);

% Diagnostic: complex transmission angle in this medium
n = sqrt(eps_r);
theta_t = asin(kx ./ (k0 .* n));
end
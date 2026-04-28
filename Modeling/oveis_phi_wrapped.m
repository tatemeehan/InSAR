function phi_wrapped = oveis_phi_wrapped(deltaSWE_mm, theta_deg, f_Hz, rho_kgm3)
% Canonical refraction-only wrapped phase vs SWE change.
% deltaSWE_mm : SWE change in mm water equivalent (numerically = kg/m^2)
% theta_deg   : incidence angle in air (deg)
% f_Hz        : radar frequency (Hz)
% rho_kgm3    : snow density (kg/m^3)

c0 = 299792458;
k0 = 2*pi*f_Hz/c0;                     % 2π/λ0  [1/m]

rho_gcm3 = rho_kgm3/1000;              % kg/m^3 -> g/cm^3
eps_s = eps_dry_maetzler1987(rho_gcm3);

th = deg2rad(theta_deg);

C = cos(th) - sqrt(eps_s - sin(th).^2); % unitless

% Convert SWE(mm) -> SWE(m)
deltaSWE_m = deltaSWE_mm / 1000;

rho_w = 1000;                          % kg/m^3

% Δd = (ρw/ρ) ΔSWE
delta_d = (rho_w/rho_kgm3) * deltaSWE_m;

phi_unwrapped = -2 * k0 .* C .* delta_d;   % radians
phi_wrapped   = wrapToPi(phi_unwrapped);

phi_wrapped = phi_wrapped(:);
end

function eps_r = eps_dry_maetzler1987(rho_gcm3)
rho = rho_gcm3;
eps_r = zeros(size(rho));

m = (rho < 0.4);
eps_r(m)  = 1 + 1.5995*rho(m) + 1.861*rho(m).^3;
eps_r(~m) = ((1 - rho(~m)./0.917) + 1.4759*(rho(~m)./0.917)).^3;
end

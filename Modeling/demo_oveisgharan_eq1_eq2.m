function dphi2 = demo_oveisgharan_eq1_eq2()
% Demo for Oveisgharan et al. (2024) Eq. (1) and Eq. (2)

% Sentinel-1 C-band center frequency used in the paper
% f0 = 5.405e9;                 % Hz
f0 = 1.3e9;                 % Hz

theta = 35;                   % deg incidence angle
rho_gcm3 = 0.30;              % g/cm^3 (300 kg/m^3)
dSWE_cm = 2.0;                % cm SWE change

% Compute dry-snow permittivity from Mätzler (1987) as used in the paper
eps_snow = eps_dry_maetzler1987(rho_gcm3);

% Convert ΔSWE -> Δd (snow depth change), using ΔSWE = (rho/rho_w) * Δd
rho_w = 1.0;                  % g/cm^3
dSWE_m = dSWE_cm / 100;        % m
dd_m = dSWE_m * (rho_w / rho_gcm3);

% Eq. (1): Δphi from Δd and eps
dphi1 = dphi_eq1(theta, f0, eps_snow, dd_m);

% Eq. (2): Δphi from ΔSWE and rho via C(theta,rho)
dphi2 = dphi_eq2(theta, f0, rho_gcm3, dSWE_m);

fprintf("eps_dry(rho=%.3f g/cm^3) = %.6f\n", rho_gcm3, eps_snow);
fprintf("ΔSWE = %.2f cm  =>  Δd = %.3f m\n", dSWE_cm, dd_m);
fprintf("Eq(1) Δphi = %.4f rad (%.2f deg)\n", dphi1, rad2deg(dphi1));
fprintf("Eq(2) Δphi = %.4f rad (%.2f deg)\n", dphi2, rad2deg(dphi2));
fprintf("Difference Eq(2)-Eq(1) = %.4e rad\n", dphi2-dphi1);

% SWE ambiguity (Δphi = 2π) for this theta and rho using Eq(2)
dSWE_amb_m = (2*pi) / abs(dphi_eq2(theta, f0, rho_gcm3, 1.0)); % response to 1 m SWE
fprintf("SWE ambiguity (2π) ~ %.2f cm at θ=%.1f°, ρ=%.2f g/cm^3\n", ...
        100*dSWE_amb_m, theta, rho_gcm3);

%% Let's Loop over a SWE change
dSWE = 0:0.001:1;
for kk = 1:numel(dSWE)
dphi2(kk) = dphi_eq2(theta, f0, rho_gcm3, dSWE(kk));
end
end

function dphi = dphi_eq1(theta_deg, f_Hz, eps_snow, dd_m)
% Oveisgharan et al. Eq. (1): Δphi from snow depth change Δd
% Δphi = -2 * kappa_i * (cosθ - sqrt(eps - sin^2θ)) * Δd
%
% Notes:
% - theta is incidence angle in air
% - eps_snow is (typically real) permittivity for dry snow
% - dd_m is snow depth change (m)
%
% kappa_i is an "incidence wavenumber" with units 1/m (paper states m^-1).
% A natural choice is the free-space wavenumber: k0 = 2π / λ0 = 2π f / c.

c0 = 299792458;                 % m/s
kappa_i = 2*pi*f_Hz/c0;         % 1/m (free-space wavenumber)

th = deg2rad(theta_deg);
term = cos(th) - sqrt(eps_snow - sin(th).^2);

dphi = -2 * kappa_i .* term .* dd_m;
end

function dphi = dphi_eq2(theta_deg, f_Hz, rho_gcm3, dSWE_m)
% Oveisgharan et al. Eq. (2): Δphi from ΔSWE using C(theta,rho)
% Δphi = -2 * kappa_i * C(theta,rho) * (rho_w/rho) * ΔSWE
% where C(theta,rho) = cosθ - sqrt(eps(rho) - sin^2θ)
%
% eps(rho) follows Mätzler (1987) dry-snow model as written in the paper.

c0 = 299792458;                 % m/s
kappa_i = 2*pi*f_Hz/c0;         % 1/m

rho_w = 1.0;                    % g/cm^3
eps_snow = eps_dry_maetzler1987(rho_gcm3);

th = deg2rad(theta_deg);
C = cos(th) - sqrt(eps_snow - sin(th).^2);

dphi = -2 * kappa_i .* C .* (rho_w./rho_gcm3) .* dSWE_m;
end

function eps_r = eps_dry_maetzler1987(rho_gcm3)
% Dry-snow real permittivity ε'(ρ) from Mätzler (1987) as written in the paper:
% For rho < 0.4 g/cm^3:
%   eps = 1 + 1.5995*rho + 1.861*rho^3
% For rho >= 0.4 g/cm^3:
%   eps = ((1 - rho/0.917) + 1.4759*(rho/0.917))^3
%
% (ρ in g/cm^3)

rho = rho_gcm3;

if any(rho < 0.0)
    error("rho must be nonnegative.");
end

eps_r = zeros(size(rho));

mask = (rho < 0.4);
eps_r(mask) = 1 + 1.5995*rho(mask) + 1.861*rho(mask).^3;

mask2 = ~mask;
eps_r(mask2) = ((1 - rho(mask2)./0.917) + 1.4759*(rho(mask2)./0.917)).^3;
end

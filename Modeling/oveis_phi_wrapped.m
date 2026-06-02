
function phi_wrapped = oveis_phi_wrapped(deltaSWE_mm, theta_deg, f_Hz, rho_kgm3, varargin)
%OVEIS_PHI_WRAPPED Canonical dry-snow refraction-only wrapped phase.
%
% Inputs
%   deltaSWE_mm : SWE change [mm water equivalent]
%                 Numerically equivalent to kg/m^2.
%   theta_deg   : incidence angle in air [deg]
%   f_Hz        : radar frequency [Hz]
%   rho_kgm3    : snow density [kg/m^3]
%
% Name-value options
%   'Sign'       : +1 or -1, default +1
%   'WrapOutput' : true/false, default true
%
% Output
%   phi_wrapped : wrapped phase [rad]

p = inputParser;
p.addParameter('Sign', +1, @isnumeric);
p.addParameter('WrapOutput', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

sgn = p.Results.Sign;
wrapOutput = logical(p.Results.WrapOutput);

c0 = 299792458;
k0 = 2*pi*f_Hz/c0;

rho_gcm3 = rho_kgm3 ./ 1000;
eps_s = eps_dry_maetzler1987(rho_gcm3);

th = deg2rad(theta_deg);

C = cos(th) - sqrt(eps_s - sin(th).^2);

% SWE mm -> SWE m
deltaSWE_m = deltaSWE_mm ./ 1000;

rho_w = 1000;

% Equivalent snow-depth change [m]
delta_d = (rho_w ./ rho_kgm3) .* deltaSWE_m;

phi_unwrapped = sgn .* (-2 .* k0 .* C .* delta_d);

if wrapOutput
    phi_wrapped = angle(exp(1i .* phi_unwrapped));
else
    phi_wrapped = phi_unwrapped;
end

phi_wrapped = single(phi_wrapped);

end


function eps_r = eps_dry_maetzler1987(rho_gcm3)
rho = rho_gcm3;
eps_r = zeros(size(rho), 'like', rho);

m = rho < 0.4;

eps_r(m)  = 1 + 1.5995 .* rho(m) + 1.861 .* rho(m).^3;

eps_r(~m) = ((1 - rho(~m)./0.917) + ...
    1.4759 .* (rho(~m)./0.917)).^3;

end

% function phi_wrapped = oveis_phi_wrapped(deltaSWE_mm, theta_deg, f_Hz, rho_kgm3)
% % Canonical refraction-only wrapped phase vs SWE change.
% % deltaSWE_mm : SWE change in mm water equivalent (numerically = kg/m^2)
% % theta_deg   : incidence angle in air (deg)
% % f_Hz        : radar frequency (Hz)
% % rho_kgm3    : snow density (kg/m^3)
% 
% c0 = 299792458;
% k0 = 2*pi*f_Hz/c0;                     % 2π/λ0  [1/m]
% 
% rho_gcm3 = rho_kgm3/1000;              % kg/m^3 -> g/cm^3
% eps_s = eps_dry_maetzler1987(rho_gcm3);
% 
% th = deg2rad(theta_deg);
% 
% C = cos(th) - sqrt(eps_s - sin(th).^2); % unitless
% 
% % Convert SWE(mm) -> SWE(m)
% deltaSWE_m = deltaSWE_mm / 1000;
% 
% rho_w = 1000;                          % kg/m^3
% 
% % Δd = (ρw/ρ) ΔSWE
% delta_d = (rho_w/rho_kgm3) * deltaSWE_m;
% 
% phi_unwrapped = -2 * k0 .* C .* delta_d;   % radians
% phi_wrapped   = wrapToPi(phi_unwrapped);
% 
% phi_wrapped = phi_wrapped(:);
% end
% 
% function eps_r = eps_dry_maetzler1987(rho_gcm3)
% rho = rho_gcm3;
% eps_r = zeros(size(rho));
% 
% m = (rho < 0.4);
% eps_r(m)  = 1 + 1.5995*rho(m) + 1.861*rho(m).^3;
% eps_r(~m) = ((1 - rho(~m)./0.917) + 1.4759*(rho(~m)./0.917)).^3;
% end

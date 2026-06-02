function phi_wrapped = oveis_phi_from_passes(H1, rho1, H2, rho2, theta_deg, f_Hz, varargin)
%OVEIS_PHI_FROM_PASSES Dry-snow refraction phase using both pass states.
%
% H1,H2       : snow depth [m]
% rho1,rho2   : snow density [kg/m^3]
% theta_deg   : incidence angle [deg]
% f_Hz        : frequency [Hz]

p = inputParser;
p.addParameter('Sign', +1, @isnumeric);
p.addParameter('WrapOutput', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

sgn = p.Results.Sign;
wrapOutput = logical(p.Results.WrapOutput);

c0 = 299792458;
k0 = 2*pi*f_Hz/c0;

theta = deg2rad(theta_deg);

eps1 = eps_dry_maetzler1987(rho1 ./ 1000);
eps2 = eps_dry_maetzler1987(rho2 ./ 1000);

C1 = cos(theta) - sqrt(eps1 - sin(theta).^2);
C2 = cos(theta) - sqrt(eps2 - sin(theta).^2);

% Snow refraction contribution for each pass
phi1 = -2 .* k0 .* C1 .* H1;
phi2 = -2 .* k0 .* C2 .* H2;

phi_unwrapped = sgn .* (phi2 - phi1);

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
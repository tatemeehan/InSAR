function err = phase_bias_to_oveis_error(biasRad, theta_deg, f_Hz, rho_kgm3)
%PHASE_BIAS_TO_OVEIS_ERROR Convert refraction phase bias to equivalent dH/SWE error.
%
% biasRad    : scalar phase bias [rad]
% theta_deg  : incidence angle raster [deg]
% f_Hz       : radar frequency [Hz]
% rho_kgm3   : density scalar or raster [kg/m^3]

c0 = 299792458;
k0 = 2*pi*f_Hz/c0;

rho_gcm3 = rho_kgm3 ./ 1000;
eps_s = eps_dry_maetzler1987_local(rho_gcm3);

theta = deg2rad(theta_deg);
C = cos(theta) - sqrt(eps_s - sin(theta).^2);

K = -2 .* k0 .* C;  % rad/m snow-depth change

dH_m = biasRad ./ K;
SWE_mm = rho_kgm3 .* dH_m;

err = struct();
err.K_rad_per_m = K;
err.dH_m = dH_m;
err.SWE_mm = SWE_mm;
err.abs_dH_m_median = median(abs(dH_m(:)), 'omitnan');
err.abs_SWE_mm_median = median(abs(SWE_mm(:)), 'omitnan');
err.dH_m_median = median(dH_m(:), 'omitnan');
err.SWE_mm_median = median(SWE_mm(:), 'omitnan');

end

function eps_r = eps_dry_maetzler1987_local(rho_gcm3)
rho = rho_gcm3;
eps_r = zeros(size(rho), 'like', rho);

m = rho < 0.4;
eps_r(m)  = 1 + 1.5995 .* rho(m) + 1.861 .* rho(m).^3;
eps_r(~m) = ((1 - rho(~m)./0.917) + 1.4759 .* (rho(~m)./0.917)).^3;
end
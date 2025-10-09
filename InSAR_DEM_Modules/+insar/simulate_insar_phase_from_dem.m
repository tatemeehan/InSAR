function [phi_unwrapped, phi_wrapped] = simulate_insar_phase_from_dem(dem, Bperp, r_slant, lambda, theta)
% simulate_insar_phase_from_dem
% Computes simulated interferometric phase from a DEM
%
% INPUTS:
%   dem     - 2D DEM (elevation in meters)
%   Bperp   - 2D map of perpendicular baseline (m)
%   r_slant - 2D or scalar slant range map (m)
%   lambda  - radar wavelength (m)
%   theta   - optional: incidence angle in radians (same size as dem or scalar)
%
% OUTPUT:
%   phi_topo - Simulated InSAR phase (radians)

if nargin < 5 || isempty(theta)
    cos_theta = 1;  % Assume vertical projection
else
    if any(theta(:) > pi)
        % Convert to Radians
        theta = deg2rad(theta);
    end
    cos_theta = cos(theta);
end


% Phase due to terrain height
% phi_unwrapped = (4 * pi / lambda) .* (Bperp ./ r_slant) .* dem .* cos_theta;
phi_unwrapped = -(4 * pi / lambda) .* (Bperp ./ r_slant) .* dem ./ cos_theta;
% phi_unwrapped = -(4 * pi / lambda) .* (Bperp ./ r_slant) .* 1 ./ cos_theta;
% phi_unwrapped = -(4 * pi / lambda) .* (Bperp ./ r_slant) .* (dem-min(dem(:))) ./ cos_theta;



% Wrapped Phase
phi_wrapped = mod(phi_unwrapped + pi, 2*pi) - pi;

end

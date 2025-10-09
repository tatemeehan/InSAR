function dem_est = insar_phase_to_dem(phi, lambda, r, Bperp, theta)
% Convert unwrapped phase to elevation using InSAR geometry
%
% Inputs:
%   phi    - Unwrapped phase (radians), 2D matrix
%   lambda - Radar wavelength (m), scalar
%   r      - Slant range (m), same size as phi or scalar
%   Bperp  - Perpendicular baseline (m), same size or scalar
%   theta  - Incidence angle (rad), same size or scalar
%
% Output:
%   dem_est - Estimated elevation (m)

% Use defaults if any are scalar
if isscalar(r),      r = r * ones(size(phi)); end
if isscalar(Bperp),  Bperp = Bperp * ones(size(phi)); end
if isscalar(theta),  theta = theta * ones(size(phi)); end
if any(theta(:) > pi)
    % Convert to Radians
    theta = deg2rad(theta);
end

% Elevation estimation
% NOTE: Negative sign convention is consistent with LOS vector computation
dem_est = -((phi .* lambda .* r) ./ (4 * pi .* Bperp .* cos(theta)));

end

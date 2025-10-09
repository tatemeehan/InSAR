function [kz] = compute_vertical_wavenumber(lambda,baseline,slantRange,incidence)
% compute_vertical_wavenumber uses the geometry from a baseline insar pair
% to compute the vertical component of the wave scattering for coherence
% analysis.
% Inputs:       lambda - the wavelength (m)
%               baseline - the perpendicular baseline
%               slantRange - the slantRange of the master SLC
%               incidence - the incidence angle of the master SLC
%
% Outputs:      kz - the vertical wavenumber
if any(incidence(:)>3.2)
    % Convert to Radians
    incidence = deg2rad(incidence);
end
kz = (4*pi/lambda) .* (baseline ./ (slantRange .* max(sin(incidence),1e-3)));

end
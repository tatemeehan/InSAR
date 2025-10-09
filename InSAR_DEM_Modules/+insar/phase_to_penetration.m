function penetration = phase_to_penetration(phase_unwrapped, lambda, Bperp, r_slant, incidence)
%PHASE_TO_PENETRATION Estimate penetration depth or surface bias from InSAR phase
%
%   penetration = phase_to_penetration(phase_unwrapped, lambda, Bperp, r_slant, incidence)
%
%   Converts interferometric phase (typically unwrapped) to apparent vertical
%   displacement (e.g., signal penetration or bias in scattering surface), assuming
%   a DEM was used to flatten the phase during SAR backprojection.
%
%   INPUTS:
%       phase_unwrapped - Unwrapped interferometric phase (radians)
%       lambda           - Radar wavelength (meters)
%       Bperp            - Perpendicular baseline (meters)
%       r_slant          - Slant range distance from antenna to pixel (meters)
%       incidence        - Local incidence angle (degrees)
%
%   OUTPUT:
%       penetration      - Apparent penetration depth or surface bias (meters)
%
%   NOTES:
%       This function assumes the SAR processor used terrain flattening via a DEM
%       (e.g., time-domain backprojection), such that the phase no longer encodes
%       topographic variation but residual path delays. These can arise from
%       scattering beneath the DEM surface (e.g., in snow or vegetation).

    % Ensure all inputs are the same size
    if ~isequal(size(phase_unwrapped), size(Bperp), size(r_slant), size(incidence))
        error('All input arrays must be the same size.');
    end

    % Convert to radians if incidence is degrees
    if any(incidence(:) > pi)
        % Convert to Radians
        incidence = deg2rad(incidence);
    end
    cos_theta = cos(incidence);

    % Phase to path delay to vertical penetration estimate
    penetration = (lambda .* r_slant) ./ (4 * pi .* Bperp .* cos_theta) .* phase_unwrapped;
end

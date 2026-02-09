
function penetration = phase_to_penetration(phase_unwrapped, lambda, Bperp, r_slant, incidence)
%PHASE_TO_PENETRATION Apparent vertical phase-center offset (m) from InSAR phase
%
%   penetration is an apparent vertical bias of the effective phase center
%   relative to the DEM reference surface (sign depends on your convention).

    if ~isequal(size(phase_unwrapped), size(Bperp), size(r_slant), size(incidence))
        error('All input arrays must be the same size.');
    end

    % degrees -> radians (your check is fine)
    if any(incidence(:) > pi)
        incidence = deg2rad(incidence);
    end

    sin_theta = sin(incidence);

    % Optional masks to avoid blow-ups
    Bmin = 0.1;      % m, tune for your platform
    stmin = 0.05;    % avoid near-grazing issues
    mask = isfinite(phase_unwrapped) & isfinite(Bperp) & isfinite(r_slant) & isfinite(sin_theta) ...
         & abs(Bperp) > Bmin & sin_theta > stmin;

    penetration = nan(size(phase_unwrapped));
    penetration(mask) = (lambda .* r_slant(mask) .* sin_theta(mask)) ./ (4*pi .* Bperp(mask)) .* phase_unwrapped(mask);
end
% function penetration = phase_to_penetration(phase_unwrapped, lambda, Bperp, r_slant, incidence)
% %PHASE_TO_PENETRATION Estimate penetration depth or surface bias from InSAR phase
% %
% %   penetration = phase_to_penetration(phase_unwrapped, lambda, Bperp, r_slant, incidence)
% %
% %   Converts interferometric phase (typically unwrapped) to apparent vertical
% %   displacement (e.g., signal penetration or bias in scattering surface), assuming
% %   a DEM was used to flatten the phase during SAR backprojection.
% %
% %   INPUTS:
% %       phase_unwrapped - Unwrapped interferometric phase (radians)
% %       lambda           - Radar wavelength (meters)
% %       Bperp            - Perpendicular baseline (meters)
% %       r_slant          - Slant range distance from antenna to pixel (meters)
% %       incidence        - Local incidence angle (degrees)
% %
% %   OUTPUT:
% %       penetration      - Apparent penetration depth or surface bias (meters)
% %
% %   NOTES:
% %       This function assumes the SAR processor used terrain flattening via a DEM
% %       (e.g., time-domain backprojection), such that the phase no longer encodes
% %       topographic variation but residual path delays. These can arise from
% %       scattering beneath the DEM surface (e.g., in snow or vegetation).
% 
%     % Ensure all inputs are the same size
%     if ~isequal(size(phase_unwrapped), size(Bperp), size(r_slant), size(incidence))
%         error('All input arrays must be the same size.');
%     end
% 
%     % Convert to radians if incidence is degrees
%     if any(incidence(:) > pi)
%         % Convert to Radians
%         incidence = deg2rad(incidence);
%     end
%     % cos_theta = cos(incidence);
%         sin_theta = sin(incidence);
% 
% 
%     % Phase to path delay to vertical penetration estimate
%     % penetration = (lambda .* r_slant) ./ (4 * pi .* Bperp .* cos_theta) .* phase_unwrapped;
%         penetration = (lambda .* r_slant) ./ (4 * pi .* Bperp .* sin_theta) .* phase_unwrapped;
% 
% end

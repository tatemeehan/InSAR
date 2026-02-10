function Bcrit = compute_critical_baseline(lambda, r_slant, incidence_deg, groundRangeResolution, n)
% COMPUTE_CRITICAL_BASELINE  Topographic critical baseline (range spectral shift)
%
% Bcrit = λ * r_slant / (n * groundRangeResolution * cos(θ))
%
% n = 2 is common in practice; allow override.

    if nargin < 5 || isempty(n), n = 2; end
    if any(incidence_deg(:) > pi)
    theta = deg2rad(incidence_deg);
    end
    Bcrit = (lambda .* r_slant) ./ (n .* groundRangeResolution .* cos(theta));
end
% function Bcrit = compute_critical_baseline(lambda, r_slant, incidence_deg, groundRangeResolution)
% %COMPUTE_CRITICAL_BASELINE Computes the topographic critical baseline
% %
% %   Bcrit = compute_critical_baseline(lambda, r_slant, incidence_deg, groundRangeResolution)
% %
% %   Inputs:
% %       lambda               - Radar wavelength [m]
% %       r_slant              - Slant range [m]
% %       incidence_deg        - Incidence angle [degrees]
% %       groundRangeResolution - Ground range pixel spacing [m]
% %
% %   Output:
% %       Bcrit - Critical baseline [m]
% 
%     theta_rad = deg2rad(incidence_deg);
%     Bcrit = (lambda .* r_slant .* cos(theta_rad)) ./ (2 .* groundRangeResolution);
% end
% function B_crit = compute_critical_baseline(lambda, r_slant, L_az)
% % compute_critical_baseline
% % Computes critical perpendicular baseline beyond which InSAR decorrelation occurs
% %
% % INPUTS:
% %   lambda   - Radar wavelength (m)
% %   r_slant  - Slant range to ground target (m)
% %   L_az     - Azimuth antenna length (m)
% %
% % OUTPUT:
% %   B_crit   - Critical perpendicular baseline (m)
% 
% % Check inputs
% if nargin < 3
%     error('Inputs: lambda, r_slant, L_az are required.');
% end
% 
% % Critical baseline from azimuth decorrelation
% B_crit = (r_slant .* lambda) ./ L_az;
% 
% end

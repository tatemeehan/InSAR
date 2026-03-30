function sigma_phi = coherence_to_phase_std(coh, nLooks)
%COHERENCE_TO_PHASE_STD Approximate interferometric phase std (rad)
%
% coh    : coherence magnitude [0,1]
% nLooks : effective number of looks

    if nargin < 2 || isempty(nLooks)
        nLooks = 1;
    end
    mask = isnan(coh);
    coh = max(min(coh, 0.9999), 1e-6);  % avoid division blow-up
    sigma_phi = sqrt((1 - coh.^2) ./ (2 * nLooks .* coh.^2));
    sigma_phi(mask) = NaN;
end
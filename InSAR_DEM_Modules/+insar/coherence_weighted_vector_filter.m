function phi_filtered = coherence_weighted_vector_filter(phi, gamma, winSize, sigma, alpha)
%COHERENCE_WEIGHTED_VECTOR_FILTER Apply coherence-weighted vector filter to wrapped phase
%
%   phi_filtered = coherence_weighted_vector_filter(phi, gamma, sigma, alpha, winSize)
%
%   Inputs:
%       phi     - Wrapped phase matrix (in radians, [-pi, pi])
%       gamma   - Magnitude coherence matrix (same size as phi)
%       sigma   - Standard deviation for Gaussian spatial weights [pixels]
%       alpha   - Exponent for coherence weighting (>= 1)
%       winSize - Window size (odd integer, e.g. 5 for a 5x5 kernel)
%
%   Output:
%       phi_filtered - Filtered wrapped phase matrix

% Ensure input validity
assert(isequal(size(phi), size(gamma)), 'phi and gamma must be the same size');
assert(mod(winSize, 2) == 1, 'winSize must be an odd integer');
assert(alpha >= 0, 'alpha must be non-negative');

% Precompute Gaussian spatial weights
halfWin = floor(winSize / 2);
[xg, yg] = meshgrid(-halfWin:halfWin, -halfWin:halfWin);
G = exp(-(xg.^2 + yg.^2) / (2 * sigma^2));
G = G / sum(G(:)); % normalize

% Preallocate output
phi_filtered = nan(size(phi));
[m, n] = size(phi);

% Pad phase and coherence to handle borders
phi_pad = padarray(phi, [halfWin, halfWin], 'replicate');
gamma_pad = padarray(gamma, [halfWin, halfWin], 'replicate');

% Loop over each pixel
for ii = 1:m
    for jj = 1:n
        % Skip if center pixel is NaN
        if isnan(phi(ii, jj))
            continue;
        end
        % Extract window patches
        patch_phi = phi_pad(ii:ii+2*halfWin, jj:jj+2*halfWin);
        patch_gamma = gamma_pad(ii:ii+2*halfWin, jj:jj+2*halfWin);
        patch_gamma = min(max(patch_gamma, 0), 1); % clip gamma

        % Convert to complex
        z = exp(1i * patch_phi);

        % Compute weights: spatial × coherence^α
        W = G .* (patch_gamma .^ alpha);
        wSum = sum(W(:));

        % SAFETY NET: handle zero total weight
        if wSum == 0
            phi_filtered(ii,jj) = patch_phi(halfWin+1, halfWin+1); % fallback to center value
            continue;
        end
        % Weighted vector sum
        W = W / wSum; % normalize
        z_sum = sum(W(:) .* z(:));
        phi_filtered(ii,jj) = angle(z_sum);
        
    end
end

end

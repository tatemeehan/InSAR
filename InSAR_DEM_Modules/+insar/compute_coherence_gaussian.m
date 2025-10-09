function coherence = compute_coherence_gaussian(slc1, slc2, sigma)
%COMPUTE_COHERENCE_GAUSSIAN Estimate coherence using Gaussian kernel
%
%   coherence = compute_coherence_gaussian(slc1, slc2, sigma)
%
%   Inputs:
%       slc1, slc2 - complex SLC images (must be same size)
%       sigma      - standard deviation of Gaussian kernel
%
%   Output:
%       coherence  - real-valued coherence image in [0, 1]

    % Check size
    if ~isequal(size(slc1), size(slc2))
        error('SLC images must be the same size.');
    end

    % Interferometric product
    cross_prod = slc1 .* conj(slc2);

    % Amplitudes
    amp1 = abs(slc1);
    amp2 = abs(slc2);

    % Apply Gaussian filtering
    win_cross = imgaussfilt(real(cross_prod), sigma) + 1i * imgaussfilt(imag(cross_prod), sigma);
    win_amp1 = imgaussfilt(amp1.^2, sigma);
    win_amp2 = imgaussfilt(amp2.^2, sigma);

    % Compute coherence
    coherence = abs(win_cross) ./ sqrt(win_amp1 .* win_amp2);

    % Clamp to [0,1]
    coherence = min(max(coherence, 0), 1);
end

function coh = compute_coherence_box(slc1, slc2, winSize)
% COMPUTE_COHERENCE computes local interferometric coherence
% Inputs:
%   slc1, slc2 - complex SLC images (same size)
%   winSize    - window size for coherence estimation (odd number)
% Output:
%   coh        - coherence image (real, 0 to 1)

    if nargin < 3
        winSize = 5; % default window
    end

    % Numerator: cross-correlation
    crossProd = slc1 .* conj(slc2);
    num = abs(movmean(real(crossProd), winSize, 1, 'Endpoints','discard'));
    num = movmean(num, winSize, 2, 'Endpoints','discard');

    % Denominator: geometric mean of powers
    pow1 = movmean(abs(slc1).^2, winSize, 1, 'Endpoints','discard');
    pow1 = movmean(pow1, winSize, 2, 'Endpoints','discard');

    pow2 = movmean(abs(slc2).^2, winSize, 1, 'Endpoints','discard');
    pow2 = movmean(pow2, winSize, 2, 'Endpoints','discard');

    denom = sqrt(pow1 .* pow2);

    % Coherence
    coh = num ./ (denom + eps); % avoid division by zero
end

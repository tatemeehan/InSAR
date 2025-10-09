function [phz, cor, ccor, cphz] = compute_interferogram(slc1, slc2, qualityThresh, filterSize, useRobustWeighting, sigma, alpha)
%COMPUTE_INTERFEROGRAM Compute wrapped phase and phase quality
%
%   [phz, cor] = compute_interferogram(slc1, slc2, qualityThresh, filterSize, ...
%                                      useCoherenceFilter, sigma, alpha)
%
%   Inputs:
%       slc1, slc2           - complex SLC data
%       qualityThresh        - threshold for coherence (default = 0.75)
%       filterSize           - window size for filtering and quality (default = 5)
%       useCoherenceFilter   - true to use coherence-weighted vector filter (default = false)
%       sigma                - Gaussian kernel std. dev. (default = filterSize / 2)
%       alpha                - coherence weighting exponent (default = 2)
%
%   Outputs:
%       phz - filtered wrapped interferometric phase (NaNs where low quality)
%       cor - coherence-like magnitude metric

    if nargin < 3 || isempty(qualityThresh), qualityThresh = 0.75; end
    if nargin < 4 || isempty(filterSize), filterSize = 5; end
    if nargin < 5 || isempty(useRobustWeighting), useRobustWeighting = false; end
    if nargin < 6 || isempty(sigma), sigma = filterSize / 2; end
    if nargin < 7 || isempty(alpha), alpha = 2; end

    if useRobustWeighting
        [cor, ccor, ~] = insar.compute_coherence_robust(slc1, slc2, 'gaussian', filterSize, sigma);
        phz = angle(ccor);
        qualityThresh = qualityThresh;
        phz(cor<qualityThresh) = NaN;
    else
        [cor, ccor, ~] = insar.compute_coherence(slc1, slc2, 'gaussian', filterSize, sigma);
        phz = angle(ccor);
        phz(cor<qualityThresh) = NaN;
    end

    [cphz, ~, ~] = insar.insar_phase_complex(slc1, slc2, filterSize, sigma);

    % if useCoherenceFilter
    %     [phz, cor, ccor] = insar.compute_interferogram_adaptive(slc1, slc2, ...
    %         qualityThresh, filterSize, sigma, alpha);
    % else
    %     [phz, cor, ccor] = insar.compute_interferogram_lockstep(slc1, slc2, ...
    %         qualityThresh, filterSize, sigma);
    % end

    % % Compute raw wrapped interferogram
    % interferogram = slc1 .* conj(slc2);
    % phzRaw = angle(interferogram);
    % [cor, ccor] = insar.compute_coherence(slc1, slc2, 'gaussian', filterSize, sigma);
    % if exist("ccor","var")
    %     if useCoherenceFilter
    %         [ccor, phz] = insar.coherence_weighted_complex_interferogram( ...
    %             ccor, filterSize, sigma, alpha);
    %     else
    %         % no extra coherence-weighted filtering: just use complex coherence
    %         phz   = angle(ccor);
    %         cor     = abs(ccor);
    %     end
    % else
    %     if useCoherenceFilter
    %         % Use coherence-weighted vector filtering
            % phz = insar.coherence_weighted_vector_filter(phzRaw, cor, filterSize, sigma, alpha);
    %     else
    %         % Use standard vector median filtering
    %         phz = insar.vector_median_filter(phzRaw, filterSize);
    %     end
    % end
    % 
    % % Apply threshold mask
    % phz(cor < qualityThresh) = NaN;
end
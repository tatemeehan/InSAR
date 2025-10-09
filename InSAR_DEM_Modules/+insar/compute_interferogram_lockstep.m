function [phz, cor, ccor] = compute_interferogram_lockstep( ...
    slc1, slc2, qualityThresh, winSize, sigma)

if nargin < 3 || isempty(qualityThresh), qualityThresh = 0.75; end
if nargin < 4 || isempty(winSize),       winSize       = 5;    end
if nargin < 5 || isempty(sigma),         sigma         = winSize/2; end

% One pass (your compute_coherence)
[~, ccor] = insar.compute_coherence(slc1, slc2, 'gaussian', winSize, sigma);

phz   = angle(ccor);
cor   = abs(ccor);    % == coh

phz(cor < qualityThresh) = NaN;
end

function [coh,ccoh, meta] = compute_coherence(slc1, slc2, kernelType, winSize, sigma)
%COMPUTE_COHERENCE Compute interferometric coherence between two SLCs.
%
%   [coh, meta] = compute_coherence(slc1, slc2, kernelType, winSize, sigma)
%
%   Inputs:
%       slc1, slc2   - complex SLC images (same size)
%       kernelType   - 'gaussian' (default) or 'boxcar'
%       winSize      - window size (odd integer, default = 7)
%       sigma        - standard deviation for Gaussian (default = winSize/6)
%
%   Output:
%       coh          - magnitude of coherence (0 to 1)
%       ccoh         - complex coherence
%       meta         - struct with parameters used

% Defaults
if nargin < 3 || isempty(kernelType), kernelType = 'gaussian'; end
if nargin < 4 || isempty(winSize), winSize = 7; end
if nargin < 5 || isempty(sigma), sigma = winSize / 6; end

if (mod(winSize,2)==0)
    winSize = winSize+1;
end

% Cross product and norms
crossProd = slc1 .* conj(slc2);
power1 = abs(slc1).^2;
power2 = abs(slc2).^2;

% Choose filter kernel
switch lower(kernelType)
    case 'gaussian'
        kernel = fspecial('gaussian', winSize, sigma);
    case 'boxcar'
        kernel = ones(winSize) / winSize^2;
    otherwise
        error('Unknown kernelType: use ''gaussian'' or ''boxcar''');
end

% Filter numerator and denominator terms
num = imfilter(crossProd, kernel, 'symmetric');
den1 = imfilter(power1, kernel, 'symmetric');
den2 = imfilter(power2, kernel, 'symmetric');

% Complex coherence (filtered complex interferogram)
ccoh = num ./ sqrt(den1 .* den2);
ccoh(den1==0 | den2==0) = NaN;

% Coherence magnitude
coh = abs(ccoh);

% Metadata
meta.kernelType = kernelType;
meta.winSize = winSize;
meta.sigma = sigma;
meta.kernel = kernel;

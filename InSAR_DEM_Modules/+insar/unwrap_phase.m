function insarData = unwrap_phase(insarData, method, qualityThresh)
%UNWRAP_PHASE Unwraps interferometric phase in insarData using given method.
%   insarData = unwrap_phase(insarData, method, qualityThresh)
%
%   method: 'multiseed' (default) | 'goldstein'
%   qualityThresh: minimum coherence for valid unwrapping (default 0.75)

if nargin < 2 || isempty(method), method = 'multiseed'; end
if nargin < 3 || isempty(qualityThresh), qualityThresh = 0.75; end

for k = 1:numel(insarData)
    phz = insarData(k).phase;
    cor = insarData(k).coherence;
    mask = cor >= qualityThresh;

    switch lower(method)
        case 'multiseed'
            % phzUnw = insar.unwrap_multiseed(phz, mask);
            phzUnw = insar.region_grow_unwrap_multiseed(phz, cor, qualityThresh);
        case 'goldstein'
            phzUnw = insar.unwrap_goldstein(phz);
        otherwise
            error('Unknown unwrap method: %s', method);
    end

    insarData(k).unwrapped = phzUnw;
end
end
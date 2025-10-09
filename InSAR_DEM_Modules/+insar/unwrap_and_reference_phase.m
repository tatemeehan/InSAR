function sarData = unwrap_and_reference_phase(sarData, crResult, unwrapOptions)
%UNWRAP_AND_REFERENCE_PHASE Unwraps and phase-references interferometric phase.
%
% Inputs:
%   sarData       - Struct array with fields:
%                     .interferogram.phase : Wrapped phase (radians)
%                     .interferogram.coherence : Coherence-like quality
%   crResult      - Struct from CRpower() with CR positions and weights
%   unwrapOptions - Struct with optional fields:
%                     .method     : 'multiseed' (default), 'goldstein', etc.
%                     .crMode     : 'powerWeighted' (default), 'unweighted', 'sceneMean', 'none'
%                     .qualityThresh : Minimum coherence to unwrap (default 0.75)
%                     .filterSize    : Filtering kernel size (default 5)

if nargin < 3, unwrapOptions = struct(); end
if ~isfield(unwrapOptions, 'method'), unwrapOptions.method = 'multiseed'; end
if ~isfield(unwrapOptions, 'crMode'), unwrapOptions.crMode = 'powerWeighted'; end
if ~isfield(unwrapOptions, 'qualityThresh'), unwrapOptions.qualityThresh = 0.75; end
if ~isfield(unwrapOptions, 'filterSize'), unwrapOptions.filterSize = 5; end

for i = 1:numel(sarData)
    for j = 1:numel(sarData(i).interferogram)
        phz = sarData(i).interferogram(j).phase;
        cor = sarData(i).interferogram(j).coherence;

        % Unwrap phase
        switch lower(unwrapOptions.method)
            case 'multiseed'
                phzUnw = insar.unwrap_multiseed(phz, cor >= unwrapOptions.qualityThresh);
            case 'goldstein'
                phzUnw = insar.unwrap_goldstein(phz);
            otherwise
                error('Unknown unwrap method: %s', unwrapOptions.method);
        end

        % Reference unwrapped phase
        switch lower(unwrapOptions.crMode)
            case 'powerweighted'
                refPhases = [];
                weights = [];
                for c = 1:numel(crResult)
                    w = crResult(c).Weights{i}(j);
                    if isempty(w) || ~isfield(w, 'indices'), continue; end
                    px = w.indices;
                    wt = w.value;
                    if any(isnan(phzUnw(px))), continue; end
                    refPhases(end+1) = sum(wt .* phzUnw(px));
                    weights(end+1) = sum(wt);
                end
                if ~isempty(refPhases)
                    referencePhase = sum(refPhases .* weights) / sum(weights);
                else
                    referencePhase = mean(phzUnw(cor >= unwrapOptions.qualityThresh), 'omitnan');
                end
            case 'unweighted'
                refPx = []; 
                for c = 1:numel(crResult)
                    w = crResult(c).Weights{i}(j);
                    if isempty(w) || ~isfield(w, 'indices'), continue; end
                    refPx = [refPx; w.indices(:)];
                end
                refPx = unique(refPx);
                referencePhase = mean(phzUnw(refPx), 'omitnan');
            case 'scenemean'
                referencePhase = mean(phzUnw(cor >= unwrapOptions.qualityThresh), 'omitnan');
            case 'none'
                referencePhase = 0;
            otherwise
                error('Unknown crMode: %s', unwrapOptions.crMode);
        end

        sarData(i).interferogram(j).unwrapped = phzUnw - referencePhase;
        sarData(i).interferogram(j).refPhase = referencePhase;
    end
end
end

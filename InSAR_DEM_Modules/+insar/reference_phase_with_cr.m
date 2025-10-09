function phzRef = reference_phase_with_cr(phzUnwrapped, crResult, sarData, pairIdx)
%REFERENCE_PHASE_WITH_CR Subtracts CR-based mean phase from unwrapped phase
%
% phzRef = reference_phase_with_cr(phzUnwrapped, crResult, sarData, pairIdx)
%
% Inputs:
%   phzUnwrapped : Unwrapped interferometric phase map (2D)
%   crResult     : Struct array of CR metadata with .Weights and .Index
%   sarData      : SAR data array (for accessing SLCs)
%   pairIdx      : [i, j] index into sarData(i).slc{j}
%
% Output:
%   phzRef       : Phase map with CR reference phase removed

% Validate inputs
if nargin < 4 || isempty(pairIdx) || ~isnumeric(pairIdx) || numel(pairIdx) ~= 2
    error('pairIdx must be a two-element array [i j] for SAR pair index.');
end

refPhases = nan(numel(crResult), 1);

for c = 1:numel(crResult)
    weightsEntry = crResult(c).Weights{pairIdx(1)};
    if numel(weightsEntry) < pairIdx(2) || isempty(weightsEntry(pairIdx(2)).value)
        continue;
    end

    weights = weightsEntry(pairIdx(2)).value;
    idxs = weightsEntry(pairIdx(2)).indices;
    slc1 = sarData(pairIdx(1)).slc{pairIdx(2)};
    phz = angle(slc1(idxs));

    % Compute weighted average phase via complex phasors
    refPhases(c) = angle(sum(weights .* exp(1i * phz)));
end

% Remove NaNs and take average
validPhz = refPhases(~isnan(refPhases));
if isempty(validPhz)
    warning('No valid CR reference phases found for this SLC pair. Returning original phase.');
    phzRef = phzUnwrapped;
    return;
end

refPhase = angle(mean(exp(1i * validPhz)));
phzRef = phzUnwrapped - refPhase;

end

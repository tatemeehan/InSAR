function [phzCorrected, cycleMap, correctionMap, regionLabels] = correct_phase_cycles_global(phzUnwrapped, coherence, minRegionSize)
%CORRECT_PHASE_CYCLES_GLOBAL Aligns phase regions using dominant cycle value (no propagation)
%
%   INPUTS:
%       phzUnwrapped   - Referenced unwrapped phase [radians]
%       coherence      - Coherence map to mask invalid pixels
%       minRegionSize  - Minimum region size to consider (default: 30 pixels)
%
%   OUTPUTS:
%       phzCorrected   - Phase map with global 2pi cycle alignment
%       cycleMap       - Map of cycle numbers before correction (floor(phz / 2pi))
%       correctionMap  - Map of 2pi cycles subtracted per pixel
%       regionLabels   - Connected component labels for regions

if nargin < 3, minRegionSize = 30; end

% Step 1: Compute cycle map
validMask = isfinite(phzUnwrapped) & coherence > 0.05;
cycleMap = floor(phzUnwrapped / (2*pi));
cycleMap(~validMask) = 0;

% Step 2: Label connected regions
[regionLabels, numRegions] = bwlabel(validMask, 8);

% Step 3: Compute mode cycle value per region
regionCycles = zeros(1, numRegions);
for i = 1:numRegions
    pix = (regionLabels == i);
    if nnz(pix) <1% minRegionSize./2
        regionCycles(i) = NaN;
        continue;
    end
    vals = cycleMap(pix);
    regionCycles(i) = mode(vals(:));
end

% Step 4: Determine reference cycle value
refCycle = mode(regionCycles(~isnan(regionCycles)));

% Step 5: Apply correction
phzCorrected = phzUnwrapped;
correctionMap = zeros(size(phzUnwrapped));

for i = 1:numRegions
    mask = (regionLabels == i);
    if isnan(regionCycles(i)), continue; end
    dCycle = regionCycles(i) - refCycle;
    correctionMap(mask) = dCycle;
    phzCorrected(mask) = phzUnwrapped(mask) - 2*pi*dCycle;
end

end

function [phzCorrected, cycleMap, correctionMap, regionLabels] = correct_phase_cycles_regionwise(phzUnwrapped, coherence, minRegionSize)
%CORRECT_PHASE_CYCLES_REGIONWISE Aligns 2π-offset phase islands into a continuous map.
%
%   INPUTS:
%       phzUnwrapped   - Referenced unwrapped phase [radians]
%       coherence      - Coherence map (used to mask invalid pixels)
%       minRegionSize  - Minimum pixel count to retain a region
%
%   OUTPUTS:
%       phzCorrected   - Phase map with region-wise 2π cycle alignment applied
%       cycleMap       - Initial integer cycle map from floor(phzUnwrapped / 2π)
%       correctionMap  - Map of number of 2π cycles subtracted
%       regionLabels   - Labeled regions used in the correction

if nargin < 3, minRegionSize = 30; end

% ------------------- Step 1: Prepare Masks ------------------- %
validMask = isfinite(phzUnwrapped) & coherence > 0.05;
cycleMap = floor(phzUnwrapped / (2*pi));
cycleMap(~validMask) = 0;

% ------------------- Step 2: Label Regions ------------------- %
[L_all, num_all] = bwlabel(cycleMap ~= 0, 8);
props_all = regionprops(L_all, phzUnwrapped, 'PixelIdxList', 'PixelValues', 'Area');

% Filter out small regions
isValid = arrayfun(@(r) r.Area > minRegionSize, props_all);
props = props_all(isValid);
L_filtered = zeros(size(L_all));
for i = 1:numel(props)
    L_filtered(props(i).PixelIdxList) = i;
end
regionLabels = L_filtered;
numValid = numel(props);

% If no valid regions, return early
if numValid == 0
    phzCorrected = phzUnwrapped;
    correctionMap = zeros(size(phzUnwrapped));
    return;
end

% Median phase per region
medPhase = arrayfun(@(r) median(r.PixelValues, 'omitnan'), props);

% ------------------- Step 3: Flood from All Seeds ------------------- %
phzCorrected = phzUnwrapped;
correctionMap = zeros(size(phzUnwrapped));
visited = false(1, numValid);
phzOffsets = zeros(1, numValid);

% 8-connectivity neighborhood
connKernel = strel('square', 11).Neighborhood;

for seed = 1:numValid
    if visited(seed), continue; end

    % Start new flood
    queue = seed;
    visited(seed) = true;

    while ~isempty(queue)
        curr = queue(1); queue(1) = [];

        % Binary mask for current region
        maskCurr = (regionLabels == curr);

        % Dilate to find neighbor pixels
        neighborMask = imdilate(maskCurr, connKernel) & ~maskCurr & validMask;
        neighborLabels = unique(regionLabels(neighborMask));
        neighborLabels = neighborLabels(neighborLabels > 0 & neighborLabels <= numValid);
        neighborLabels = neighborLabels(~visited(neighborLabels));

        for n = neighborLabels(:)'
            dPhi = medPhase(n) - medPhase(curr);
            % cycleOffset = floor(dPhi / (2*pi));
            cycleOffset = floor((dPhi + pi) / (2*pi));
            phzOffsets(n) = phzOffsets(curr) + cycleOffset;
            visited(n) = true;
            queue(end+1) = n;
        end
    end
end

% ------------------- Step 4: Apply Corrections ------------------- %
for i = 1:numValid
    pix = props(i).PixelIdxList;
    phzCorrected(pix) = phzUnwrapped(pix) - 2*pi*phzOffsets(i);
    correctionMap(pix) = phzOffsets(i);
end

end

% function [phzCorrected, cycleMap, correctionMap, regionLabels] = correct_phase_cycles_regionwise(phzUnwrapped, coherence, minRegionSize)
% %CORRECT_PHASE_CYCLES_REGIONWISE Aligns 2π-offset phase islands into a continuous map.
% %
% %   INPUTS:
% %       phzUnwrapped   - Referenced unwrapped phase [radians]
% %       coherence      - Coherence map (used to mask invalid pixels)
% %       minRegionSize  - Minimum pixel count to retain a region
% %
% %   OUTPUTS:
% %       phzCorrected   - Phase map with region-wise 2π cycle alignment applied
% %       cycleMap       - Initial integer cycle map from floor(phzUnwrapped / 2π)
% %       correctionMap  - Map of number of 2π cycles subtracted
% %       regionLabels   - Labeled regions used in the correction
% 
% if nargin < 3, minRegionSize = 30; end
% 
% % ------------------- Step 1: Prepare Masks ------------------- %
% validMask = isfinite(phzUnwrapped) & coherence > 0.05;
% cycleMap = floor(phzUnwrapped / (2*pi));
% cycleMap(~validMask) = 0;
% 
% % ------------------- Step 2: Label Regions ------------------- %
% [L_all, num_all] = bwlabel(cycleMap ~= 0, 8);
% props_all = regionprops(L_all, phzUnwrapped, 'PixelIdxList', 'PixelValues', 'Area');
% 
% % Filter out small regions
% isValid = arrayfun(@(r) r.Area > minRegionSize, props_all);
% props = props_all(isValid);
% L_filtered = zeros(size(L_all));
% for i = 1:numel(props)
%     L_filtered(props(i).PixelIdxList) = i;
% end
% regionLabels = L_filtered;
% numValid = numel(props);
% 
% % Median phase per region
% medPhase = arrayfun(@(r) median(r.PixelValues, 'omitnan'), props);
% 
% % ------------------- Step 3: Flood from All Seeds ------------------- %
% phzCorrected = phzUnwrapped;
% correctionMap = zeros(size(phzUnwrapped));
% visited = false(1, numValid);
% phzOffsets = zeros(1, numValid);
% 
% % 8-connectivity neighborhood
% connKernel = strel('square', 11).Neighborhood;
% 
% for seed = 1:numValid
%     if visited(seed), continue; end
% 
%     % Start new flood
%     queue = seed;
%     visited(seed) = true;
% 
%     while ~isempty(queue)
%         curr = queue(1); queue(1) = [];
% 
%         % Binary mask for current region
%         maskCurr = (regionLabels == curr);
% 
%         % Dilate to find neighbor pixels
%         neighborMask = imdilate(maskCurr, connKernel) & ~maskCurr & validMask;
%         neighborLabels = unique(regionLabels(neighborMask));
%         neighborLabels = neighborLabels(neighborLabels > 0 & neighborLabels <= numValid);
%         neighborLabels = neighborLabels(~visited(neighborLabels));
% 
%         for n = neighborLabels(:)'
%             % dPhi = medPhase(n) - (medPhase(curr) + 2*pi*phzOffsets(curr));
%             dPhi = medPhase(n) - medPhase(curr);
%             % cycleOffset = round(dPhi / (2*pi));
%             % cycleOffset = floor((dPhi + pi) / (2*pi));
%             % phzOffsets(n) = phzOffsets(curr) + cycleOffset;
%             cycleOffset = round((medPhase(n) - medPhase(curr)) / (2*pi));
%             phzOffsets(n) = phzOffsets(curr) + cycleOffset;
%             visited(n) = true;
%             queue(end+1) = n;
%         end
%     end
% end
% 
% % ------------------- Step 4: Apply Corrections ------------------- %
% for i = 1:numValid
%     pix = props(i).PixelIdxList;
%     phzCorrected(pix) = phzUnwrapped(pix) - 2*pi*phzOffsets(i);
%     correctionMap(pix) = phzOffsets(i);
% end
% 
% end

% function [phzCorrected, cycleMap, correctionMap, regionLabels] = correct_phase_cycles_regionwise(phzUnwrapped, coherence, minRegionSize)
% %CORRECT_PHASE_CYCLES_REGIONWISE Aligns 2pi-offset phase islands into a continuous map.
% %
% %   INPUTS:
% %       phzUnwrapped    - Unwrapped (but locally inconsistent) phase [radians]
% %       coherence       - Coherence map to mask invalid pixels
% %       minRegionSize   - Minimum size (pixels) for valid region (default: 30)
% %
% %   OUTPUTS:
% %       phzCorrected    - Phase map with corrected regional 2π cycle alignment
% %       cycleMap        - Initial cycle class (floor(phz / 2π))
% %       correctionMap   - Map of applied cycle offsets (in integer multiples of 2π)
% %       regionLabels    - Labeled map of regions used for region-wise correction
% 
%     if nargin < 3, minRegionSize = 30; end
%     phzCorrected = phzUnwrapped;
% 
%     % Mask invalid pixels
%     validMask = isfinite(phzUnwrapped) & coherence > 0.05;
% 
%     % Initial cycle map (for diagnostics)
%     cycleMap = floor(phzUnwrapped / (2*pi));
%     cycleMap(~validMask) = 0;
% 
%     % Label all connected regions
%     [regionLabels, numRegions] = bwlabel(validMask, 8);
%     if numRegions == 0
%         correctionMap = zeros(size(phzUnwrapped));
%         return
%     end
% 
%     % Compute region properties
%     props = regionprops(regionLabels, phzUnwrapped, 'PixelIdxList', 'PixelValues', 'Area');
% 
%     % Filter out small regions
%     validRegion = arrayfun(@(r) r.Area > minRegionSize, props);
%     props = props(validRegion);
%     numValid = numel(props);
% 
%     regionMap = zeros(size(regionLabels));
%     for i = 1:numValid
%         regionMap(props(i).PixelIdxList) = i;
%     end
%     regionLabels = regionMap;
% 
%     % Compute median phase of each region
%     medPhase = arrayfun(@(r) median(r.PixelValues, 'omitnan'), props);
% 
%     % Initialize visit state and cycle offsets
%     visited = false(1, numValid);
%     phzOffsets = zeros(1, numValid);
% 
%     % Start from largest region
%     [~, rootRegion] = max([props.Area]);
%     visited(rootRegion) = true;
%     queue = rootRegion;
% 
%     % Connectivity kernel
%     connKernel = strel('square', 25).Neighborhood;
% 
%     % Correction map (same size as input)
%     correctionMap = zeros(size(phzUnwrapped));
% 
%     % Region-wise flood traversal
%     while ~isempty(queue)
%         curr = queue(1); queue(1) = [];
% 
%         % Binary mask of current region
%         maskCurr = (regionLabels == curr);
% 
%         % Dilate to find neighboring candidate pixels
%         neighborMask = imdilate(maskCurr, connKernel) & validMask & (regionLabels ~= curr);
%         neighborLabels = unique(regionLabels(neighborMask));
%         % neighborLabels = neighborLabels(neighborLabels > 0 & ~visited(neighborLabels));
%         neighborLabels = neighborLabels(neighborLabels > 0 & neighborLabels <= numValid);
%         neighborLabels = neighborLabels(~visited(neighborLabels));
% 
%         for n = neighborLabels(:)'
%             % Estimate phase offset to current region
%             dPhi = medPhase(n) - (medPhase(curr) + 2*pi*phzOffsets(curr));
%             cycleOffset = round(dPhi / (2*pi));
%             phzOffsets(n) = phzOffsets(curr) + cycleOffset;
%             visited(n) = true;
%             queue(end+1) = n;
%         end
%     end
% 
%     % Apply correction to each region
%     for i = 1:numValid
%         idx = props(i).PixelIdxList;
%         phzCorrected(idx) = phzUnwrapped(idx) - 2*pi*phzOffsets(i);
%         correctionMap(idx) = phzOffsets(i);
%     end
% end
% 
% % function [phzCorrected, cycleMap, correctionMap] = correct_phase_cycles_regionwise(phzUnwrapped, coherence, minRegionSize)
% % %CORRECT_PHASE_CYCLES_REGIONWISE Aligns 2pi-offset phase islands using local floods
% % %
% % %   INPUTS:
% % %       phzUnwrapped    - Referenced (but locally inconsistent) unwrapped phase [rad]
% % %       coherence       - Coherence map to mask invalid pixels
% % %       minRegionSize   - Minimum size to keep a region (default: 30 pixels)
% % %
% % %   OUTPUTS:
% % %       phzCorrected    - Region-corrected phase map
% % %       cycleMap        - Initial floor(phz/2pi) map (helpful for diagnostics)
% % %       correctionMap   - Map of applied integer 2pi cycle offsets per pixel
% % 
% % if nargin < 3, minRegionSize = 30; end
% % phzCorrected = phzUnwrapped;
% % 
% % % --- Mask invalid regions
% % validMask = isfinite(phzUnwrapped) & coherence > 0.05;
% % 
% % % --- Initial cycle estimate
% % cycleMap = floor(phzUnwrapped / (2*pi));
% % cycleMap(~validMask) = 0;
% % 
% % % --- Label all regions
% % [L_all, ~] = bwlabel(cycleMap ~= 0, 8);
% % if max(L_all(:)) == 0
% %     correctionMap = zeros(size(phzUnwrapped));
% %     return;
% % end
% % 
% % % --- Identify connected components in L_all
% % connKernel = strel('square', 50).Neighborhood;
% % L_clustered = bwlabel(L_all > 0, 4);  % Clusters of connected regions (not just individual ones)
% % numClusters = max(L_clustered(:));
% % 
% % % --- Initialize correction map
% % correctionMap = zeros(size(phzUnwrapped));
% % 
% % for c = 1:numClusters
% %     mask_cluster = (L_clustered == c);
% %     [L, num] = bwlabel(mask_cluster & L_all > 0, 8);
% %     if num <= 1, continue; end
% % 
% %     props = regionprops(L, phzUnwrapped, 'PixelIdxList', 'PixelValues', 'Area');
% % 
% %     % Filter small islands
% %     validRegion = arrayfun(@(r) r.Area > minRegionSize, props);
% %     props = props(validRegion);
% %     if isempty(props), continue; end
% % 
% %     % Relabel region IDs
% %     L_filtered = zeros(size(L));
% %     for i = 1:numel(props)
% %         L_filtered(props(i).PixelIdxList) = i;
% %     end
% %     L = L_filtered;
% % 
% %     % Compute medians
% %     medPhase = arrayfun(@(r) median(r.PixelValues, 'omitnan'), props);
% %     visited = false(1, numel(props));
% %     phzOffsets = zeros(1, numel(props));
% % 
% %     % --- Start from largest region in this cluster
% %     [~, rootRegion] = max([props.Area]);
% %     visited(rootRegion) = true;
% %     queue = rootRegion;
% % 
% %     while ~isempty(queue)
% %         curr = queue(1); queue(1) = [];
% % 
% %         maskCurr = false(size(L));
% %         maskCurr(props(curr).PixelIdxList) = true;
% % 
% %         neighborMask = imdilate(maskCurr, connKernel) & ~maskCurr & mask_cluster;
% %         neighborLabels = unique(L(neighborMask));
% %         neighborLabels = neighborLabels(neighborLabels > 0 & ~visited(neighborLabels));
% % 
% %         for n = neighborLabels(:)'
% %             dPhi = medPhase(n) - (medPhase(curr) + 2*pi*phzOffsets(curr));
% %             cycleOffset = round(dPhi / (2*pi));
% %             phzOffsets(n) = phzOffsets(curr) + cycleOffset;
% %             visited(n) = true;
% %             queue(end+1) = n;
% %         end
% %     end
% % 
% %     % --- Apply local corrections
% %     for i = 1:numel(props)
% %         idx = props(i).PixelIdxList;
% %         phzCorrected(idx) = phzUnwrapped(idx) - 2*pi*phzOffsets(i);
% %         correctionMap(idx) = phzOffsets(i);
% %     end
% % end
% % end
% % 
% % % function [phzCorrected, cycleMap, correctionMap, regionLabels] = correct_phase_cycles_regionwise(phzUnwrapped, coherence, minRegionSize)
% % % %CORRECT_PHASE_CYCLES_REGIONWISE Aligns 2pi-offset phase islands into a continuous map.
% % % %
% % % %   INPUTS:
% % % %       phzUnwrapped   - Referenced (but locally inconsistent) unwrapped phase [radians]
% % % %       coherence      - Coherence map to mask invalid pixels
% % % %       minRegionSize  - Minimum region size to consider (default: 30 pixels)
% % % %
% % % %   OUTPUTS:
% % % %       phzCorrected   - Phase map with region-wise 2π cycle alignment applied
% % % %       cycleMap       - Integer cycles per pixel (floor(phz/2π))
% % % %       correctionMap  - Pixel-wise 2π offset correction map (in integer units)
% % % %       regionLabels   - Labeled regions used during correction
% % % 
% % % if nargin < 3, minRegionSize = 30; end
% % % phzCorrected = phzUnwrapped;
% % % 
% % % % Mask invalid
% % % validMask = isfinite(phzUnwrapped) & coherence > 0.05;
% % % 
% % % % Integer cycle map
% % % cycleMap = floor(phzUnwrapped / (2*pi));
% % % cycleMap(~validMask) = 0;
% % % 
% % % % Initial region labeling
% % % [L, num] = bwlabel(cycleMap ~= 0, 8);
% % % if num == 0
% % %     correctionMap = zeros(size(phzUnwrapped));
% % %     regionLabels = zeros(size(phzUnwrapped));
% % %     return;
% % % end
% % % 
% % % % Region properties
% % % props = regionprops(L, phzUnwrapped, 'PixelIdxList', 'PixelValues', 'Area');
% % % validRegion = arrayfun(@(r) r.Area > minRegionSize, props);
% % % 
% % % % Remap valid regions
% % % remap = zeros(1, num); % Old → new
% % % newIdx = 0;
% % % L_filtered = zeros(size(L));
% % % for i = 1:num
% % %     if validRegion(i)
% % %         newIdx = newIdx + 1;
% % %         remap(i) = newIdx;
% % %         L_filtered(props(i).PixelIdxList) = newIdx;
% % %     end
% % % end
% % % 
% % % regionLabels = L_filtered;
% % % numValid = newIdx;
% % % 
% % % % Update props to only valid
% % % props = props(validRegion);
% % % medPhase = arrayfun(@(r) median(r.PixelValues, 'omitnan'), props);
% % % 
% % % % Initialize correction tracking
% % % visited = false(1, numValid);
% % % phzOffsets = zeros(1, numValid);
% % % 
% % % % Start from largest region
% % % [~, rootRegion] = max([props.Area]);
% % % visited(rootRegion) = true;
% % % queue = rootRegion;
% % % 
% % % % Use connectivity kernel
% % % connKernel = strel('square',11).Neighborhood;
% % % 
% % % % BFS-style traversal to align regions
% % % while ~isempty(queue)
% % %     curr = queue(1); queue(1) = [];
% % % 
% % %     % Binary mask of current region
% % %     maskCurr = regionLabels == curr;
% % % 
% % %     % Neighboring pixels
% % %     neighborMask = imdilate(maskCurr, connKernel) & ~maskCurr & validMask;
% % %     neighborLabels = unique(regionLabels(neighborMask));
% % %     % neighborLabels = neighborLabels(neighborLabels > 0 & ~visited(neighborLabels));
% % %     % Ensure neighbor labels are in valid range
% % %     neighborLabels = neighborLabels(neighborLabels > 0 & neighborLabels <= numValid);
% % %     neighborLabels = neighborLabels(~visited(neighborLabels));
% % % 
% % % 
% % %     for n = neighborLabels(:)'
% % %         % Offset relative to current region
% % %         dPhi = medPhase(n) - (medPhase(curr) + 2*pi*phzOffsets(curr));
% % %         % cycleOffset = round(dPhi / (2*pi));
% % %         cycleOffset = floor((dPhi + pi) / (2*pi));  % Shifted floor to detect offsets better
% % %         phzOffsets(n) = phzOffsets(curr) + cycleOffset;
% % %         visited(n) = true;
% % %         queue(end+1) = n;
% % %     end
% % % end
% % % 
% % % % Build correction map
% % % correctionMap = zeros(size(phzUnwrapped));
% % % for i = 1:numValid
% % %     idx = props(i).PixelIdxList;
% % %     phzCorrected(idx) = phzUnwrapped(idx) - 2*pi * phzOffsets(i);
% % %     correctionMap(idx) = phzOffsets(i);
% % % end
% % % 
% % % end
% % % 
% % % % function [phzCorrected, cycleMap, cycleCorrectionMap, regionLabels] = ...
% % % %     correct_phase_cycles_regionwise(phzUnwrapped, coherence, minRegionSize)
% % % % %CORRECT_PHASE_CYCLES_REGIONWISE Aligns 2pi-offset phase islands into a continuous map.
% % % % %
% % % % %   INPUTS:
% % % % %       phzUnwrapped       - Referenced (but locally inconsistent) unwrapped phase [radians]
% % % % %       coherence          - Coherence map to mask invalid pixels
% % % % %       minRegionSize      - Minimum size to keep a region (default: 30 pixels)
% % % % %
% % % % %   OUTPUTS:
% % % % %       phzCorrected       - Phase map with region-wise 2π cycle alignment applied
% % % % %       cycleMap           - floor(phz / 2π) phase cycles before correction
% % % % %       cycleCorrectionMap - Number of 2π cycles removed per pixel
% % % % %       regionLabels       - Connected region ID (0 = background)
% % % % 
% % % %     if nargin < 3, minRegionSize = 30; end
% % % %     phzCorrected = phzUnwrapped;
% % % % 
% % % %     % Mask NaNs and low coherence
% % % %     validMask = isfinite(phzUnwrapped) & coherence > 0.05;
% % % % 
% % % %     % Build phase-cycle map
% % % %     cycleMap = floor(phzUnwrapped / (2*pi));
% % % %     cycleMap(~validMask) = 0;
% % % % 
% % % %     % Label regions
% % % %     [L, num] = bwlabel(cycleMap ~= 0, 8);
% % % %     if num == 0
% % % %         regionLabels = zeros(size(phzUnwrapped));
% % % %         cycleCorrectionMap = zeros(size(phzUnwrapped));
% % % %         return;
% % % %     end
% % % % 
% % % %     % Compute region properties
% % % %     props = regionprops(L, phzUnwrapped, 'PixelIdxList', 'PixelValues', 'Area');
% % % % 
% % % %     % Filter small junk
% % % %     validRegion = arrayfun(@(r) r.Area > minRegionSize, props);
% % % %     props = props(validRegion);
% % % %     L_filtered = zeros(size(L));
% % % %     for i = 1:numel(props)
% % % %         L_filtered(props(i).PixelIdxList) = i;
% % % %     end
% % % %     L = L_filtered;
% % % %     regionLabels = L;  % <-- Return region label map
% % % % 
% % % %     % Median phase of each region
% % % %     medPhase = arrayfun(@(r) median(r.PixelValues, 'omitnan'), props);
% % % % 
% % % %     % Track visited + cycle offset
% % % %     visited = false(1, numel(props));
% % % %     phzOffsets = zeros(1, numel(props));  % in # of 2π cycles
% % % % 
% % % %     % Start from largest region
% % % %     [~, rootRegion] = max([props.Area]);
% % % %     visited(rootRegion) = true;
% % % %     queue = rootRegion;
% % % % 
% % % %     % Neighbor offsets
% % % %     connKernel = strel('square', 3).Neighborhood;
% % % % 
% % % %     while ~isempty(queue)
% % % %         curr = queue(1); queue(1) = [];
% % % % 
% % % %         % Binary mask of current region
% % % %         maskCurr = false(size(L));
% % % %         maskCurr(props(curr).PixelIdxList) = true;
% % % % 
% % % %         % Dilate to find neighbors
% % % %         neighborMask = imdilate(maskCurr, connKernel) & ~maskCurr & validMask;
% % % %         neighborLabels = unique(L(neighborMask));
% % % %         % neighborLabels = neighborLabels(neighborLabels > 0 & ~visited(neighborLabels));
% % % %         % Ensure neighbor labels are in valid range
% % % %         neighborLabels = neighborLabels(neighborLabels > 0 & neighborLabels <= numValid);
% % % %         neighborLabels = neighborLabels(~visited(neighborLabels));
% % % % 
% % % % 
% % % %         for n = neighborLabels(:)'
% % % %             % Estimate cycle offset using region medians
% % % %             dPhi = medPhase(n) - (medPhase(curr) + 2*pi*phzOffsets(curr));
% % % %             cycleOffset = round(dPhi / (2*pi));
% % % %             phzOffsets(n) = phzOffsets(curr) + cycleOffset;
% % % %             visited(n) = true;
% % % %             queue(end+1) = n;
% % % %         end
% % % %     end
% % % % 
% % % %     % Build full correction map (pixelwise)
% % % %     cycleCorrectionMap = zeros(size(phzUnwrapped));
% % % %     for i = 1:numel(props)
% % % %         cycleCorrectionMap(props(i).PixelIdxList) = phzOffsets(i);
% % % %         phzCorrected(props(i).PixelIdxList) = ...
% % % %             phzUnwrapped(props(i).PixelIdxList) - 2*pi * phzOffsets(i);
% % % %     end
% % % % end
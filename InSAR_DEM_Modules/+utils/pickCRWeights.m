function wEntry = pickCRWeights(crC, dirIdx, burstIdx, policy)
%PICKCRWEIGHTS Select a CR weight entry for (directory, burst) with layout-compat.
%   wEntry has fields .indices (linear DEM idx) and .value (weights).
%   policy: 'burst' (use burstIdx), 'best' (use BestSLCPerDir(dirIdx)),
%           'mean' (combine all SLCs in dir, power-weighted if available).
%
% Returns [] if nothing valid is found.

if nargin < 4 || isempty(policy), policy = 'burst'; end
wEntry = [];

% NEW layout: Weights is 1xNumDirs cell; Weights{dir} is 1xNSLC cell of structs
if iscell(crC.Weights)
    if dirIdx > numel(crC.Weights) || isempty(crC.Weights{dirIdx}), return; end

    switch lower(policy)
        case 'burst'
            jSel = burstIdx;

        case 'best'
            jSel = NaN;
            if isfield(crC,'BestSLCPerDir') && dirIdx <= numel(crC.BestSLCPerDir)
                jSel = crC.BestSLCPerDir(dirIdx);
            end
            if ~isfinite(jSel) || jSel < 1 || jSel > numel(crC.Weights{dirIdx})
                jSel = burstIdx; % fallback
            end

        case 'mean'
            % Combine all SLCs' weights in this dir (power-weighted if available)
            wCells = crC.Weights{dirIdx};
            pVec   = [];
            if isfield(crC,'Power') && dirIdx <= numel(crC.Power)
                pVec = crC.Power{dirIdx};
            end
            [idxUnion, valAccum, denom] = deal([], [], 0);
            for jj = 1:numel(wCells)
                ent = wCells{jj};
                if ~isstruct(ent) || ~isfield(ent,'indices') || ~isfield(ent,'value'), continue; end
                idxs = ent.indices(:); vals = ent.value(:);
                if isempty(idxs) || isempty(vals), continue; end

                % Build union of indices, accumulate (optionally power-weighted)
                if isempty(idxUnion)
                    idxUnion = idxs; valAccum = zeros(size(idxs));
                end
                [inA, ~] = ismember(idxs, idxUnion);
                newIdxs  = idxs(~inA);
                idxUnion = [idxUnion; newIdxs];
                valAccum = [valAccum; zeros(size(newIdxs))];

                [~, locB] = ismember(idxs, idxUnion);
                pw = 1;
                if ~isempty(pVec) && jj <= numel(pVec) && isfinite(pVec(jj)), pw = pVec(jj); end
                valAccum(locB) = valAccum(locB) + pw * vals;
                denom = denom + pw;
            end
            if ~isempty(idxUnion) && denom > 0
                wEntry.indices = idxUnion;
                wEntry.value   = valAccum / denom;
            end
            return

        otherwise
            jSel = burstIdx;
    end

    if jSel >= 1 && jSel <= numel(crC.Weights{dirIdx})
        ent = crC.Weights{dirIdx}{jSel};
        if isstruct(ent) && isfield(ent,'indices') && isfield(ent,'value')
            wEntry = ent;
        end
    end

else
    % LEGACY layout: Weights is a 2-D cell array {dir, slc}
    if size(crC.Weights,1) >= dirIdx && size(crC.Weights,2) >= burstIdx
        ent = crC.Weights{dirIdx, burstIdx};
        if isstruct(ent) && isfield(ent,'indices') && isfield(ent,'value')
            wEntry = ent;
        end
    end
end
end
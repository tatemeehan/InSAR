function insarData = reference_phase(insarData, crResult, mode)
%REFERENCE_PHASE Applies phase referencing to unwrapped phase in insarData
%   insarData = reference_phase(insarData, crResult, mode)
%   mode: 'powerWeighted' (default), 'unweighted', 'sceneMean', 'none'

if nargin < 3 || isempty(mode), mode = 'powerWeighted'; end

for k = 1:numel(insarData)
    phzUnw = insarData(k).unwrapped;
    cor = insarData(k).coherence;
    i = insarData(k).pair(1);
    j = insarData(k).pair(2);

    switch lower(mode)
        case 'powerweighted'
            refPhases = [];
            weights = [];
            for c = 1:numel(crResult)
                if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
                wEntry = crResult(c).Weights{i,j};
                if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
                px = wEntry.indices;
                wt = wEntry.value;
                if any(isnan(phzUnw(px))), continue; end
                refPhases(end+1) = sum(wt .* phzUnw(px));
                weights(end+1) = sum(wt);
            end
            if ~isempty(refPhases)
                referencePhase = sum(refPhases .* weights) / sum(weights);
            else
                referencePhase = mean(phzUnw(cor >= 0.75), 'omitnan');
            end

        case 'unweighted'
            refPx = [];
            for c = 1:numel(crResult)
                if size(crResult(c).Weights,1) < i || size(crResult(c).Weights,2) < j, continue; end
                wEntry = crResult(c).Weights{i,j};
                if ~isfield(wEntry, 'indices'), continue; end
                refPx = [refPx; wEntry.indices(:)];
            end
            refPx = unique(refPx);
            referencePhase = mean(phzUnw(refPx), 'omitnan');

        case 'scenemean'
            referencePhase = mean(phzUnw(cor >= 0.75), 'omitnan');

        case 'none'
            referencePhase = 0;

        otherwise
            error('Unknown crMode: %s', mode);
    end

    insarData(k).unwrapped = phzUnw - referencePhase;
    insarData(k).refPhase = referencePhase;
end
end

% function insarData = reference_phase(insarData, crResult, mode)
% %REFERENCE_PHASE Applies phase referencing to unwrapped phase in insarData
% %   insarData = reference_phase(insarData, crResult, mode)
% %   mode: 'powerWeighted' (default), 'unweighted', 'sceneMean', 'none'
% 
% if nargin < 3 || isempty(mode), mode = 'powerWeighted'; end
% 
% for k = 1:numel(insarData)
%     phzUnw = insarData(k).unwrapped;
%     cor = insarData(k).coherence;
%     i = insarData(k).pair(1);
%     j = insarData(k).pair(2);
% 
%     switch lower(mode)
%         case 'powerweighted'
%             refPhases = [];
%             weights = [];
%             for c = 1:numel(crResult)
%                 if i > numel(crResult(c).Weights), continue; end
%                 wList = crResult(c).Weights{i};
%                 if isempty(wList) || j > numel(wList), continue; end
%                 wEntry = wList(j);
%                 if ~isstruct(wEntry) || ~isfield(wEntry, 'indices') || ~isfield(wEntry, 'value'), continue; end
%                 px = wEntry.indices;
%                 wt = wEntry.value;
%                 if any(isnan(phzUnw(px))), continue; end
%                 refPhases(end+1) = sum(wt .* phzUnw(px));
%                 weights(end+1) = sum(wt);
%             end
%             if ~isempty(refPhases)
%                 referencePhase = sum(refPhases .* weights) / sum(weights);
%             else
%                 referencePhase = mean(phzUnw(cor >= 0.75), 'omitnan');
%             end
% 
%         case 'unweighted'
%             refPx = [];
%             for c = 1:numel(crResult)
%                 if i > numel(crResult(c).Weights), continue; end
%                 wList = crResult(c).Weights{i};
%                 if isempty(wList) || j > numel(wList), continue; end
%                 wEntry = wList(j);
%                 if ~isfield(wEntry, 'indices'), continue; end
%                 refPx = [refPx; wEntry.indices(:)];
%             end
%             refPx = unique(refPx);
%             referencePhase = mean(phzUnw(refPx), 'omitnan');
% 
%         case 'scenemean'
%             referencePhase = mean(phzUnw(cor >= 0.75), 'omitnan');
% 
%         case 'none'
%             referencePhase = 0;
% 
%         otherwise
%             error('Unknown crMode: %s', mode);
%     end
% 
%     insarData(k).unwrapped = phzUnw - referencePhase;
%     insarData(k).refPhase = referencePhase;
% end
% end
% 
% 
% % function insarData = reference_phase(insarData, crResult, mode)
% % %REFERENCE_PHASE Applies phase referencing to unwrapped phase in insarData
% % %   insarData = reference_phase(insarData, crResult, mode, sarData)
% % %   mode: 'powerWeighted' (default), 'unweighted', 'sceneMean', 'none'
% % 
% % if nargin < 3 || isempty(mode), mode = 'powerWeighted'; end
% % 
% % for k = 1:numel(insarData)
% %     phzUnw = insarData(k).unwrapped;
% %     cor = insarData(k).coherence;
% %     i = insarData(k).pair(1);
% %     j = insarData(k).pair(2);
% % 
% %     switch lower(mode)
% %         case 'powerweighted'
% %             refPhases = [];
% %             weights = [];
% %             for c = 1:numel(crResult)
% %                 if i > numel(crResult(c).Weights), continue; end
% %                 w = crResult(c).Weights{i};
% %                 if isempty(w) || j > numel(w) || ~isfield(w(j), 'indices'), continue; end
% %                 px = w(j).indices;
% %                 wt = w(j).value;
% %                 if any(isnan(phzUnw(px))), continue; end
% %                 refPhases(end+1) = sum(wt .* phzUnw(px));
% %                 weights(end+1) = sum(wt);
% %             end
% %             if ~isempty(refPhases)
% %                 referencePhase = sum(refPhases .* weights) / sum(weights);
% %             else
% %                 referencePhase = mean(phzUnw(cor >= 0.75), 'omitnan');
% %             end
% % 
% %         case 'unweighted'
% %             refPx = [];
% %             for c = 1:numel(crResult)
% %                 if i > numel(crResult(c).Weights), continue; end
% %                 w = crResult(c).Weights{i};
% %                 if isempty(w) || j > numel(w) || ~isfield(w(j), 'indices'), continue; end
% %                 refPx = [refPx; w(j).indices(:)];
% %             end
% %             refPx = unique(refPx);
% %             referencePhase = mean(phzUnw(refPx), 'omitnan');
% % 
% %         case 'scenemean'
% %             referencePhase = mean(phzUnw(cor >= 0.75), 'omitnan');
% % 
% %         case 'none'
% %             referencePhase = 0;
% % 
% %         otherwise
% %             error('Unknown crMode: %s', mode);
% %     end
% % 
% %     insarData(k).unwrapped = phzUnw - referencePhase;
% %     insarData(k).refPhase = referencePhase;
% % end
% % end

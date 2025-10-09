function pairList = generatePairList(sarData, modeOrOpts, customPairs)
%GENERATEPAIRLIST Generate [i1 b1 i2 b2] pairs using build_pairs semantics.
%
% Usage:
%   pairList = generatePairList(sarData)                     % defaults to 'all'
%   pairList = generatePairList(sarData, 'all')
%   pairList = generatePairList(sarData, 'sequential')
%   pairList = generatePairList(sarData, 'custom', customPairsNx4)
%   pairList = generatePairList(sarData, optsStructLikeBuildPairs)
%
% The returned pairList is Px4: [i1 j1 i2 j2], canonicalized.

    if nargin < 2 || isempty(modeOrOpts), modeOrOpts = 'all'; end

    % --- CUSTOM explicit list (Nx4) ---
    if ischar(modeOrOpts) || isstring(modeOrOpts)
        mode = lower(string(modeOrOpts));
        switch mode
            case "custom"
                if nargin < 3 || isempty(customPairs)
                    error('Custom mode requires customPairs (Nx4).');
                end
                validateattributes(customPairs, {'numeric'}, ...
                    {'ncols',4,'finite','integer','>=',1}, mfilename, 'customPairs');
                pairList = utils.canonicalizePairList(customPairs);
                return
        end
    end

    % --- Map modes/opts -> build_pairs opts ---
    if isstruct(modeOrOpts)
        % Pass-through: accept same fields as build_pairs (includeInter, includeIntra, etc.)
        opts = modeOrOpts;
    else
        % Back-compat modes
        switch mode
            case "all"
                % All inter- and intra-dir combinations
                opts = struct( ...
                    'includeInter',true, 'includeIntra',true, ...
                    'burstModeInter','allCombos', 'burstModeIntra','allCombos', ...
                    'ensureIntraAllCombos',false);  % already all combos
            case "sequential"
                % Inter: same-index; Intra: neighbors only (classic)
                opts = struct( ...
                    'includeInter',true, 'includeIntra',true, ...
                    'burstModeInter','sameIndex', 'burstModeIntra','neighbors', ...
                    'ensureIntraAllCombos',false);
            case "inter"
                % Only inter-directory pairs, same-index bursts
                opts = struct( ...
                    'includeInter',true, 'includeIntra',false, ...
                    'burstModeInter','sameIndex');
            case "intra"
                % Only intra-directory pairs, all burst combos (b1<b2)
                opts = struct( ...
                    'includeInter',false, 'includeIntra',true, ...
                    'burstModeIntra','allCombos', 'ensureIntraAllCombos',false);
            otherwise
                error('Unknown mode: %s. Use ''all'',''sequential'',''inter'',''intra'',''custom''.', mode);
        end
    end

    % Sensible defaults if they weren't provided
    if ~isfield(opts,'minBurstsPerDir'),      opts.minBurstsPerDir = 1; end
    if ~isfield(opts,'ensureIntraAllCombos'), opts.ensureIntraAllCombos = true; end
    % (filterFn, burstModeInter/Intra can be passed through as in build_pairs)

    % --- Delegate to build_pairs, then convert to Px4 ---
    P = utils.build_pairs(sarData, opts);
    if isempty(P)
        pairList = zeros(0,4);
    else
        pairList = [[P.i1]' [P.b1]' [P.i2]' [P.b2]'];
    end

    % Canonicalize (ensures (i1,b1) < (i2,b2) ordering if your helper does that)
    pairList = utils.canonicalizePairList(pairList);
end

% Legacy Version
% function pairList = generatePairList(sarData, mode, customPairs)
% %GENERATE_CROSS_FRAME_PAIR_LIST Generate list of (ii,jj) vs (ii,jj) trajectory pairs
% %
% %   pairList = generate_cross_frame_pair_list(sarData, mode)
% %   pairList = generate_cross_frame_pair_list(..., 'custom', customPairs)
% %
% %   Inputs:
% %       sarData     - Struct array of SAR directories with traj{jj}
% %       mode        - 'all', 'sequential', or 'custom'
% %       customPairs - Optional Mx4 array: [ii1 jj1 ii2 jj2]
% %
% %   Output:
% %       pairList    - Px4 array of index pairs: [ii1 jj1 ii2 jj2]
% 
% pairIdx = [];
% for ii = 1:numel(sarData)
%     for jj = 1:numel(sarData(ii).traj)
%         pairIdx(end+1,:) = [ii, jj]; %#ok<AGROW>
%     end
% end
% 
% switch lower(mode)
%     case 'all'
%         combs = nchoosek(1:size(pairIdx,1), 2);
%         pairList = [pairIdx(combs(:,1),:)  pairIdx(combs(:,2),:)];
% 
%     case 'sequential'
%         pairList = [];
%         for k = 1:size(pairIdx,1)-1
%             pairList(end+1,:) = [pairIdx(k,:) pairIdx(k+1,:)]; %#ok<AGROW>
%         end
% 
%     case 'custom'
%         if nargin < 3 || isempty(customPairs)
%             error('Custom mode requires a customPairs input (Nx4 array).');
%         end
%         if size(customPairs,2) ~= 4
%             error('customPairs must be an Nx4 array.');
%         end
%         pairList = customPairs;
% 
%     otherwise
%         error('Unknown mode: %s. Choose ''all'', ''sequential'', or ''custom''.', mode);
% end
% % Enure Sorted Order of Pairs
% pairList = utils.canonicalizePairList(pairList);
% end
function geomData = computeSARgeometry(sarData, demData, opts)
% Flexible geometry builder with 3 modes: 'all' | 'selected' | 'trajtol'
% Outputs:
%   geomData.sarGeometry   : cell of unique geoms (Bperp, slant, inc, mask, ...)
%   geomData.pairResults   : struct per requested oriented pair with .geomIndex
%   geomData.pairList      : [i1 b1 i2 b2] oriented pairs actually wired
%   geomData.pair2geom     : containers.Map 'i1_b1_i2_b2' -> idx
%   geomData.node2geom     : cell{ii,jj} list of all geoms touching that SLC
%   geomData.slcGeomIndex  : legacy last-index per (ii,jj) w/ .pass

% ---- opts defaults
if nargin < 3, opts = struct; end
opts = defaults(opts, struct( ...
    'mode', 'all', ...            % 'all' | 'selected' | 'trajtol'
    'pairs', [], ...              % Nx4 [i1 b1 i2 b2] when mode='selected'
    'trajTol', 0.5, ...           % meters (used in 'trajtol' mode)
    'geomMatchTol', 1.0, ...      % meters (re-use threshold for isSameGeometry)
    'look', 'right' ...           % 'right' or 'left' for your engine
));

% ---- prep DEM
Xi = demData.X(:);  Yi = demData.Y(:);
DEM = demData.dem;
normalVec = reshape(demData.surfaceNormal, [], 3);
maxBursts = max(cellfun(@numel,{sarData.slc}));

% ---- choose the pair list by mode
switch lower(opts.mode)
    case 'all'
        reqPairs = utils.generatePairList(sarData, 'all');   % oriented [i1 b1 i2 b2]
    case 'selected'
        validatePairs(opts.pairs, sarData);
        reqPairs = opts.pairs;
    case 'trajtol'
        % Only make pairs whose two trajectories are within 'trajTol'
        reqPairs = pairsByTrajectoryTolerance(sarData, opts.trajTol);
    otherwise
        error('computeSARgeometry: unknown opts.mode=%s', opts.mode);
end

% ---- storages
geomStore      = {};                                      % unique geoms
pairResults    = repmat(struct('ii1',[], 'jj1',[], 'ii2',[], 'jj2',[], 'geomIndex',[]), 0, 1);
pair2geom      = containers.Map('KeyType','char','ValueType','int32'); % oriented
geomKeyMap     = containers.Map();                        % undirected reuse
node2geom      = cell(numel(sarData), maxBursts);         % per (ii,jj) list
slcGeomIndex   = repmat(struct('idx', NaN, 'pass', ''), numel(sarData), maxBursts);

% ---- main loop
for p = 1:size(reqPairs,1)
    ii1 = reqPairs(p,1); jj1 = reqPairs(p,2);
    ii2 = reqPairs(p,3); jj2 = reqPairs(p,4);

    r1 = sarData(ii1).traj{jj1};
    r2 = sarData(ii2).traj{jj2};

    % undirected reuse key => same geom for [i1,b1;i2,b2] and [i2,b2;i1,b1]
    undirectedKey = undirectedPairKey(ii1,jj1,ii2,jj2);

    if isKey(geomKeyMap, undirectedKey)
        geomIdx = geomKeyMap(undirectedKey);
    else
        % try to match an existing unique geometry by trajectory similarity
        matchIdx = [];
        for k = 1:numel(geomStore)
            if utils.isSameGeometry(r1, r2, geomStore{k}.r1, geomStore{k}.r2, opts.geomMatchTol)
                matchIdx = k; break
            end
        end

        if isempty(matchIdx)
            [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
                Xi, Yi, DEM, r1, r2, normalVec, opts.look);

            geomIdx = numel(geomStore) + 1;
            geomStore{geomIdx} = struct( ...
                'Bperp', Bperp, ...
                'slant', rslant, ...
                'incidence', inc, ...
                'lookMask', mask, ...
                'slant2', r2slant, ...
                'incidence2', r2inc, ...
                'r1', r1, 'r2', r2, ...
                'pair', [ii1 jj1 ii2 jj2]);
        else
            geomIdx = matchIdx;
        end
        geomKeyMap(undirectedKey) = geomIdx;
    end

    % oriented keys → same geom index
    kFwd = sprintf('%d_%d_%d_%d', ii1,jj1, ii2,jj2);
    kRev = sprintf('%d_%d_%d_%d', ii2,jj2, ii1,jj1);
    pair2geom(kFwd) = int32(geomIdx);
    pair2geom(kRev) = int32(geomIdx);

    % legacy per-node last index + role
    slcGeomIndex(ii1,jj1).idx  = geomIdx; slcGeomIndex(ii1,jj1).pass = 'master';
    slcGeomIndex(ii2,jj2).idx  = geomIdx; slcGeomIndex(ii2,jj2).pass = 'slave';

    % full list per node (for diagnostics or alt lookups)
    node2geom{ii1,jj1}(end+1) = struct('idx',geomIdx,'role','master','mate',[ii2 jj2]);
    node2geom{ii2,jj2}(end+1) = struct('idx',geomIdx,'role','slave', 'mate',[ii1 jj1]);

    % record
    pairResults(end+1) = struct('ii1',ii1,'jj1',jj1,'ii2',ii2,'jj2',jj2,'geomIndex',geomIdx);
end

% ---- pack
geomData.sarGeometry  = geomStore;
geomData.pairResults  = pairResults;
geomData.pairList     = reqPairs;
geomData.slcGeomIndex = slcGeomIndex;
geomData.pair2geom    = pair2geom;
geomData.node2geom    = node2geom;
end

% ---------- helpers ----------
function S = defaults(S, D)
f = fieldnames(D);
for k=1:numel(f), if ~isfield(S,f{k}) || isempty(S.(f{k})), S.(f{k}) = D.(f{k}); end, end
end

function key = undirectedPairKey(i1,b1,i2,b2)
kk = sortrows([i1 b1; i2 b2]);
key = sprintf('%d_%d_%d_%d', kk.');
end

function validatePairs(P, sarData)
if isempty(P) || size(P,2)~=4
    error('opts.pairs must be Nx4 [i1 b1 i2 b2].');
end
if any(P(:) < 1), error('Pair indices must be positive.'); end
imax = numel(sarData);
for r=1:size(P,1)
    if P(r,1)>imax || P(r,3)>imax
        error('Pair row %d references SLC index out of range.', r);
    end
    if P(r,2)>numel(sarData(P(r,1)).slc) || P(r,4)>numel(sarData(P(r,3)).slc)
        error('Pair row %d references burst out of range.', r);
    end
end
end

function reqPairs = pairsByTrajectoryTolerance(sarData, tolMeters)
% Build oriented pairs only when the two SLC trajectories are within tol
idx = [];  % [ii jj] rows
for ii = 1:numel(sarData)
    for jj = 1:numel(sarData(ii).slc)
        idx(end+1,:) = [ii jj]; %#ok<AGROW>
    end
end
% precompute simple trajectory signatures (center + mean altitude)
cent = zeros(size(idx,1),3);
for k = 1:size(idx,1)
    ii = idx(k,1); jj = idx(k,2);
    r  = sarData(ii).traj{jj};
    cent(k,:) = mean(r,1);
end

reqPairs = zeros(0,4);
for a = 1:size(idx,1)
    for b = 1:size(idx,1)
        if a==b, continue; end
        da = norm(cent(a,:)-cent(b,:));
        if da <= tolMeters
            reqPairs(end+1,:) = [idx(a,:) idx(b,:)]; %#ok<AGROW>
        end
    end
end
end

% function [geomData] = computeSARgeometry(sarData, demData, tol)
% %COMPUTESARGEOMETRY Compute and cache SAR interferometric geometry.
% %
% % This function computes geometric InSAR parameters (e.g., perpendicular
% % baseline, slant range, incidence angle, and visibility mask) for all
% % valid trajectory pairs across one or more SAR frames.
% %
% % It reuses geometry computations when repeat-pass (identical) trajectories
% % are encountered, avoiding redundant calculations.
% %
% % INPUTS:
% %   sarData   : Structure array containing SAR trajectory data. Each
% %               sarData(ii).traj{jj} must be a Tx3 matrix of positions.
% %   demData   : Structure containing DEM and surface normal info:
% %                  - demData.dem            : MxN elevation grid
% %                  - demData.X, demData.Y   : MxN meshgrid coordinates
% %                  - demData.surfaceNormal  : MxN x 3 surface normal vectors
% %   tol       : Tolerance for detecting matching geometries (default: 1.0 meter)
% %
% % OUTPUT:
% %   geomData.sarGeometry : Cell array of unique geometry results (cached)
% %   geomData.pairResults : Struct array mapping pairs to geometry indices
% %   geomData.pairList    : Raw list of (ii1, jj1, ii2, jj2) pair indices
% %   geomData.slcGeomIndex: Matrix of geometry index and pass type per SLC
% 
% if nargin < 3, tol = 1.0; end
% 
% pairList = utils.generatePairList(sarData, 'intra');   % [i1 b1 i2 b2] (oriented)
% [m, n] = size(demData.dem);
% Xi = demData.X(:);  Yi = demData.Y(:);
% DEM = demData.dem;
% normalVec = reshape(demData.surfaceNormal, [], 3);
% 
% geomStore   = {};                                        % cell of unique geoms
% pairResults = struct([]);                                % per-row mapping
% geomKeyMap  = containers.Map();                          % undirected reuse key
% maxBursts   = max(cellfun(@numel,{sarData.slc}));
% slcGeomIndex = repmat(struct('idx', NaN, 'pass', ''), numel(sarData), maxBursts);
% 
% % pair->geometry direct lookup (oriented), node->list of geoms
% pair2geom = containers.Map('KeyType','char','ValueType','int32');
% node2geom = cell(numel(sarData), maxBursts);             % each is [] or struct array
% 
% for p = 1:size(pairList, 1)
%     ii1 = pairList(p,1); jj1 = pairList(p,2);
%     ii2 = pairList(p,3); jj2 = pairList(p,4);
% 
%     r1 = sarData(ii1).traj{jj1};
%     r2 = sarData(ii2).traj{jj2};
% 
%     % undirected cache key so [i1,b1;i2,b2] == [i2,b2;i1,b1]
%     keyVec    = [ii1, jj1, ii2, jj2];
%     sortedKey = sortrows([keyVec(1:2); keyVec(3:4)]);
%     undirectedKey = sprintf('%d_%d_%d_%d', sortedKey.');
% 
%     if isKey(geomKeyMap, undirectedKey)
%         geomIdx = geomKeyMap(undirectedKey);
%     else
%         matchIdx = [];
%         for k = 1:numel(geomStore)
%             ru1 = geomStore{k}.r1; ru2 = geomStore{k}.r2;
%             if utils.isSameGeometry(r1, r2, ru1, ru2, tol)
%                 matchIdx = k; break
%             end
%         end
%         if isempty(matchIdx)
%             [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
%                 Xi, Yi, DEM, r1, r2, normalVec, 'right');  % <- your function
% 
%             geomIdx = numel(geomStore) + 1;
%             geomStore{geomIdx} = struct( ...
%                 'Bperp', Bperp, ...
%                 'slant', rslant, ...        
%                 'incidence', inc, ...
%                 'lookMask', mask, ...
%                 'slant2', r2slant, ...       
%                 'incidence2', r2inc, ...
%                 'r1', r1, 'r2', r2, ...
%                 'pair', [ii1 jj1 ii2 jj2]);  
%         else
%             geomIdx = matchIdx;
%         end
%         geomKeyMap(undirectedKey) = geomIdx;
%     end
% 
%     % oriented pair lookup (both directions -> same geomIdx)
%     kFwd = sprintf('%d_%d_%d_%d', ii1,jj1, ii2,jj2);
%     kRev = sprintf('%d_%d_%d_%d', ii2,jj2, ii1,jj1);
%     pair2geom(kFwd) = int32(geomIdx);
%     pair2geom(kRev) = int32(geomIdx);
% 
%     % Keep your original per-node last-index (legacy / not unique)
%     slcGeomIndex(ii1,jj1).idx  = geomIdx; slcGeomIndex(ii1,jj1).pass = 'master';
%     slcGeomIndex(ii2,jj2).idx  = geomIdx; slcGeomIndex(ii2,jj2).pass = 'slave';
% 
%     % NEW: append full list per node (master/slave role, mate)
%     node2geom{ii1,jj1}(end+1) = struct('idx',geomIdx,'role','master','mate',[ii2 jj2]);
%     node2geom{ii2,jj2}(end+1) = struct('idx',geomIdx,'role','slave', 'mate',[ii1 jj1]);
% 
%     % Keep your pairResults
%     pairResults(p).ii1 = ii1; pairResults(p).jj1 = jj1;
%     pairResults(p).ii2 = ii2; pairResults(p).jj2 = jj2;
%     pairResults(p).geomIndex = geomIdx;
% end
% 
% % Add self-pairs for any orphan nodes (optional, as before)
% usedPairs = unique([pairList(:,1:2); pairList(:,3:4)], 'rows');
% for ii = 1:numel(sarData)
%     for jj = 1:numel(sarData(ii).slc)
%         if ismember([ii jj], usedPairs, 'rows'), continue; end
%         r = sarData(ii).traj{jj};
%         sortedKey   = sortrows([ii jj; ii jj]);
%         undirectedKey = sprintf('%d_%d_%d_%d', sortedKey.');
%         if ~isKey(geomKeyMap, undirectedKey)
%             matchIdx = [];
%             for k = 1:numel(geomStore)
%                 if utils.isSameGeometry(r, r, geomStore{k}.r1, geomStore{k}.r2, tol)
%                     matchIdx = k; break
%                 end
%             end
%             if isempty(matchIdx)
%                 [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
%                     Xi, Yi, DEM, r, r, normalVec, 'right');
%                 matchIdx = numel(geomStore) + 1;
%                 geomStore{matchIdx} = struct( ...
%                     'Bperp', Bperp, ...
%                     'slant', rslant, 'incidence', inc, 'lookMask', mask, ...
%                     'slant2', r2slant,'incidence2', r2inc, ...
%                     'r1', r, 'r2', r, 'pair',[ii jj ii jj]);
%             end
%             geomKeyMap(undirectedKey) = matchIdx;
%             pairList = [pairList; ii, jj, ii, jj];
%             pairResults(end+1).ii1 = ii; pairResults(end).jj1 = jj;
%             pairResults(end).ii2 = ii; pairResults(end).jj2 = jj;
%             pairResults(end).geomIndex = matchIdx;
% 
%             kSelf = sprintf('%d_%d_%d_%d', ii,jj, ii,jj);
%             pair2geom(kSelf) = int32(matchIdx);
%             node2geom{ii,jj}(end+1) = struct('idx',matchIdx,'role','self','mate',[ii jj]);
%         end
%         slcGeomIndex(ii,jj).idx = geomKeyMap(undirectedKey);
%         slcGeomIndex(ii,jj).pass = 'master';
%     end
% end
% 
% % Package outputs
% geomData.sarGeometry   = geomStore;
% geomData.pairResults   = pairResults;
% geomData.pairList      = pairList;
% geomData.slcGeomIndex  = slcGeomIndex;   % legacy (last-writer wins)
% geomData.pair2geom     = pair2geom;      % NEW: authoritative lookup
% geomData.node2geom     = node2geom;      % NEW: all geoms per node
% end
% 
% %% Previous Version
% % if nargin < 3, tol = 1.0; end
% % 
% % pairList = utils.generatePairList(sarData, 'all');
% % [m, n] = size(demData.dem);
% % Xi = demData.X(:);
% % Yi = demData.Y(:);
% % DEM = demData.dem;
% % normalVec = reshape(demData.surfaceNormal, [], 3);
% % 
% % geomStore = {};
% % pairResults = struct([]);
% % geomKeyMap = containers.Map();
% % slcGeomIndex = repmat(struct('idx', NaN, 'pass', ''), numel(sarData), max(cellfun(@numel,{sarData.slc})));
% % 
% % for p = 1:size(pairList, 1)
% %     ii1 = pairList(p,1); jj1 = pairList(p,2);
% %     ii2 = pairList(p,3); jj2 = pairList(p,4);
% % 
% %     r1 = sarData(ii1).traj{jj1};
% %     r2 = sarData(ii2).traj{jj2};
% % 
% %     keyVec = [ii1, jj1, ii2, jj2];
% %     sortedKey = sortrows([keyVec(1:2); keyVec(3:4)]);
% %     key = sprintf('%d_%d_%d_%d', sortedKey');
% % 
% %     if isKey(geomKeyMap, key)
% %         geomIdx = geomKeyMap(key);
% %     else
% %         matchIdx = [];
% %         for k = 1:numel(geomStore)
% %             ru1 = geomStore{k}.r1;
% %             ru2 = geomStore{k}.r2;
% %             if utils.isSameGeometry(r1, r2, ru1, ru2, tol)
% %                 matchIdx = k;
% %                 break
% %             end
% %         end
% % 
% %         if isempty(matchIdx)
% %             [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
% %                 Xi, Yi, DEM, r1, r2, normalVec, 'right');
% % 
% %             geomIdx = numel(geomStore) + 1;
% %             geomStore{geomIdx} = struct( ...
% %                 'Bperp', Bperp, ...
% %                 'slant', rslant, ...
% %                 'incidence', inc, ...
% %                 'lookMask', mask, ...
% %                 'slant2', r2slant, ...
% %                 'incidence2', r2inc, ...
% %                 'r1', r1, ...
% %                 'r2', r2);
% %         else
% %             geomIdx = matchIdx;
% %         end
% % 
% %         geomKeyMap(key) = geomIdx;
% %     end
% % 
% %     slcGeomIndex(ii1,jj1).idx = geomIdx;
% %     slcGeomIndex(ii1,jj1).pass = 'master';
% %     slcGeomIndex(ii2,jj2).idx = geomIdx;
% %     slcGeomIndex(ii2,jj2).pass = 'slave';
% % 
% %     pairResults(p).ii1 = ii1;
% %     pairResults(p).jj1 = jj1;
% %     pairResults(p).ii2 = ii2;
% %     pairResults(p).jj2 = jj2;
% %     pairResults(p).geomIndex = geomIdx;
% % end
% % 
% % usedPairs = unique([pairList(:,1:2); pairList(:,3:4)], 'rows');
% % 
% % for ii = 1:numel(sarData)
% %     for jj = 1:numel(sarData(ii).slc)
% %         if ismember([ii jj], usedPairs, 'rows')
% %             continue
% %         end
% % 
% %         r = sarData(ii).traj{jj};
% %         sortedKey = sortrows([ii jj; ii jj]);
% %         key = sprintf('%d_%d_%d_%d', sortedKey');
% % 
% %         if ~isKey(geomKeyMap, key)
% %             matchIdx = [];
% %             for k = 1:numel(geomStore)
% %                 if utils.isSameGeometry(r, r, geomStore{k}.r1, geomStore{k}.r2, tol)
% %                     matchIdx = k;
% %                     break
% %                 end
% %             end
% % 
% %             if isempty(matchIdx)
% %                 [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
% %                     Xi, Yi, DEM, r, r, normalVec, 'right');
% % 
% %                 matchIdx = numel(geomStore) + 1;
% %                 geomStore{matchIdx} = struct( ...
% %                     'Bperp', Bperp, ...
% %                     'slant', rslant, ...
% %                     'incidence', inc, ...
% %                     'lookMask', mask, ...
% %                     'slant2', r2slant, ...
% %                     'incidence2', r2inc, ...
% %                     'r1', r, ...
% %                     'r2', r);
% %             end
% % 
% %             geomKeyMap(key) = matchIdx;
% %             pairList = [pairList; ii, jj, ii, jj];
% %             pairResults(end+1).ii1 = ii;
% %             pairResults(end).jj1 = jj;
% %             pairResults(end).ii2 = ii;
% %             pairResults(end).jj2 = jj;
% %             pairResults(end).geomIndex = matchIdx;
% %         end
% % 
% %         slcGeomIndex(ii,jj).idx = geomKeyMap(key);
% %         slcGeomIndex(ii,jj).pass = 'master';
% %     end
% % end
% % 
% % geomData.sarGeometry = geomStore;
% % geomData.pairResults = pairResults;
% % geomData.pairList = pairList;
% % geomData.slcGeomIndex = slcGeomIndex;
% % end
% 
% % function [geomData] = computeSARgeometry(sarData, demData, tol)
% % %COMPUTESARGEOMETRY Compute and cache SAR interferometric geometry.
% % %
% % % This function computes geometric InSAR parameters (e.g., perpendicular
% % % baseline, slant range, incidence angle, and visibility mask) for all
% % % valid trajectory pairs across one or more SAR frames.
% % %
% % % It reuses geometry computations when repeat-pass (identical) trajectories
% % % are encountered, avoiding redundant calculations.
% % %
% % % INPUTS:
% % %   sarData   : Structure array containing SAR trajectory data. Each
% % %               sarData(ii).traj{jj} must be a Tx3 matrix of positions.
% % %   demData   : Structure containing DEM and surface normal info:
% % %                  - demData.dem            : MxN elevation grid
% % %                  - demData.X, demData.Y   : MxN meshgrid coordinates
% % %                  - demData.surfaceNormal  : MxN x 3 surface normal vectors
% % %   tol       : Tolerance for detecting matching geometries (default: 1.0 meter)
% % %
% % % OUTPUT:
% % %   geomData.sarGeometry : Cell array of unique geometry results (cached)
% % %   geomData.pairResults : Struct array mapping pairs to geometry indices
% % %   geomData.pairList    : Raw list of (ii1, jj1, ii2, jj2) pair indices
% % %   geomData.slcGeomIndex: Matrix of geometry indices for each SLC
% % %
% % % Calls:
% % %   - utils.generatePairList()
% % %   - insar.InSARgeometry2()
% % %   - utils.isSameGeometry()
% % %
% % % Written by the Surfing SARonauts ∆, 2025
% % 
% % if nargin < 3, tol = 1.0; end  % Default geometry matching tolerance
% % 
% % % Generate list of all SAR trajectory pairs
% % pairList = utils.generatePairList(sarData, 'all');
% % 
% % % Flatten DEM and coordinate grids
% % [m, n] = size(demData.dem);
% % Xi = demData.X(:);
% % Yi = demData.Y(:);
% % DEM = demData.dem;
% % normalVec = reshape(demData.surfaceNormal, [], 3);
% % 
% % % Initialize geometry cache and results
% % geomStore = {};
% % pairResults = struct([]);
% % geomKeyMap = containers.Map();
% % slcGeomIndex = nan(numel(sarData), max(cellfun(@numel,{sarData.slc})));
% % 
% % % Loop through each pair of SAR trajectories
% % for p = 1:size(pairList, 1)
% %     ii1 = pairList(p,1); jj1 = pairList(p,2);
% %     ii2 = pairList(p,3); jj2 = pairList(p,4);
% % 
% %     r1 = sarData(ii1).traj{jj1};
% %     r2 = sarData(ii2).traj{jj2};
% % 
% %     keyVec = [ii1, jj1, ii2, jj2];
% %     sortedKey = sortrows([keyVec(1:2); keyVec(3:4)]);
% %     key = sprintf('%d_%d_%d_%d', sortedKey');
% % 
% %     if isKey(geomKeyMap, key)
% %         geomIdx = geomKeyMap(key);
% %     else
% %         matchIdx = [];
% %         for k = 1:numel(geomStore)
% %             ru1 = geomStore{k}.r1;
% %             ru2 = geomStore{k}.r2;
% %             if utils.isSameGeometry(r1, r2, ru1, ru2, tol)
% %                 matchIdx = k;
% %                 break
% %             end
% %         end
% % 
% %         if isempty(matchIdx)
% %             [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
% %                 Xi, Yi, DEM, r1, r2, normalVec, 'right');
% % 
% %             geomIdx = numel(geomStore) + 1;
% %             geomStore{geomIdx} = struct( ...
% %                 'Bperp', Bperp, ...
% %                 'slant', rslant, ...
% %                 'incidence', inc, ...
% %                 'lookMask', mask, ...
% %                 'slant2', r2slant, ...
% %                 'incidence2', r2inc, ...
% %                 'r1', r1, ...
% %                 'r2', r2);
% %         else
% %             geomIdx = matchIdx;
% %         end
% % 
% %         geomKeyMap(key) = geomIdx;
% %     end
% % 
% %     slcGeomIndex(ii1,jj1) = geomIdx;
% %     slcGeomIndex(ii2,jj2) = geomIdx;
% % 
% %     pairResults(p).ii1 = ii1;
% %     pairResults(p).jj1 = jj1;
% %     pairResults(p).ii2 = ii2;
% %     pairResults(p).jj2 = jj2;
% %     pairResults(p).geomIndex = geomIdx;
% % end
% % 
% % % Get all used (ii,jj) combinations
% % usedPairs = unique([pairList(:,1:2); pairList(:,3:4)], 'rows');
% % 
% % % Add repeat-pass (self) pairs if missing
% % for ii = 1:numel(sarData)
% %     for jj = 1:numel(sarData(ii).slc)
% %         if ismember([ii jj], usedPairs, 'rows')
% %             continue
% %         end
% % 
% %         r = sarData(ii).traj{jj};
% %         sortedKey = sortrows([ii jj; ii jj]);
% %         key = sprintf('%d_%d_%d_%d', sortedKey');
% % 
% %         if ~isKey(geomKeyMap, key)
% %             matchIdx = [];
% %             for k = 1:numel(geomStore)
% %                 if utils.isSameGeometry(r, r, geomStore{k}.r1, geomStore{k}.r2, tol)
% %                     matchIdx = k;
% %                     break
% %                 end
% %             end
% % 
% %             if isempty(matchIdx)
% %                 [Bperp, rslant, inc, mask, r2slant, r2inc] = insar.InSARgeometry2( ...
% %                     Xi, Yi, DEM, r, r, normalVec, 'right');
% % 
% %                 matchIdx = numel(geomStore) + 1;
% %                 geomStore{matchIdx} = struct( ...
% %                     'Bperp', Bperp, ...
% %                     'slant', rslant, ...
% %                     'incidence', inc, ...
% %                     'lookMask', mask, ...
% %                     'slant2', r2slant, ...
% %                     'incidence2', r2inc, ...
% %                     'r1', r, ...
% %                     'r2', r);
% %             end
% % 
% %             geomKeyMap(key) = matchIdx;
% %             pairList = [pairList; ii, jj, ii, jj];
% %             pairResults(end+1).ii1 = ii;
% %             pairResults(end).jj1 = jj;
% %             pairResults(end).ii2 = ii;
% %             pairResults(end).jj2 = jj;
% %             pairResults(end).geomIndex = matchIdx;
% %         end
% % 
% %         slcGeomIndex(ii,jj) = geomKeyMap(key);
% %     end
% % end
% % 
% % geomData.sarGeometry = geomStore;
% % geomData.pairResults = pairResults;
% % geomData.pairList = pairList;
% % geomData.slcGeomIndex = slcGeomIndex;
% % end

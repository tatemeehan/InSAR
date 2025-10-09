function insarClosureData = compute_insar_closure_bias(insarData, opts)
% Closure via complex products: angle(C_uv * C_vw * C_wu)
% Orientation is handled by conjugating reverse edges.
%
% opts.fieldComplex  : name of complex field (default 'complexCoherence')
% opts.requireFiniteAll (default true)
% opts.sameBurstOnly : if true, require u,v,w all have the same burst (default false)
% opts.sortBurstsStrict (default false)

if nargin < 2, opts = struct(); end
if ~isfield(opts,'fieldComplex'),     opts.fieldComplex = 'complexPhase'; end
if ~isfield(opts,'requireFiniteAll'), opts.requireFiniteAll = true; end
if ~isfield(opts,'sameBurstOnly'),    opts.sameBurstOnly = false; end
if ~isfield(opts,'sortBurstsStrict'), opts.sortBurstsStrict = false; end

insarClosureData = struct([]);
N = numel(insarData);
if N==0, return; end

% Normalize pair/burst
pairs  = vertcat(insarData.pair);      % [N x 2]
bursts = vertcat(insarData.burst);     % [N x 1] or [N x 2]
if size(bursts,2)==1, bursts = [bursts bursts]; end

% --- Node maps: (dir,burst)->nodeId and reverse ---
nodeId = containers.Map('KeyType','char','ValueType','int32');
rev    = {}; nid = 0;
nodeKey = @(d,b) sprintf('%d_%d',d,b);

function id = add_node(d,b)
    keystr = nodeKey(d,b);           % <-- was 'k'
    if nodeId.isKey(keystr)
        id = nodeId(keystr);
    else
        nid = nid + 1;
        nodeId(keystr) = nid;
        rev{nid} = [d,b];
        id = nid;
    end
end

% Oriented complex edge map: forward stores C, reverse stores conj(C)
% Build oriented complex edge map
ek    = @(u,v) sprintf('%d_%d',u,v);
edgeC = containers.Map('KeyType','char','ValueType','any');
adj   = cell(0);

for k = 1:N
    i = pairs(k,1); j = pairs(k,2);
    b1 = bursts(k,1); b2 = bursts(k,2);
    u = add_node(i,b1); v = add_node(j,b2);

    % Pick a complex field, or synthesize from wrapped phase
    if isfield(insarData(k),'complexPhase') && ~isempty(insarData(k).complexPhase)
        C = insarData(k).complexCoherence;
    elseif isfield(insarData(k),'phzWrapped') && ~isempty(insarData(k).phzWrapped)
        C = exp(1i * insarData(k).phzWrapped);
    else
        continue; % nothing to use
    end

    edgeC(ek(u,v)) = C;           % forward
    edgeC(ek(v,u)) = conj(C);     % reverse

    m = max(u,v); if numel(adj)<m, adj{m} = []; end
    adj{u}(end+1) = v; adj{v}(end+1) = u;
end

numNodes = numel(rev);
if numNodes < 3, return; end

% Enumerate triangles u<v<w using adjacency and presence in edgeC
tri = [];
for u = 1:numel(rev)-2
    if u>numel(adj) || isempty(adj{u}), continue; end
    Nu = unique(adj{u}); Nu = Nu(Nu>u);
    for v = Nu
        if v>numel(adj) || isempty(adj{v}), continue; end
        Nv = unique(adj{v}); Nv = Nv(Nv>v);
        for w = Nv
            if edgeC.isKey(ek(u,v)) && edgeC.isKey(ek(v,w)) && edgeC.isKey(ek(w,u))
                tri(end+1,:) = [u v w]; %#ok<AGROW>
            end
        end
    end
end

out = struct([]); tcount = 0;
for t = 1:size(tri,1)
    u = tri(t,1); v = tri(t,2); w = tri(t,3);
    Cuv = edgeC(ek(u,v)); Cvw = edgeC(ek(v,w)); Cwu = edgeC(ek(w,u));
    if ~isequal(size(Cuv),size(Cvw),size(Cwu)), continue; end
    closureMap = angle(Cuv .* Cvw .* Cwu);   % (-pi,pi]

    finite = isfinite(closureMap);
    tcount = tcount+1;
    out(tcount).triplet       = [rev{u}(1) rev{v}(1) rev{w}(1)];
    out(tcount).burst         = [rev{u}(2) rev{v}(2) rev{w}(2)];
    out(tcount).meanBias      = mean(closureMap(finite),'omitnan');
    out(tcount).stdBias       = std(closureMap(finite),'omitnan');
    out(tcount).validFraction = sum(finite(:))/numel(finite);
    out(tcount).closureMap    = closureMap;
end
insarClosureData = out;


% Sort results (canonical by dir triplet; bursts secondary)
if ~isempty(insarClosureData)
    M = numel(insarClosureData);
    tripCanon = zeros(M,3);
    burstKey  = nan(M,1);
    for m = 1:M
        d = insarClosureData(m).triplet(:).';
        tripCanon(m,:) = sort(d);
        b = insarClosureData(m).burst;
        if isscalar(b), burstKey(m) = b; else, burstKey(m) = NaN; end
    end
    [~,idx] = sortrows([tripCanon, burstKey]);
    insarClosureData = insarClosureData(idx);
end
end

% function insarClosureData = compute_insar_closure_bias(insarData, opts)
% % Compute closure bias over all 3-cycles in the (dir,burst) graph.
% % Each interferogram defines a directed edge between nodes:
% %   node = (dir, burstIndex)
% % We require 3 directed edges (u->v), (v->w), (w->u).
% %
% % opts.field            (default 'phzReferenced')
% % opts.requireFiniteAll (default true)
% % opts.sortBurstsStrict (default false)
% 
% if nargin < 2, opts = struct(); end
% if ~isfield(opts,'field'),            opts.field = 'complexCoherence'; end
% if ~isfield(opts,'requireFiniteAll'), opts.requireFiniteAll = true; end
% if ~isfield(opts,'sortBurstsStrict'), opts.sortBurstsStrict = false; end
% 
% insarClosureData = struct([]);
% Nifg = numel(insarData);
% if Nifg == 0, return; end
% 
% % --- Normalize pair/burst fields into arrays ---
% pairs  = vertcat(insarData.pair);      % [Nifg x 2], [i j]
% bursts = vertcat(insarData.burst);     % [Nifg x 1] or [Nifg x 2]
% if size(bursts,2) == 1
%     % legacy: same burst on both ends
%     bursts = [bursts bursts];
% end
% 
% % --- Node maps: (dir,burst)->nodeId and reverse ---
% nodeId = containers.Map('KeyType','char','ValueType','int32');
% rev    = {}; nid = 0;
% nodeKey = @(d,b) sprintf('%d_%d',d,b);
%     function id = add_node(d,b)
%         k = nodeKey(d,b);
%         if nodeId.isKey(k)
%             id = nodeId(k);
%         else
%             nid = nid + 1;
%             nodeId(k) = nid;
%             rev{nid} = [d,b];
%             id = nid;
%         end
%     end
% 
% % --- Edge index map built directly from insarData (both orientations) ---
% ek = @(a,b) sprintf('%d_%d',a,b);
% 
% % IMPORTANT: make this a fresh map with ValueType 'double'
% edgeIdx = containers.Map('KeyType','char','ValueType','double');
% 
% for k = 1:Nifg
%     i  = pairs(k,1);  j  = pairs(k,2);
%     b1 = bursts(k,1); b2 = bursts(k,2);
% 
%     u = add_node(i,b1);
%     v = add_node(j,b2);
% 
%     edgeIdx(ek(u,v)) = +double(k);   % forward
%     edgeIdx(ek(v,u)) = -double(k);   % reverse
% end
% 
% numNodes = numel(rev);
% if numNodes < 3, return; end
% 
% % --- Enumerate triangles u<v<w and require all three directed edges ---
% triangles = [];  % [u v w]
% for u = 1:numNodes-2
%     for v = u+1:numNodes-1
%         for w = v+1:numNodes
%             has_uv = edgeIdx.isKey(sprintf('%d_%d',u,v));
%             has_vw = edgeIdx.isKey(sprintf('%d_%d',v,w));
%             has_wu = edgeIdx.isKey(sprintf('%d_%d',w,u));
%             if has_uv && has_vw && has_wu
%                 triangles(end+1,:) = [u v w]; %#ok<AGROW>
%             end
%         end
%     end
% end
% if isempty(triangles), return; end
% 
% % --- Compute closure for each triangle ---
% out = struct([]);
% tcount = 0;
% 
% for t = 1:size(triangles,1)
%     u = triangles(t,1); v = triangles(t,2); w = triangles(t,3);
% 
%     kv = edgeIdx(ek(u,v));  s_uv = sign(kv);  k_uv = abs(kv);
%     kv = edgeIdx(ek(v,w));  s_vw = sign(kv);  k_vw = abs(kv);
%     kv = edgeIdx(ek(w,u));  s_wu = sign(kv);  k_wu = abs(kv);
% 
%     % Extra safety: verify indices are valid
%     if any([k_uv k_vw k_wu] < 1) || any([k_uv k_vw k_wu] > Nifg)
%         % Skip malformed mapping quietly
%         continue;
%     end
% 
%     A = getfield_safe(insarData(k_uv), opts.field);
%     B = getfield_safe(insarData(k_vw), opts.field);
%     C = getfield_safe(insarData(k_wu), opts.field);
% 
%     if isempty(A) || isempty(B) || isempty(C)
%         closureMap = [];
%         validFraction = 0;
%         meanBias = NaN; stdBias = NaN;
%     else
%         % Resize guard (skip if sizes differ)
%         if ~isequal(size(A), size(B), size(C))
%             warning('Closure skip: map size mismatch.');
%             continue;
%         end
% 
%         if ~isreal(A)  % complex interferograms -> conjugate for reverse edges, then product
%             if s_uv < 0, A = conj(A); end
%             if s_vw < 0, B = conj(B); end
%             if s_wu < 0, C = conj(C); end
%             closureMap = angle(A .* B .* C);           % (-pi, pi]
%         else           % phase maps (radians) -> oriented sum, then wrap
%             phi = s_uv*A + s_vw*B + s_wu*C;
%             closureMap = angle(exp(1i * phi));         % (-pi, pi]
%         end
% 
%         if opts.requireFiniteAll
%             validMask = isfinite(A) & isfinite(B) & isfinite(C);
%             closureMap(~validMask) = NaN;
%         end
% 
%         finiteMask    = isfinite(closureMap);
%         validFraction = sum(finiteMask(:)) / numel(finiteMask);
%         meanBias      = mean(closureMap(finiteMask), 'omitnan');
%         stdBias       = std(closureMap(finiteMask), 'omitnan');
%     end
% 
% 
%     dirs   = [rev{u}(1), rev{v}(1), rev{w}(1)];
%     burstsThis = [rev{u}(2), rev{v}(2), rev{w}(2)];
%     if numel(unique(burstsThis)) == 1
%         burstOut = burstsThis(1);
%     else
%         burstOut = burstsThis; % mixed-burst triangle
%     end
% 
%     tcount = tcount + 1;
%     out(tcount).triplet       = dirs;
%     out(tcount).burst         = burstOut;
%     out(tcount).meanBias      = meanBias;
%     out(tcount).stdBias       = stdBias;
%     out(tcount).validFraction = validFraction;
%     out(tcount).closureMap    = closureMap;
% end
% 
% insarClosureData = out;
% 
% % --- Canonical sort by triplet (and optionally bursts) ---
% if ~isempty(insarClosureData)
%     N = numel(insarClosureData);
%     strict = isfield(opts,'sortBurstsStrict') && opts.sortBurstsStrict;
% 
%     tripletCanon = zeros(N,3);
%     if strict, burstKey = nan(N,3); else, burstKey = nan(N,1); end
% 
%     for n = 1:N
%         dirs = insarClosureData(n).triplet(:).';
%         [dirsSorted, p] = sort(dirs);
%         tripletCanon(n,:) = dirsSorted;
% 
%         b = insarClosureData(n).burst;
%         if isscalar(b)
%             if strict, burstKey(n,:) = [b b b];
%             else,     burstKey(n,1) = b;
%             end
%         else
%             b = b(:).'; b = b(p);
%             if strict
%                 m = min(3,numel(b));
%                 burstKey(n,1:m) = b(1:m);
%             else
%                 burstKey(n,1) = min(b);
%             end
%         end
%     end
% 
%     sortMat = [tripletCanon, burstKey];
%     [~, idx] = sortrows(sortMat);
%     insarClosureData = insarClosureData(idx);
% end
% 
%     function val = getfield_safe(S, fieldname)
%         if isfield(S, fieldname), val = S.(fieldname); else, val = []; end
%     end
% end
% 
% % function insarClosureData = compute_insar_closure_bias(insarData, opts)
% % % Computes closure bias for all triangles present in insarData,
% % % with support for mixed bursts across pairs.
% % 
% % if nargin < 2, opts = struct(); end
% % if ~isfield(opts,'field'), opts.field = 'phzReferenced'; end
% % if ~isfield(opts,'requireFiniteAll'), opts.requireFiniteAll = true; end
% % if ~isfield(opts,'sortBurstsStrict'), opts.sortBurstsStrict = false; end
% % 
% % insarClosureData = struct([]);
% % if isempty(insarData), return; end
% % 
% % pairs  = vertcat(insarData.pair);      % [N×2]
% % bursts = vertcat(insarData.burst);     % [N×1] or [N×2]
% % if size(bursts,2)==1, bursts = [bursts bursts]; end
% % 
% % % Node indexing: (dir,burst) -> node id
% % nodes = containers.Map('KeyType','char','ValueType','int32');
% % rev   = {}; nid = 0;
% % key = @(d,b) sprintf('%d_%d',d,b);
% %     function idx = add_node_impl(d,b)
% %         k = key(d,b);
% %         if nodes.isKey(k)
% %             idx = nodes(k);
% %         else
% %             nid = nid + 1;
% %             nodes(k) = nid;
% %             rev{nid} = [d,b];
% %             idx = nid;
% %         end
% %     end
% % 
% % % Build directed edges list E = [u v k sgn]
% % E = [];
% % for k = 1:numel(insarData)
% %     i = pairs(k,1); j = pairs(k,2);
% %     b1 = bursts(k,1); b2 = bursts(k,2);
% %     u = add_node_impl(i,b1);
% %     v = add_node_impl(j,b2);
% %     E = [E; u, v, k, +1]; %#ok<AGROW>
% %     E = [E; v, u, k, -1]; %#ok<AGROW>
% % end
% % 
% % numNodes = numel(rev);
% % if numNodes < 3, return; end
% % 
% % % Adjacency and edge lookup
% % adj = cell(numNodes,1);
% % edgeKey = containers.Map('KeyType','char','ValueType','any');
% % for r = 1:size(E,1)
% %     u = E(r,1); v = E(r,2); k = E(r,3); sgn = E(r,4);
% %     adj{u}(end+1) = v; %#ok<AGROW>
% %     edgeKey(sprintf('%d_%d',u,v)) = [k, sgn];
% % end
% % 
% % % Find triangles (undirected uniqueness via u<v<w)
% % triangles = [];
% % for u = 1:numNodes-2
% %     if isempty(adj{u}), continue; end
% %     Nu = unique(adj{u});  Nu = Nu(Nu > u);
% %     for ii = 1:numel(Nu)
% %         v = Nu(ii);
% %         if isempty(adj{v}), continue; end
% %         Nv = unique(adj{v});  Nv = Nv(Nv > v);
% %         for jj = 1:numel(Nv)
% %             w = Nv(jj);
% %             if edgeKey.isKey(sprintf('%d_%d',w,u))
% %                 triangles = [triangles; u v w]; %#ok<AGROW>
% %             end
% %         end
% %     end
% % end
% % if isempty(triangles), return; end
% % 
% % % Process each triangle
% % for t = 1:size(triangles,1)
% %     u = triangles(t,1); v = triangles(t,2); w = triangles(t,3);
% % 
% %     kv = edgeKey(sprintf('%d_%d',u,v)); k_uv = kv(1); s_uv = kv(2);
% %     kv = edgeKey(sprintf('%d_%d',v,w)); k_vw = kv(1); s_vw = kv(2);
% %     kv = edgeKey(sprintf('%d_%d',w,u)); k_wu = kv(1); s_wu = kv(2);
% % 
% %     A = getfield_safe(insarData(k_uv), opts.field);
% %     B = getfield_safe(insarData(k_vw), opts.field);
% %     C = getfield_safe(insarData(k_wu), opts.field);
% % 
% %     if isempty(A) || isempty(B) || isempty(C)
% %         closureMap = [];
% %     else
% %         phi = s_uv*A + s_vw*B + s_wu*C;
% %         closureMap = angle(exp(1i * phi));
% %         if opts.requireFiniteAll
% %             validMask = isfinite(A) & isfinite(B) & isfinite(C);
% %             closureMap(~validMask) = NaN;
% %         end
% %     end
% % 
% %     finiteMask    = isfinite(closureMap);
% %     validFraction = sum(finiteMask(:)) / numel(finiteMask);
% %     meanBias      = mean(closureMap(finiteMask), 'omitnan');
% %     stdBias       = std(closureMap(finiteMask), 'omitnan');
% % 
% %     dirs        = [rev{u}(1), rev{v}(1), rev{w}(1)];
% %     burstsThis  = [rev{u}(2), rev{v}(2), rev{w}(2)];
% %     if numel(unique(burstsThis)) == 1
% %         burstOut = burstsThis(1);
% %     else
% %         burstOut = burstsThis; % mixed-burst triangle
% %     end
% % 
% %     insarClosureData(t).triplet       = dirs;
% %     insarClosureData(t).burst         = burstOut;
% %     insarClosureData(t).meanBias      = meanBias;
% %     insarClosureData(t).stdBias       = stdBias;
% %     insarClosureData(t).validFraction = validFraction;
% %     insarClosureData(t).closureMap    = closureMap;
% % end
% % 
% % % Sort output (canonical by triplet, then bursts)
% % if ~isempty(insarClosureData)
% %     N = numel(insarClosureData);
% %     strict = isfield(opts,'sortBurstsStrict') && opts.sortBurstsStrict;
% %     tripletCanon = zeros(N,3);
% %     if strict, burstKey = nan(N,3); else, burstKey = nan(N,1); end
% % 
% %     for n = 1:N
% %         dirs = insarClosureData(n).triplet(:).';
% %         [dirsSorted, p] = sort(dirs);
% %         tripletCanon(n,:) = dirsSorted;
% % 
% %         b = insarClosureData(n).burst;
% %         if isscalar(b)
% %             if strict, burstKey(n,:) = [b b b]; else, burstKey(n,1) = b; end
% %         else
% %             b = b(:).'; b = b(p);
% %             if strict
% %                 m = min(3, numel(b));
% %                 burstKey(n,1:m) = b(1:m);
% %             else
% %                 burstKey(n,1) = min(b);
% %             end
% %         end
% %     end
% % 
% %     sortMat = [tripletCanon, burstKey];
% %     [~, idx] = sortrows(sortMat);
% %     insarClosureData = insarClosureData(idx);
% % end
% % 
% %     function val = getfield_safe(S, fieldname)
% %         if isfield(S, fieldname), val = S.(fieldname);
% %         else, val = [];
% %         end
% %     end
% % end
% % 
% % function insarClosureData = compute_insar_closure_bias(insarData, opts)
% % % Computes closure bias for all triangles present in insarData,
% % % keeping your original output struct format and sorting by triplet+burst.
% % 
% % if nargin < 2, opts = struct(); end
% % if ~isfield(opts,'field'), opts.field = 'phzReferenced'; end
% % if ~isfield(opts,'requireFiniteAll'), opts.requireFiniteAll = true; end
% % if ~isfield(opts,'sortBurstsStrict'), opts.sortBurstsStrict = false; end
% % 
% % insarClosureData = struct([]);
% % 
% % if isempty(insarData)
% %     return;
% % end
% % 
% % pairs  = vertcat(insarData.pair);    
% % bursts = vertcat(insarData.burst);   
% % if size(bursts,2)==1, bursts = [bursts bursts]; end
% % 
% % % Map from (dir,burst) to node ID
% % nodes = containers.Map('KeyType','char','ValueType','int32');
% % rev   = {}; nid = 0;
% % key = @(d,b) sprintf('%d_%d',d,b);
% % addNode = @(d,b) add_node_impl(d,b);
% % 
% %     function idx = add_node_impl(d,b)
% %         k = key(d,b);
% %         if nodes.isKey(k)
% %             idx = nodes(k);
% %         else
% %             nid = nid + 1;
% %             nodes(k) = nid;
% %             rev{nid} = [d,b];
% %             idx = nid;
% %         end
% %     end
% % 
% % % Build directed edges
% % E = [];
% % for k = 1:numel(insarData)
% %     i = pairs(k,1); j = pairs(k,2);
% %     b1 = bursts(k,1); b2 = bursts(k,2);
% %     u = addNode(i,b1);
% %     v = addNode(j,b2);
% %     % forward
% %     E = [E; u, v, k, +1];
% %     % reverse
% %     E = [E; v, u, k, -1];
% % end
% % 
% % numNodes = numel(rev);
% % if numNodes < 3
% %     return;
% % end
% % 
% % % Build adjacency and edge lookup
% % adj = cell(numNodes,1);
% % edgeKey = containers.Map('KeyType','char','ValueType','any');
% % for r = 1:size(E,1)
% %     u = E(r,1); v = E(r,2); k = E(r,3); sgn = E(r,4);
% %     adj{u}(end+1) = v;
% %     edgeKey(sprintf('%d_%d',u,v)) = [k, sgn];
% % end
% % 
% % % Find triangles
% % triangles = [];
% % for u = 1:numNodes-2
% %     Nu = unique(adj{u});
% %     Nu = Nu(Nu > u);
% %     for ii = 1:numel(Nu)
% %         v = Nu(ii);
% %         Nv = unique(adj{v});
% %         Nv = Nv(Nv > v);
% %         for jj = 1:numel(Nv)
% %             w = Nv(jj);
% %             if edgeKey.isKey(sprintf('%d_%d',w,u))
% %                 triangles = [triangles; u v w];
% %             end
% %         end
% %     end
% % end
% % 
% % % Process each triangle
% % for t = 1:size(triangles,1)
% %     u = triangles(t,1); v = triangles(t,2); w = triangles(t,3);
% % 
% %     [k_uv, s_uv] = edgeKey(sprintf('%d_%d',u,v));
% %     [k_vw, s_vw] = edgeKey(sprintf('%d_%d',v,w));
% %     [k_wu, s_wu] = edgeKey(sprintf('%d_%d',w,u));
% % 
% %     A = getfield_safe(insarData(k_uv), opts.field);
% %     B = getfield_safe(insarData(k_vw), opts.field);
% %     C = getfield_safe(insarData(k_wu), opts.field);
% % 
% %     if isempty(A) || isempty(B) || isempty(C)
% %         closureMap = [];
% %     else
% %         phi = s_uv*A + s_vw*B + s_wu*C;
% %         closureMap = angle(exp(1i * phi));
% %         if opts.requireFiniteAll
% %             validMask = isfinite(A) & isfinite(B) & isfinite(C);
% %             closureMap(~validMask) = NaN;
% %         end
% %     end
% % 
% %     finiteMask = isfinite(closureMap);
% %     validFraction = sum(finiteMask(:)) / numel(finiteMask);
% %     meanBias = mean(closureMap(finiteMask), 'omitnan');
% %     stdBias  = std(closureMap(finiteMask), 'omitnan');
% % 
% %     dirs   = [rev{u}(1), rev{v}(1), rev{w}(1)];
% %     burstsThis = [rev{u}(2), rev{v}(2), rev{w}(2)];
% %     if numel(unique(burstsThis)) == 1
% %         burstOut = burstsThis(1);
% %     else
% %         burstOut = burstsThis;
% %     end
% % 
% %     insarClosureData(t).triplet       = dirs;
% %     insarClosureData(t).burst         = burstOut;
% %     insarClosureData(t).meanBias      = meanBias;
% %     insarClosureData(t).stdBias       = stdBias;
% %     insarClosureData(t).validFraction = validFraction;
% %     insarClosureData(t).closureMap    = closureMap;
% % end
% % 
% % % --- Canonical sort of results by triplet (ascending) and bursts ---
% % if ~isempty(insarClosureData)
% %     N = numel(insarClosureData);
% % 
% %     % Flag once (default false if missing)
% %     strict = isfield(opts,'sortBurstsStrict') && opts.sortBurstsStrict;
% % 
% %     % Preallocate sort keys
% %     tripletCanon = zeros(N,3);
% %     if strict
% %         burstKey = nan(N,3);   % [b1 b2 b3] aligned to sorted triplet
% %     else
% %         burstKey = nan(N,1);   % min(burst) for grouping
% %     end
% % 
% %     for n = 1:N
% %         % Canonical triplet order
% %         dirs = insarClosureData(n).triplet(:).';
% %         [dirsSorted, p] = sort(dirs);     % ascending
% %         tripletCanon(n,:) = dirsSorted;
% % 
% %         % Bursts aligned to the same permutation (if strict)
% %         b = insarClosureData(n).burst;
% %         if isscalar(b)
% %             if strict
% %                 burstKey(n,:) = [b b b];
% %             else
% %                 burstKey(n,1) = b;
% %             end
% %         else
% %             b = b(:).';
% %             b = b(p);  % align bursts with sorted dir order
% %             if strict
% %                 % pad/truncate defensively to 3 columns
% %                 m = min(3, numel(b));
% %                 burstKey(n,1:m) = b(1:m);
% %             else
% %                 burstKey(n,1) = min(b);
% %             end
% %         end
% %     end
% % 
% %     % One clear sort matrix (shape differs by 'strict')
% %     sortMat = [tripletCanon, burstKey];
% % 
% %     [~, idx] = sortrows(sortMat);
% %     insarClosureData = insarClosureData(idx);
% % end
% % 
% %     function val = getfield_safe(S, fieldname)
% %         if isfield(S, fieldname), val = S.(fieldname);
% %         else, val = []; end
% %     end
% % end
% % % 
% % % 
% % % % function insarClosureData = compute_insar_closure_bias(insarData)
% % % % %COMPUTE_INSAR_CLOSURE_BIAS Compute closure bias for valid SLC triplets and bursts
% % % % %   Handles nonuniform burst counts across pairs.
% % % % %
% % % % %   INPUT:
% % % % %     insarData - Struct array with fields:
% % % % %                   .pair           [i j] SLC indices
% % % % %                   .burst          burst index
% % % % %                   .phzReferenced  unwrapped & referenced phase map
% % % % %
% % % % %   OUTPUT:
% % % % % insarClosureData - Struct array with fields:
% % % % %               .triplet        [i j k]
% % % % %               .burst          burst index
% % % % %               .meanBias       mean closure error (rad)
% % % % %               .stdBias        std deviation of closure (rad)
% % % % %               .validFraction  fraction of pixels used in closure
% % % % %               .closureMap     closure phase map (rad)
% % % % 
% % % % fprintf('Computing phase closure bias from %d interferograms...\n', numel(insarData));
% % % % 
% % % % % Build map of all available (i, j, b) combinations
% % % % pairMap = containers.Map;
% % % % burstMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
% % % % 
% % % % for idx = 1:numel(insarData)
% % % %     i = insarData(idx).pair(1);
% % % %     j = insarData(idx).pair(2);
% % % %     b = insarData(idx).burst;
% % % %     key = sprintf('%d_%d_%d', i, j, b);
% % % %     pairMap(key) = idx;
% % % % 
% % % %     % Track available bursts per pair
% % % %     shortKey = sprintf('%d_%d', i, j);
% % % %     if ~isKey(burstMap, shortKey)
% % % %         burstMap(shortKey) = b;
% % % %     else
% % % %         burstMap(shortKey) = unique([burstMap(shortKey), b]);
% % % %     end
% % % % end
% % % % 
% % % % % Auto-detect number of unique SLCs from .pair values
% % % % allPairs = vertcat(insarData.pair);
% % % % allDirs = unique(allPairs(:));
% % % % numDirs = max(allDirs);
% % % % 
% % % % triplets = nchoosek(1:numDirs, 3);
% % % % closureData = struct([]);
% % % % c = 0;
% % % % 
% % % % for t = 1:size(triplets,1)
% % % %     i = triplets(t,1);
% % % %     j = triplets(t,2);
% % % %     k = triplets(t,3);
% % % % 
% % % %     % Retrieve available bursts for each pair
% % % %     key_ij = sprintf('%d_%d', i, j);
% % % %     key_jk = sprintf('%d_%d', j, k);
% % % %     key_ik = sprintf('%d_%d', i, k);
% % % % 
% % % %     if isKey(burstMap, key_ij) && isKey(burstMap, key_jk) && isKey(burstMap, key_ik)
% % % %         bursts_ij = burstMap(key_ij);
% % % %         bursts_jk = burstMap(key_jk);
% % % %         bursts_ik = burstMap(key_ik);
% % % % 
% % % %         % Find common bursts for which all three interferograms exist
% % % %         commonBursts = intersect(intersect(bursts_ij, bursts_jk), bursts_ik);
% % % % 
% % % %         for b = commonBursts(:)'
% % % %             key1 = sprintf('%d_%d_%d', i, j, b);
% % % %             key2 = sprintf('%d_%d_%d', j, k, b);
% % % %             key3 = sprintf('%d_%d_%d', i, k, b);
% % % % 
% % % %             if isKey(pairMap, key1) && isKey(pairMap, key2) && isKey(pairMap, key3)
% % % %                 idx_ij = pairMap(key1);
% % % %                 idx_jk = pairMap(key2);
% % % %                 idx_ik = pairMap(key3);
% % % % 
% % % %                 phi_ij = insarData(idx_ij).phzReferenced;
% % % %                 phi_jk = insarData(idx_jk).phzReferenced;
% % % %                 phi_ik = insarData(idx_ik).phzReferenced;
% % % % 
% % % %                 valid = ~isnan(phi_ij) & ~isnan(phi_jk) & ~isnan(phi_ik);
% % % %                 closure = NaN(size(phi_ij));
% % % %                 closure(valid) = wrapToPi(phi_ij(valid) + phi_jk(valid) - phi_ik(valid));
% % % % 
% % % %                 c = c + 1;
% % % %                 closureData(c).triplet = [i j k];
% % % %                 closureData(c).burst = b;
% % % %                 closureData(c).meanBias = mean(closure(valid), 'omitnan');
% % % %                 closureData(c).stdBias = std(closure(valid), 'omitnan');
% % % %                 closureData(c).validFraction = sum(valid(:)) / numel(valid);
% % % %                 closureData(c).closureMap = closure;
% % % %             end
% % % %         end
% % % %     end
% % % % end
% % % % 
% % % % insarClosureData = closureData;
% % % % fprintf('Closure bias complete: %d valid triplet-bursts.\n', numel(insarClosureData));
% % % % end

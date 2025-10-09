function idx = lookup_pair_geom(geomData, i1,b1,i2,b2)
% Return geometry index for the exact oriented pair; NaN if missing.
idx = NaN;

% --- Fast exact lookup via pair2geom
if isfield(geomData,'pair2geom') && isa(geomData.pair2geom,'containers.Map')
    keyExact = sprintf('%d_%d_%d_%d', i1,b1,i2,b2);
    if geomData.pair2geom.isKey(keyExact)
        idx = double(geomData.pair2geom(keyExact));
        return;
    end
end

% --- Slow exact fallback over pairList
if isfield(geomData,'pairList') && ~isempty(geomData.pairList)
    pl = geomData.pairList; % [i1 b1 i2 b2]
    r = find(pl(:,1)==i1 & pl(:,2)==b1 & pl(:,3)==i2 & pl(:,4)==b2, 1);
    if ~isempty(r), idx = r; return; end
end

% --- NEW 1: try reversed orientation (if geometry is symmetric in your use)
if isfield(geomData,'pair2geom') && isa(geomData.pair2geom,'containers.Map')
    keyRev = sprintf('%d_%d_%d_%d', i2,b2,i1,b1);
    if geomData.pair2geom.isKey(keyRev)
        idx = double(geomData.pair2geom(keyRev));
        return;
    end
end

% --- NEW 2: directory-agnostic reuse using self-pair representatives
% If the pair is a repeat-trajectory cross (i1≠i2) and we only computed
% self-pairs like (i1,b1)-(i1,b1) or (i2,b2)-(i2,b2), reuse one of those.
if isfield(geomData,'pair2geom') && isa(geomData.pair2geom,'containers.Map')
    keySelf1 = sprintf('%d_%d_%d_%d', i1,b1,i1,b1);
    if geomData.pair2geom.isKey(keySelf1)
        idx = double(geomData.pair2geom(keySelf1));
        return;
    end
    keySelf2 = sprintf('%d_%d_%d_%d', i2,b2,i2,b2);
    if geomData.pair2geom.isKey(keySelf2)
        idx = double(geomData.pair2geom(keySelf2));
        return;
    end
end
if isfield(geomData,'pairList') && ~isempty(geomData.pairList)
    % Same self-pair logic when only pairList exists
    r = find(pl(:,1)==i1 & pl(:,2)==b1 & pl(:,3)==i1 & pl(:,4)==b1, 1);
    if ~isempty(r), idx = r; return; end
    r = find(pl(:,1)==i2 & pl(:,2)==b2 & pl(:,3)==i2 & pl(:,4)==b2, 1);
    if ~isempty(r), idx = r; return; end
end

% --- NEW 3: optional UID fallback (works if you assign same UID to repeats)
if isfield(geomData,'idx2uid') && isfield(geomData,'pairList') && ~isempty(geomData.pairList)
    uid1 = geomData.idx2uid(i1); uid2 = geomData.idx2uid(i2);
    if ~isempty(uid1) && ~isempty(uid2) && uid1==uid2
        % Find any row in pairList with same UID on both ends and same bursts
        pl = geomData.pairList;
        sameB = (pl(:,2)==b1) & (pl(:,4)==b2);
        sameU = (geomData.idx2uid(pl(:,1))==uid1) & (geomData.idx2uid(pl(:,3))==uid2);
        r = find(sameB & sameU, 1);
        if ~isempty(r), idx = r; return; end
    end
end

% Miss
idx = NaN;
end

% function idx = lookup_pair_geom(geomData, i1,b1,i2,b2)
% % Return geometry index for the exact oriented pair; NaN if missing.
% idx = NaN;
% 
% % 0) Fast exact lookup via pair2geom (as you already had)
% if isfield(geomData,'pair2geom') && isa(geomData.pair2geom,'containers.Map')
%     keyExact = sprintf('%d_%d_%d_%d', i1,b1,i2,b2);
%     if geomData.pair2geom.isKey(keyExact)
%         idx = double(geomData.pair2geom(keyExact)); 
%         return;
%     end
% end
% 
% % 1) Slow exact fallback over pairList (as you already had)
% if isfield(geomData,'pairList') && ~isempty(geomData.pairList)
%     pl = geomData.pairList; % [i1 b1 i2 b2] per row
%     r = find(pl(:,1)==i1 & pl(:,2)==b1 & pl(:,3)==i2 & pl(:,4)==b2, 1);
%     if ~isempty(r), idx = r; return; end
% end
% 
% % 2) NEW: UID-based reuse when trajectories were repeated (directory-agnostic)
% % Requires: geomData.idx2uid(i) gives a stable UID for SLC index i.
% if isfield(geomData,'idx2uid') && isfield(geomData,'pairList') && ~isempty(geomData.pairList)
%     uid1 = geomData.idx2uid(i1);
%     uid2 = geomData.idx2uid(i2);
% 
%     % Build a lazy cache from UID pairs -> geometry row index (undirected)
%     if ~isfield(geomData,'uidpair2geom') || ~isa(geomData.uidpair2geom,'containers.Map')
%         % Create it on the struct (Matlab passes by value; return via assignin if needed)
%         uidpair2geom = containers.Map('KeyType','char','ValueType','double');
% 
%         pl = geomData.pairList; % rows: [i1 b1 i2 b2]
%         for r = 1:size(pl,1)
%             iA = pl(r,1); bA = pl(r,2);
%             iB = pl(r,3); bB = pl(r,4);
%             % Only catalog pairs with same bursts so reuse is consistent across b
%             if ~isfield(geomData,'idx2uid'), continue; end
%             uA = geomData.idx2uid(iA);
%             uB = geomData.idx2uid(iB);
% 
%             % Use undirected canonical key + bursts (keeps b-consistency)
%             keyU = sprintf('U%06d|U%06d|B%02d|B%02d', min(uA,uB), max(uA,uB), bA, bB);
%             if ~uidpair2geom.isKey(keyU)
%                 uidpair2geom(keyU) = r; % prefer first seen representative
%             end
%         end
% 
%         % Stash back onto geomData (caller keeps latest copy)
%         try
%             % If this function is nested, this preserves state outside:
%             assignin('caller','geomData', setfield(geomData,'uidpair2geom',uidpair2geom)); %#ok<SFLD>
%         catch
%             % If assignin fails (e.g., isolated test), at least use local variable
%         end
%     else
%         uidpair2geom = geomData.uidpair2geom;
%     end
% 
%     % Lookup by UID pair (undirected) with same bursts
%     keyTry = sprintf('U%06d|U%06d|B%02d|B%02d', min(uid1,uid2), max(uid1,uid2), b1, b2);
%     if uidpair2geom.isKey(keyTry)
%         idx = double(uidpair2geom(keyTry));
%         return;
%     end
% 
%     % Optional relaxed fallback: if bursts weren’t cataloged, try any burst combo
%     % (Uncomment if you want even more permissive reuse.)
%     %{
%     keysU = uidpair2geom.keys;
%     prefix = sprintf('U%06d|U%06d|', min(uid1,uid2), max(uid1,uid2));
%     hit = find(startsWith(keysU, prefix), 1);
%     if ~isempty(hit)
%         idx = double(uidpair2geom(keysU{hit}));
%         return;
%     end
%     %}
% end
% 
% % 3) As a last-ditch, try reversed orientation (if your geometry is symmetric)
% if isfield(geomData,'pair2geom') && isa(geomData.pair2geom,'containers.Map')
%     keyRev = sprintf('%d_%d_%d_%d', i2,b2,i1,b1);
%     if geomData.pair2geom.isKey(keyRev)
%         idx = double(geomData.pair2geom(keyRev)); 
%         return;
%     end
% end
% 
% % Still missing
% idx = NaN;
% end
% 
% % function idx = lookup_pair_geom(geomData, i1,b1,i2,b2)
% % % Return geometry index for the exact oriented pair; NaN if missing.
% % idx = NaN;
% % if isfield(geomData,'pair2geom') && isa(geomData.pair2geom,'containers.Map')
% %     key = sprintf('%d_%d_%d_%d', i1,b1,i2,b2);
% %     if geomData.pair2geom.isKey(key)
% %         idx = double(geomData.pair2geom(key));
% %         return;
% %     end
% % end
% % % Fallback: slow linear search (shouldn’t be needed)
% % if isfield(geomData,'pairList')
% %     pl = geomData.pairList;
% %     r = find(pl(:,1)==i1 & pl(:,2)==b1 & pl(:,3)==i2 & pl(:,4)==b2, 1);
% %     if ~isempty(r), idx = r; end
% % end
% % end

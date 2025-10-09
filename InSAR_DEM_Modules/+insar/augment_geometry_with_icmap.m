function geomData = augment_geometry_with_icmap(geomData, sarData, demData, opts)
% Adds closestIndex_master / closestIndex_slave to each geomData.sarGeometry{idx}
% Picks trajectories from slcGeomIndex that reference that idx.

if nargin<4, opts = struct(); end
if ~isfield(opts,'downsample'), opts.downsample = [8 8]; end

% convenience: get grids (X,Y,Z) once
X = demData.X; Y = demData.Y; Z = demData.dem;

usedIdx = unique([geomData.slcGeomIndex.idx]);
usedIdx = usedIdx(~isnan(usedIdx));

for uu = usedIdx(:)'
    G = geomData.sarGeometry{uu};
    % if isfield(G,'closestIndex_master') && isfield(G,'closestIndex_slave')
    %     continue;  % already done
    % end

    % mask to swath if available, else everywhere
    if isfield(G,'lookMask') && ~isempty(G.lookMask)
        mask = logical(G.lookMask);
    else
        mask = true(size(X));
    end

    % find a representative (i,b) that used this idx as 'master' and as 'slave'
    % [im, bm] = find(strcmp({geomData.slcGeomIndex.pass},'master') & [geomData.slcGeomIndex.idx]==uu, 1, 'first');
    % [isv, bs] = find(strcmp({geomData.slcGeomIndex.pass},'slave')  & [geomData.slcGeomIndex.idx]==uu, 1, 'first');
    % 
    % % fallbacks: if none found (e.g., single-trajectory case), reuse whichever exists
    % if isempty(im)
    %     [im,bm] = find([geomData.slcGeomIndex.idx]==uu, 1, 'first');
    % end
    % if isempty(isv)
    %     [isv,bs] = find([geomData.slcGeomIndex.idx]==uu, 1, 'first');
    % end
    % 
    % % trajectories for master/slave roles (these are what InSARgeometry2 used as r1/r2)
    % Pm = sarData(im).traj{bm};
    % Ps = sarData(isv).traj{bs};

    % Find all (i,b) pairs for master and slave roles matching idx == uu
    imask = strcmp({geomData.slcGeomIndex.pass}, 'master') & [geomData.slcGeomIndex.idx] == uu;
    smask = strcmp({geomData.slcGeomIndex.pass}, 'slave')  & [geomData.slcGeomIndex.idx] == uu;

    [im, bm] = find_index_pair(geomData.slcGeomIndex, imask);
    [isv, bs] = find_index_pair(geomData.slcGeomIndex, smask);

    % Fallbacks if master/slave not found
    if isempty(im)
        fallback_mask = [geomData.slcGeomIndex.idx] == uu;
        [im, bm] = find_index_pair(geomData.slcGeomIndex, fallback_mask);
    end
    if isempty(isv)
        fallback_mask = [geomData.slcGeomIndex.idx] == uu;
        [isv, bs] = find_index_pair(geomData.slcGeomIndex, fallback_mask);
    end

    % Access master/slave trajectories
    Pm = sarData(im).traj{bm};
    Ps = sarData(isv).traj{bs};


    % build fractional closest-approach maps
    icM = utils.closestIndexMap(Pm, X, Y, Z, mask, struct('downsample',opts.downsample));
    icS = utils.closestIndexMap(Ps, X, Y, Z, mask, struct('downsample',opts.downsample));

    geomData.sarGeometry{uu}.closestIndex_master = icM;
    geomData.sarGeometry{uu}.closestIndex_slave  = icS;
end
end
function [ii, jj] = find_index_pair(indexStruct, logicMask)
    ii = [];
    jj = [];
    idxFound = find(logicMask, 1, 'first');
    if ~isempty(idxFound)
        % This assumes indexStruct is organized like geomData.slcGeomIndex(ii,j)
        [ii, jj] = ind2sub(size(indexStruct), idxFound);
    end
end

% function geomData = augment_geometry_with_icmap(geomData, sarData, opts)
% % Add closest-approach index map (icMap) to each used sarGeometry{idx}.
% % opts.downsample controls speed of the nearest-segment search.
% 
% if nargin<3, opts=struct(); end
% if ~isfield(opts,'downsample'), opts.downsample=[8 8]; end
% 
% used = unique([geomData.slcGeomIndex.idx]); used = used(~isnan(used));
% for uu = used
%     G = geomData.sarGeometry{uu};
%     if isfield(G,'closestIndex'), continue; end  % already done
% 
%     % pick any SLC that references this idx to get its trajectory
%     [i1,b1] = find(geomData.slcGeomIndex.idx==uu,1,'first');
%     % you store traj per SLC:
%     P = sarData(i1).traj{b1};             % [N x 3], UTM/ENU
% 
%     % KD tree on trajectory
%     kdt = KDTreeSearcher(P);
%     [ny,nx] = size(G.slant);              % same shape as DEM grid
%     icMap = nan(ny,nx);
% 
%     ry = 1:opts.downsample(1):ny;
%     rx = 1:opts.downsample(2):nx;
%     for iy = ry
%         cols = rx(G.lookMask(iy,rx));
%         if isempty(cols), continue; end
%         Px = [G.X(iy,cols)', G.Y(iy,cols)', G.Z(iy,cols)'];  % if you keep X,Y,Z here; else use demData.*
% 
%         idx0 = knnsearch(kdt, Px);        % nearest sample indices
%         % refine on the better of the two adjacent segments, store fractional index
%         frac = nan(numel(cols),1);
%         for k = 1:numel(cols)
%             i0 = idx0(k);
%             best = struct('d2',inf,'ic',i0);
%             for s=1:2
%                 i1s = max(i0+(s-2),1); i2s = min(i1s+1, size(P,1));
%                 if i1s==i2s, continue; end
%                 A = P(i1s,:); B = P(i2s,:); AB = B-A; AP = Px(k,:)-A;
%                 t = max(0,min(1, dot(AP,AB)/max(dot(AB,AB),eps)));
%                 C = A + t*AB; d2 = sum((Px(k,:)-C).^2);
%                 if d2<best.d2, best = struct('d2',d2,'ic', i1s + t); end
%             end
%             frac(k) = best.ic;            % continuous index (e.g., 723.4)
%         end
%         icMap(iy,cols) = frac;
%     end
% 
%     % simple inpainting to full grid
%     icMap = fillmissing(icMap,'nearest');
%     geomData.sarGeometry{uu}.closestIndex = icMap;
% end
% end

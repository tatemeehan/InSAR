function icMap = closestIndexMap(P, X, Y, Z, mask, opts)
% P [N x 3] trajectory (UTM/ENU), X,Y,Z grids same size, mask logical
% opts.downsample = [dy dx] (default [8 8])

if nargin<6, opts = struct(); end
if ~isfield(opts,'downsample'), opts.downsample = [8 8]; end

[ny,nx] = size(X);
icMap = nan(ny,nx);
dy = opts.downsample(1); dx = opts.downsample(2);

% KD-tree on trajectory
kdt = KDTreeSearcher(P);

rows = 1:dy:ny;
for iy = rows
    cols = 1:dx:nx;
    cols = cols(mask(iy,cols));
    if isempty(cols), continue; end
    G = [X(iy,cols)', Y(iy,cols)', Z(iy,cols)'];

    % nearest sample index as seed
    idx0 = knnsearch(kdt, G);

    % refine on better of two adjacent segments → fractional index i+t
    frac = nan(numel(cols),1);
    N = size(P,1);
    for k = 1:numel(cols)
        i0 = idx0(k);
        best_d2 = inf; best_ic = i0;
        for s = 1:2
            i1 = max(i0+(s-2),1); i2 = min(i1+1,N);
            if i1==i2, continue; end
            A = P(i1,:); B = P(i2,:);
            AB = B-A; AP = G(k,:)-A;
            AB2 = dot(AB,AB);
            if AB2<=0, continue; end
            t = max(0,min(1, dot(AP,AB)/AB2));    % clamp to segment
            C = A + t*AB;
            d2 = sum((G(k,:)-C).^2);
            if d2 < best_d2
                best_d2 = d2;
                best_ic = i1 + t;                  % fractional index
            end
        end
        frac(k) = best_ic;
    end
    icMap(iy,cols) = frac;
end

% Inpaint to full grid (nearest works well for smooth fields)
% icMap = fillmissing(icMap,'nearest');
end

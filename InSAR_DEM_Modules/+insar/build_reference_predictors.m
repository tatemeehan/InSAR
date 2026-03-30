function pred = build_reference_predictors(G, mask)
%BUILD_REFERENCE_PREDICTORS Build pair-specific sensor-frame predictors.
%
% Outputs:
%   pred.azCoord    - normalized azimuth-like coordinate from closest-index maps
%   pred.rangeCoord - normalized slant-range-like coordinate
%   pred.mask       - valid mask used for fitting
%   pred.azRaw      - raw azimuth coordinate
%   pred.rangeRaw   - raw range coordinate

if nargin < 2 || isempty(mask)
    if isfield(G,'lookMask') && ~isempty(G.lookMask)
        mask = logical(G.lookMask);
    else
        mask = true(size(G.slant));
    end
end

% --- azimuth-like coordinate from closest indices
hasCM = isfield(G,'closestIndex_master') && ~isempty(G.closestIndex_master);
hasCS = isfield(G,'closestIndex_slave')  && ~isempty(G.closestIndex_slave);

if hasCM && hasCS
    azRaw = 0.5 * (double(G.closestIndex_master) + double(G.closestIndex_slave));
elseif hasCM
    azRaw = double(G.closestIndex_master);
elseif hasCS
    azRaw = double(G.closestIndex_slave);
else
    % fallback to row index if absolutely necessary
    [rr, ~] = ndgrid(1:size(mask,1), 1:size(mask,2));
    azRaw = double(rr);
end

% --- range-like coordinate from slant
hasS1 = isfield(G,'slant')  && ~isempty(G.slant);
hasS2 = isfield(G,'slant2') && ~isempty(G.slant2);

if hasS1 && hasS2
    rangeRaw = 0.5 * (double(G.slant) + double(G.slant2));
elseif hasS1
    rangeRaw = double(G.slant);
elseif hasS2
    rangeRaw = double(G.slant2);
else
    [~, cc] = ndgrid(1:size(mask,1), 1:size(mask,2));
    rangeRaw = double(cc);
end

% Planar fit to Azimuth
az0 = azRaw;
beamMask = mask;

% Build XY grids in image coordinates
[nr, nc] = size(az0);
[X, Y] = meshgrid(1:nc, 1:nr);

% Exclude end-of-aperture values
v = az0(beamMask & isfinite(az0));
qLow  = prctile(v, 5);
qHigh = prctile(v, 95);

fitMask = beamMask & isfinite(az0) & az0 > qLow & az0 < qHigh;

if nnz(fitMask) < 20
    azPlane = az0;
else
    A = [ones(nnz(fitMask),1), X(fitMask), Y(fitMask)];
    N = A' * A;
    if rcond(N) < 1e-10
        azPlane = az0;
    else
        b = A \ az0(fitMask);
        azPlane = nan(size(az0));
        azPlane(beamMask) = b(1) + b(2)*X(beamMask) + b(3)*Y(beamMask);
    end
end

azRaw = azPlane;
% normalize on valid mask for better conditioning
azCoord = local_normalize_field(azRaw, mask);
rangeCoord = local_normalize_field(rangeRaw, mask);

pred = struct();
pred.azCoord = single(azCoord);
pred.rangeCoord = single(rangeCoord);
pred.azRaw = single(azRaw);
pred.rangeRaw = single(rangeRaw);
pred.mask = logical(mask);

end

function out = local_normalize_field(x, mask)
x = double(x);
v = x(mask & isfinite(x));
if isempty(v)
    out = zeros(size(x));
    return;
end

mu = mean(v, 'omitnan');
sd = std(v, 0, 'omitnan');
if ~isfinite(sd) || sd < eps
    sd = 1;
end

out = (x - mu) ./ sd;
out(~isfinite(out)) = NaN;
end
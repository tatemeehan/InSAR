function fdc = doppler_from_trajectory_map(X,Y,Z, P, V, lambda, opts)
% Weighted Doppler centroid map from trajectory (UTM/ENU frame).
% X,Y,Z : DEM grids (same size)
% P,V   : [N x 3] positions & velocities (same frame as X,Y,Z)
% lambda: wavelength (m)
% opts:
%   .downsample     [dy dx] (default [4 4])
%   .method         'interp1' (default) | 'nearest'   % closest approach locator
%   .t              [N x 1] slow time (s) (optional; used if .windowTime is given)
%   .windowSamples  integer half-window on each side (default auto ~ 51)
%   .windowMeters   scalar half-window in meters (overrides samples if given)
%   .windowTime     scalar half-window in seconds (needs opts.t)
%   .weight         'gaussian' (default) | 'hann'
%   .useR2Weight    true/false (default true) multiply by 1/R^2
%   .mask           logical swath mask (optional)
%
% Returns fdc (Hz) with NaN outside mask or where undefined.

if nargin<7, opts = struct(); end
if ~isfield(opts,'downsample'),  opts.downsample  = [4 4]; end
if ~isfield(opts,'method'),      opts.method      = 'interp1'; end
if ~isfield(opts,'weight'),      opts.weight      = 'gaussian'; end
if ~isfield(opts,'useR2Weight'), opts.useR2Weight = true; end

[ny,nx] = size(X);
fdc = nan(ny,nx);

% Precompute along-track sample spacing if using meters→samples
ds = vecnorm(diff(P),2,2);
ds_med = median(ds(~isnan(ds)&isfinite(ds)));
if isempty(ds_med) || ds_med==0, ds_med = 1; end

% Decide half-window size in samples
halfW = [];
if isfield(opts,'windowMeters') && ~isempty(opts.windowMeters)
    halfW = max(3, round(opts.windowMeters / ds_med));
elseif isfield(opts,'windowTime') && ~isempty(opts.windowTime) && isfield(opts,'t') && ~isempty(opts.t)
    dt_med = median(diff(opts.t));
    if ~isfinite(dt_med) || dt_med<=0, dt_med = 1; end
    halfW = max(3, round(opts.windowTime / dt_med));
elseif isfield(opts,'windowSamples') && ~isempty(opts.windowSamples)
    halfW = max(3, round(opts.windowSamples));
else
    halfW = 51; % safe default
end

% KD-tree to find nearest sample quickly
kdt = KDTreeSearcher(P);

% Downsampled evaluation grid
ry = 1:opts.downsample(1):ny;
rx = 1:opts.downsample(2):nx;

% Optional mask
if isfield(opts,'mask') && ~isempty(opts.mask)
    mask = logical(opts.mask);
else
    mask = true(ny,nx);
end

for iy = ry
    useCols = rx(mask(iy,rx));
    if isempty(useCols), continue; end
    Px = [X(iy,useCols).', Y(iy,useCols).', Z(iy,useCols).'];

    % Seed index via nearest neighbor
    idx0 = knnsearch(kdt, Px);  % Kx1

    fd_line = nan(numel(useCols),1);

    for k = 1:numel(useCols)
        i0 = idx0(k);

        % Refine closest point on better of the two adjacent segments
        switch lower(opts.method)
            case 'nearest'
                alpha = 0; iBase = i0;  % no refinement
            otherwise % 'interp1'
                % candidates: [i0-1,i0] and [i0,i0+1]
                best = struct('d2',inf,'alpha',0,'iBase',i0);
                for s = 1:2
                    i1 = max(i0+(s-2),1); i2 = min(i1+1, size(P,1));
                    if i1==i2, continue; end
                    A = P(i1,:); B = P(i2,:);
                    AP = Px(k,:) - A;  AB = B - A; AB2 = dot(AB,AB);
                    if AB2<=0, continue; end
                    a = max(0,min(1, dot(AP,AB)/AB2));
                    C = A + a*AB; d2 = sum((Px(k,:)-C).^2);
                    if d2 < best.d2
                        best = struct('d2',d2,'alpha',a,'iBase',i1);
                    end
                end
                alpha = best.alpha; iBase = best.iBase;
        end

        % Window limits (clamped)
        i1 = max(1, iBase - halfW);
        i2 = min(size(P,1), iBase + halfW);

        % Build weights across window indices
        idxs = (i1:i2).';
        % center index (in continuous sense)
        ic = iBase + alpha;
        u = idxs - ic;  % distance in samples
        switch lower(opts.weight)
            case 'hann'
                % symmetric Hann over window length
                W = 0.5*(1 + cos(pi*u/halfW));
                W(abs(u)>halfW) = 0;
            otherwise
                % Gaussian, sigma ~ halfW/2.5
                sigma = max(1, halfW/2.5);
                W = exp(-0.5*(u./sigma).^2);
        end

        % LOS & radial velocities
        PP = P(idxs,:); VV = V(idxs,:);
        R  = Px(k,:) - PP;      Rn = sqrt(sum(R.^2,2));
        Rhat = R ./ Rn;
        vdotr = sum(VV .* Rhat, 2);

        if opts.useR2Weight
            W = W ./ max(Rn.^2, 1); % avoid div by zero
        end

        wsum = sum(W);
        if wsum <= 0 || any(~isfinite(vdotr))
            continue;
        end
        fd_line(k) = (2/lambda) * sum(W .* vdotr) / wsum;
    end

    fdc(iy,useCols) = fd_line;
end
end

% function fd = doppler_from_trajectory_map(X, Y, Z, trajPos, trajVel, lambda, opts)
% % Doppler from trajectory over a map grid.
% % X,Y,Z: grids (same size); trajPos, trajVel: [N x 3] in same Cartesian frame.
% % opts:
% %   .downsample  [dy dx] (default [4 4])
% %   .method      'nearest' (default) | 'interp1'
% %   .t           [N x 1] slow time (s) (optional, used with 'interp1')
% %
% % Returns fd (Hz), NaN where not evaluated.
% 
% if nargin<7, opts=struct(); end
% if ~isfield(opts,'downsample'), opts.downsample = [4 4]; end
% if ~isfield(opts,'method'),     opts.method     = 'interp1'; end
% 
% [ny,nx] = size(X);
% fd = nan(ny,nx);
% 
% pos = trajPos; vel = trajVel;
% N = size(pos,1);
% if size(vel,1) ~= N, error('trajPos and trajVel length mismatch'); end
% 
% useTime = isfield(opts,'t') && ~isempty(opts.t);
% if useTime
%     t = opts.t(:);
%     if numel(t) ~= N, error('opts.t must be length N'); end
% end
% 
% % Precompute for nearest/segment search
% kdt = KDTreeSearcher(pos);
% 
% % Valid pixels to evaluate (downsampled)
% ry = 1:opts.downsample(1):ny;
% rx = 1:opts.downsample(2):nx;
% 
% for iy = ry
%     Px = [X(iy,rx)', Y(iy,rx)', Z(iy,rx)'];   % Kx3 ground points
%     % 1) Nearest index (seed)
%     idx0 = knnsearch(kdt, Px);                % Kx1
% 
%     switch lower(opts.method)
%         case 'nearest'
%             Pplat = pos(idx0,:);  Vplat = vel(idx0,:);
%             R = Px - Pplat; Rhat = R ./ vecnorm(R,2,2);
%             vdotr = sum(Vplat .* Rhat, 2);
%             fd_line = (2/lambda) * vdotr;
% 
%         case 'interp1'
%             % Project onto the closest segment (idx0-1 .. idx0+1), pick the
%             % best of the two segments, then linearly interpolate P,V (and t)
%             fd_line = nan(size(idx0));
%             for k = 1:numel(idx0)
%                 i0 = idx0(k);
% 
%                 % Candidate segments: [i0-1,i0] and [i0,i0+1]
%                 segs = [max(i0-1,1), i0; i0, min(i0+1,N)];
%                 bestDist2 = inf; bestAlpha = NaN; bestI = NaN;
% 
%                 for s = 1:2
%                     i1 = segs(s,1); i2 = segs(s,2);
%                     if i1 == i2, continue; end
%                     A = pos(i1,:); B = pos(i2,:);
%                     AP = Px(k,:) - A;
%                     AB = B - A; AB2 = dot(AB,AB);
%                     if AB2 <= 0, continue; end
%                     alpha = max(0,min(1, dot(AP,AB)/AB2));  % clamp to segment
%                     C = A + alpha*AB;                        % closest point on segment
%                     d2 = sum((Px(k,:)-C).^2);
%                     if d2 < bestDist2
%                         bestDist2 = d2; bestAlpha = alpha; bestI = i1;
%                     end
%                 end
% 
%                 if ~isnan(bestAlpha)
%                     i1 = bestI; i2 = min(bestI+1,N);
%                     % Interpolate platform state at closest point
%                     Pplat = (1-bestAlpha)*pos(i1,:) + bestAlpha*pos(i2,:);
%                     Vplat = (1-bestAlpha)*vel(i1,:) + bestAlpha*vel(i2,:);
%                     % Optional time at closest approach (not used here)
%                     if useTime
%                         t_star = (1-bestAlpha)*t(i1) + bestAlpha*t(i2); %#ok<NASGU>
%                     end
%                     R = Px(k,:) - Pplat; Rhat = R / norm(R);
%                     fd_line(k) = (2/lambda) * dot(Vplat, Rhat);
%                 end
%             end
% 
%         otherwise
%             error('Unknown opts.method: %s', opts.method);
%     end
% 
%     fd(iy,rx) = fd_line;
% end
% end

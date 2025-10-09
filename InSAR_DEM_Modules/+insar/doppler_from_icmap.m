function fdc = doppler_from_icmap(P, V, lambda, geomG, opts)
% P,V : [N x 3] traj pos/vel in same frame as geom grid
% geomG.closestIndex : [ny x nx] fractional closest approach per pixel
% geomG.lookMask     : [ny x nx] logical
% opts.windowSamples / windowMeters / windowTime  (choose one)
% opts.t : slow time (needed if using windowTime)

if nargin<5, opts=struct(); end
if ~isfield(opts,'windowSamples'), opts.windowSamples = 51; end
if ~isfield(opts,'useR2Weight'),   opts.useR2Weight   = true; end

[ny,nx] = size(geomG.closestIndex);
fdc = nan(ny,nx);
mask = true(ny,nx);
if isfield(geomG,'lookMask') && ~isempty(geomG.lookMask)
    mask = geomG.lookMask;
end
% decide half-window in samples
if isfield(opts,'windowMeters') && ~isempty(opts.windowMeters)
    ds = vecnorm(diff(P),2,2); ds_med = median(ds(~isnan(ds)));
    halfW = max(3, round(opts.windowMeters / max(ds_med,1)));
elseif isfield(opts,'windowTime') && ~isempty(opts.windowTime) && isfield(opts,'t') && ~isempty(opts.t)
    dt = median(diff(opts.t)); halfW = max(3, round(opts.windowTime / max(dt,1)));
else
    halfW = max(3, round(opts.windowSamples));
end
sigma = max(1, halfW/2.5);  % for Gaussian weights

% process only masked pixels (vectorized by rows)
for iy = 1:ny
    cols = find(mask(iy,:));
    if isempty(cols), continue; end
    ic = geomG.closestIndex(iy, cols);       % fractional indices
    i0 = floor(ic);  a = ic - i0;            % linear interp weights

    % LOS at ic (interpolate P at ic, then LOS from platform to ground)
    % ground points in UTM/ENU:
    if isfield(geomG,'X'), Gx = geomG.X(iy,cols)'; Gy = geomG.Y(iy,cols)'; Gz = geomG.Z(iy,cols)'; 
    else, error('Need ground coordinates at pixel locations'); 
    end
    % interp platform position at ic
    i1 = max(1, i0); i2 = min(size(P,1), i0+1);
    P_ic = (1-a)'.*P(i1,:) + a'.*P(i2,:);           % [K x 3]
    R = [Gx Gy Gz] - P_ic;  Rn = sqrt(sum(R.^2,2)); Rhat = R ./ Rn;

    % build per-pixel window indices and weights, then average V
    f = zeros(numel(cols),1);
    for k = 1:numel(cols)
        c = ic(k);
        iC = max(1, min(size(P,1), c));
        i1w = max(1, floor(iC) - halfW);
        i2w = min(size(P,1), floor(iC) + halfW);
        idx = (i1w:i2w)';
        u = idx - c;
        w = exp(-0.5*(u./sigma).^2);
        if opts.useR2Weight
            % approximate R^2 with LOS distance at center (cheap & stable)
            w = w ./ max(Rn(k).^2, 1);
        end
        w = w / sum(w);

        Vbar = w' * V(idx,:);                    % [1 x 3]
        f(k) = (2/lambda) * dot(Vbar, Rhat(k,:));
    end
    fdc(iy, cols) = f;
end
end

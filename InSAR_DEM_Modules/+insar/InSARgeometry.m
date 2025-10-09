function [Bperp_out, r_slant_out, incident_out, lookmask_out] = ...
    InSARgeometry(Xi, Yi, DEMi, ...
    r1, r2, surfaceNormal, lookDirection, ...
    useDirectionalWeights, sigma_range, beamAngleDeg)

% Inputs:
%   Xi, Yi, DEMi       : Vectorized coordinates and elevation (Nx1)
%   r1, r2             : Tx3 transmitter positions (master/slave)
%   surfaceNormal      : Nx3 array of DEM surface normals
%   lookDirection      : 'right' or 'left' (default: 'right')
%   useDirectionalWeights : true/false to switch weighting modes (default: true)
%   sigma_range        : Gaussian weighting width for slant range [m] (default: 25)
%   beamAngleDeg       : Angular beam cutoff threshold [deg] (default: 30)

% Outputs:
%   Bperp_out          : Interpolated perpendicular baseline [m]
%   r_slant_out        : Interpolated slant range [m]
%   incident_out       : Interpolated incidence angle [deg]
%   lookmask_out       : Visibility mask (1 = valid, 0 = sidelooked)

% Handle defaults
if nargin < 7 || isempty(lookDirection), lookDirection = 'right'; end
if nargin < 8 || isempty(useDirectionalWeights), useDirectionalWeights = true; end
if nargin < 9 || isempty(sigma_range), sigma_range = 25; end
if nargin < 10 || isempty(beamAngleDeg), beamAngleDeg = 20; end

% Setup
[m,n] = size(DEMi);
DEMi = DEMi(:); Xi = Xi(:); Yi = Yi(:);
validIx = ~isnan(DEMi);
Xi = Xi(validIx); Yi = Yi(validIx); DEMi = DEMi(validIx);
surfaceNormal = surfaceNormal(validIx, :);
N = numel(Xi); T = size(r1,1);
Nout = m * n;

% Preallocate
Bperp_vals = zeros(N,T);
rslant_vals = zeros(N,T);
inc_vals = zeros(N,T);
weights = zeros(N,T);

% Heading vectors
heading = diff(r1(:,1:2),1,1);
heading(end+1,:) = heading(end,:);
vSAR = heading ./ vecnorm(heading,2,2);

% Look direction unit vectors
if strcmp(lookDirection, 'right')
    lookVec2D = [-vSAR(:,2), vSAR(:,1)]; % right-looking
else
    lookVec2D = [vSAR(:,2), -vSAR(:,1)]; % left-looking
end

% Loop through time steps
for t = 1:T
    % LOS vectors
    dX = r1(t,1) - Xi;
    dY = r1(t,2) - Yi;
    dZ = r1(t,3) - DEMi;
    losVec = [dX, dY, dZ];
    losNorm = sqrt(sum(losVec.^2, 2));
    lhat = losVec ./ losNorm;
    % Baseline and perpendicular projection
    B = r2(t,:) - r1(t,:);
    B_dot_lhat = sum(lhat .* B, 2);
    B_proj = B_dot_lhat .* lhat;
    B_perp_vec = repmat(B, N, 1) - B_proj;
    Bperp_vals(:,t) = sqrt(sum(B_perp_vec.^2,2));

    % Slant range
    rslant_vals(:,t) = losNorm;

    % Incidence angle
    dotInc = sum(lhat .* surfaceNormal, 2);
    dotInc = max(min(dotInc,1),-1);
    inc_vals(:,t) = acosd(dotInc);

    % Weighting
    if useDirectionalWeights
        w = max(B_dot_lhat, 0);  % Projected alignment
    else
        delta_r = losNorm; %- median(losNorm,'omitnan');
        w = exp(-(delta_r.^2) / (2 * sigma_range^2));
    end

    % Beam angle cutoff
    look2d = losVec(:,1:2);
    azDot = sum(look2d .* lookVec2D(t,:), 2);
    azNorm = vecnorm(look2d,2,2) * norm(lookVec2D(t,:));
    azAngle = acosd(azDot ./ azNorm);
    w(azAngle > beamAngleDeg) = 0;

    weights(:,t) = w;
end

% Normalize weights
w_sum = sum(weights,2); w_sum(w_sum==0) = eps;
normW = weights ./ w_sum;

% Weighted average
Bperp_interp = sum(Bperp_vals .* normW, 2);
rslant_interp = sum(rslant_vals .* normW, 2);
inc_interp = sum(inc_vals .* normW, 2);
lookmask_interp = any(weights > 0, 2);

% Expand to full raster
Bperp_out = nan(Nout,1); r_slant_out = nan(Nout,1);
incident_out = nan(Nout,1); lookmask_out = false(Nout,1);
Bperp_out(validIx) = Bperp_interp;
r_slant_out(validIx) = rslant_interp;
incident_out(validIx) = inc_interp;
lookmask_out(validIx) = lookmask_interp;

% Reshape
if n > 1
    Bperp_out = reshape(Bperp_out,m,n);
    r_slant_out = reshape(r_slant_out,m,n);
    incident_out = reshape(incident_out,m,n);
    lookmask_out = reshape(lookmask_out,m,n);
end
end

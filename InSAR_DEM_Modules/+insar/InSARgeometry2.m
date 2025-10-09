function [Bperp_out, r_slant_out, incident_out, lookmask_out, r2_slant_out, r2_incident_out] = ...
    InSARgeometry2(Xi, Yi, DEMi, r1, r2, surfaceNormal, lookDirection, useDirectionalWeights, sigma_range, beamAngleDeg, useBeamFractionMask, BeamFractionThreshold, incRangeDeg)

% InSAR Geometry Interpolator with Optional Beam Filtering
% ---------------------------------------------------------
% Inputs:
%   Xi, Yi, DEMi       : DEM coordinates and elevations (Nx1)
%   r1, r2             : Master/slave SAR trajectories (Tx3)
%   surfaceNormal      : DEM surface normals (Nx3)
%   lookDirection      : 'right' or 'left' (default: 'right')
%   useDirectionalWeights : If true, use directional baseline weighting
%   sigma_range        : Gaussian weighting scale (default: 250 m)
%   beamAngleDeg       : Beam angle cutoff (default: 20 deg)
%   useBeamFractionMask: If true, applies beam hit threshold (default: false)
%   BeamFractionThreshold: Percentage of Beam Hits to Include Pixel (Default: 10%)
%   incRangeDeg        : Incidence Angle Threshold for Look Mask (default: [0 90])
%
% Outputs:
%   Bperp_out          : Interpolated perpendicular baseline [m]
%   r_slant_out        : Slant range from r1 [m]
%   incident_out       : Incidence angle from r1 [deg]
%   lookmask_out       : Valid geometry mask (0 = sidelooked)
%   r2_slant_out       : Slant range from r2 [m]
%   r2_incident_out    : Incidence angle from r2 [deg]

if nargin < 7 || isempty(lookDirection), lookDirection = 'right'; end
if nargin < 8 || isempty(useDirectionalWeights), useDirectionalWeights = false; end
if nargin < 9 || isempty(sigma_range), sigma_range = 250; end
if nargin < 10 || isempty(beamAngleDeg), beamAngleDeg = 40; end
if nargin < 11 || isempty(useBeamFractionMask), useBeamFractionMask = false; end
if nargin < 12 && useBeamFractionMask == true, BeamFractionThreshold = 0.1; end
if nargin < 13 || isempty(incRangeDeg), incRangeDeg = [0 90]; end


[m,n] = size(DEMi);
DEMi = DEMi(:); Xi = Xi(:); Yi = Yi(:);
validIx = ~isnan(DEMi);
Xi = Xi(validIx); Yi = Yi(validIx); DEMi = DEMi(validIx);
surfaceNormal = surfaceNormal(validIx, :);
N = numel(Xi); T = size(r1,1); Nout = m * n;

% Preallocate
Bperp_vals = zeros(N,T); rslant_vals = zeros(N,T);
inc_vals = zeros(N,T); r2slant_vals = zeros(N,T);
r2inc_vals = zeros(N,T); weights = zeros(N,T);
azAngle_vals = zeros(N,T); beamMask_vals = false(N,T);

% Heading and look vector
heading = diff(r1(:,1:2),1,1); heading(end+1,:) = heading(end,:);
vSAR = heading ./ vecnorm(heading,2,2);
if strcmp(lookDirection, 'right')
    lookVec2D = [-vSAR(:,2), vSAR(:,1)];
else
    lookVec2D = [vSAR(:,2), -vSAR(:,1)];
end

for t = 1:T
    % Compute LOS from r1 and r2
    losVec1 = r1(t,:) - [Xi, Yi, DEMi];
    losVec2 = r2(t,:) - [Xi, Yi, DEMi];
    losNorm1 = vecnorm(losVec1, 2, 2);
    losNorm2 = vecnorm(losVec2, 2, 2);
    lhat1 = losVec1 ./ losNorm1;
    lhat2 = losVec2 ./ losNorm2;

    % Bperp and slant range
    B = r2(t,:) - r1(t,:);
    B_dot_lhat = sum(lhat1 .* B, 2);
    B_proj = B_dot_lhat .* lhat1;
    B_perp = B - B_proj;

    Bperp_vals(:,t) = vecnorm(B_perp, 2, 2);
    rslant_vals(:,t) = losNorm1;
    r2slant_vals(:,t) = losNorm2;

    % Incidence angles
    inc_vals(:,t) = acosd(min(max(sum(lhat1 .* surfaceNormal,2), -1), 1));
    r2inc_vals(:,t) = acosd(min(max(sum(lhat2 .* surfaceNormal,2), -1), 1));

    % Beam angle computation
    losXY = losVec1(:,1:2); losXY_norm = vecnorm(losXY,2,2);
    dotProd = sum(losXY .* lookVec2D(t,:),2);
    cosTheta = min(max(dotProd ./ (losXY_norm * norm(lookVec2D(t,:))), -1), 1);
    azAngle = acosd(cosTheta);
    beamMask = azAngle <= beamAngleDeg;

    azAngle_vals(:,t) = azAngle;
    beamMask_vals(:,t) = beamMask;

    % Weighting
    if useDirectionalWeights
        w = max(B_dot_lhat, 0);
    else
        w = exp(-(losNorm1.^2) / (2 * sigma_range^2));
    end
    w(~beamMask) = 0;
    weights(:,t) = w;
end

% Beam hit thresholding (optional)
if useBeamFractionMask
    minHits = ceil(BeamFractionThreshold * T);
    validBeam = sum(beamMask_vals, 2) >= minHits;
    weights(~validBeam,:) = 0;
end

% Normalize weights
wsum = sum(weights,2); wsum(wsum==0) = eps;
normW = weights ./ wsum;

% Weighted interpolation
Bperp_interp = sum(Bperp_vals .* normW, 2);
rslant_interp = sum(rslant_vals .* normW, 2);
inc_interp = sum(inc_vals .* normW, 2);
r2slant_interp = sum(r2slant_vals .* normW, 2);
r2inc_interp = sum(r2inc_vals .* normW, 2);
beamAngle_interp = sum(azAngle_vals .* normW, 2);

% Output mask based on weights
lookmask_interp = sum(weights,2) > 0.001;
% Incidence Angle Mask (0 - 90 degrees)
meanInc = (inc_interp + r2inc_interp) ./ 2;
incLo = incRangeDeg(1);
incHi = incRangeDeg(2);
lookmask_inc = isfinite(meanInc) & (meanInc >= incLo) & (meanInc <= incHi);
% Combine masks
lookmask_interp = lookmask_interp & lookmask_inc;

% Fill output grids
Bperp_out = nan(Nout,1); r_slant_out = nan(Nout,1);
incident_out = nan(Nout,1); r2_slant_out = nan(Nout,1);
r2_incident_out = nan(Nout,1); lookmask_out = false(Nout,1);

% Mask interpolated results with lookmask
Bperp_interp(~lookmask_interp)    = NaN;
rslant_interp(~lookmask_interp)   = NaN;
inc_interp(~lookmask_interp)      = NaN;
r2slant_interp(~lookmask_interp)  = NaN;
r2inc_interp(~lookmask_interp)    = NaN;

% In-fill the Raster
Bperp_out(validIx) = Bperp_interp;
r_slant_out(validIx) = rslant_interp;
incident_out(validIx) = inc_interp;
r2_slant_out(validIx) = r2slant_interp;
r2_incident_out(validIx) = r2inc_interp;
lookmask_out(validIx) = lookmask_interp;

% Reshape if grid
if n > 1
    Bperp_out = reshape(Bperp_out,m,n);
    r_slant_out = reshape(r_slant_out,m,n);
    incident_out = reshape(incident_out,m,n);
    r2_slant_out = reshape(r2_slant_out,m,n);
    r2_incident_out = reshape(r2_incident_out,m,n);
    lookmask_out = reshape(lookmask_out,m,n);
end

% Diagnostic beam angle plot
% figure(); imagesc(reshape(beamAngle_interp,m,n));
% colorbar; title('Beam Angle [deg]'); daspect([1 1 1]);
end
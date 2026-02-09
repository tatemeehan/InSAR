function [Bperp_out, r_slant_out, incident_out, lookmask_out, r2_slant_out, r2_incident_out, ...
          losBearingN_out, flatInc_out, r2_losBearingN_out, r2_flatInc_out] = ...
    InSARgeometry2_withFlatEarth(Xi, Yi, DEMi, r1, r2, surfaceNormal, lookDirection, ...
                                useDirectionalWeights, sigma_range, beamAngleDeg, ...
                                useBeamFractionMask, BeamFractionThreshold, incRangeDeg)

% InSAR Geometry Interpolator with Optional Beam Filtering + Flat-Earth Incidence + LOS Bearing-from-North
% ------------------------------------------------------------------------------------------------------
% Adds:
%   losBearingN_out     : LOS bearing (compass) ground->sensor for r1 [deg, 0..360)
%   flatInc_out         : Flat-earth incidence for r1 [deg, 0..90], using up=[0,0,1]
%   r2_losBearingN_out  : LOS bearing for r2 [deg, 0..360)
%   r2_flatInc_out      : Flat-earth incidence for r2 [deg, 0..90]
%
% Conventions:
% - Coordinates assumed local Cartesian with X=Easting, Y=Northing, Z=Up.
% - Bearing from North: 0 = North, 90 = East, increases clockwise.
% - Flat-earth incidence here is incidence-from-horizontal = asind(lhat_z) = 90 - offNadir.
%
% NOTE on bearing averaging:
% - Bearings are interpolated with a weighted circular mean (sin/cos) for robustness across 0/360 wrap.

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
Bperp_vals     = zeros(N,T);
rslant_vals    = zeros(N,T);
inc_vals       = zeros(N,T);
r2slant_vals   = zeros(N,T);
r2inc_vals     = zeros(N,T);
weights        = zeros(N,T);
azAngle_vals   = zeros(N,T);
beamMask_vals  = false(N,T);

% New preallocs (store per-time bearings and flat incidence)
losBearN_vals    = zeros(N,T);   % r1 bearing-from-North [deg 0..360)
flatInc_vals     = zeros(N,T);   % r1 flat-earth incidence [deg]
r2losBearN_vals  = zeros(N,T);   % r2 bearing-from-North [deg 0..360)
r2flatInc_vals   = zeros(N,T);   % r2 flat-earth incidence [deg]

% Heading and look vector (for beam-angle gating)
heading = diff(r1(:,1:2),1,1); heading(end+1,:) = heading(end,:);
vSAR = heading ./ vecnorm(heading,2,2);
if strcmp(lookDirection, 'right')
    lookVec2D = [-vSAR(:,2), vSAR(:,1)];
else
    lookVec2D = [vSAR(:,2), -vSAR(:,1)];
end

% Flat-earth "up" (assumes +Z up)
up = [0 0 1];

for t = 1:T
    % Compute LOS from r1 and r2 (ground -> sensor)
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

    Bperp_vals(:,t)   = vecnorm(B_perp, 2, 2);
    rslant_vals(:,t)  = losNorm1;
    r2slant_vals(:,t) = losNorm2;

    % True (terrain-normal) incidence angles
    inc_vals(:,t)   = acosd(min(max(sum(lhat1 .* surfaceNormal,2), -1), 1));
    r2inc_vals(:,t) = acosd(min(max(sum(lhat2 .* surfaceNormal,2), -1), 1));

    % NEW: LOS bearing-from-North (compass) for ground->sensor
    % For EN (X=E, Y=N): bearing = atan2(E, N) in degrees, wrapped to [0,360)
    bear1 = atan2d(losVec1(:,1), losVec1(:,2));   % atan2(E, N)
    bear2 = atan2d(losVec2(:,1), losVec2(:,2));
    losBearN_vals(:,t)   = mod(bear1, 360);
    r2losBearN_vals(:,t) = mod(bear2, 360);

    % NEW: Flat-earth incidence (incidence-from-horizontal)
    lhat1z = sum(lhat1 .* up, 2);
    lhat2z = sum(lhat2 .* up, 2);
    lhat1z = min(max(lhat1z, -1), 1);
    lhat2z = min(max(lhat2z, -1), 1);
    flatInc_vals(:,t)    = asind(lhat1z);
    r2flatInc_vals(:,t)  = asind(lhat2z);

    % Beam angle computation (azimuthal offset from nominal look direction)
    losXY = losVec1(:,1:2); losXY_norm = vecnorm(losXY,2,2);
    dotProd = sum(losXY .* lookVec2D(t,:),2);
    cosTheta = min(max(dotProd ./ (losXY_norm * norm(lookVec2D(t,:))), -1), 1);
    azAngle = acosd(cosTheta);
    beamMask = azAngle <= beamAngleDeg;

    azAngle_vals(:,t)  = azAngle;
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

% Weighted interpolation (linear quantities)
Bperp_interp     = sum(Bperp_vals    .* normW, 2);
rslant_interp    = sum(rslant_vals   .* normW, 2);
inc_interp       = sum(inc_vals      .* normW, 2);
r2slant_interp   = sum(r2slant_vals  .* normW, 2);
r2inc_interp     = sum(r2inc_vals    .* normW, 2);
beamAngle_interp = sum(azAngle_vals  .* normW, 2);

flatInc_interp    = sum(flatInc_vals    .* normW, 2);
r2flatInc_interp  = sum(r2flatInc_vals  .* normW, 2);

% NEW: Weighted circular mean for bearing (robust across wrap)
ang1 = deg2rad(losBearN_vals);      % NxT
ang2 = deg2rad(r2losBearN_vals);

sin1 = sum(sin(ang1) .* normW, 2);
cos1 = sum(cos(ang1) .* normW, 2);
losBearN_interp = mod(rad2deg(atan2(sin1, cos1)), 360);

sin2 = sum(sin(ang2) .* normW, 2);
cos2 = sum(cos(ang2) .* normW, 2);
r2losBearN_interp = mod(rad2deg(atan2(sin2, cos2)), 360);

% Output mask based on weights
lookmask_interp = sum(weights,2) > 0.001;

% Incidence Angle Mask (0 - 90 degrees) using mean TRUE incidence (your original logic)
meanInc = (inc_interp + r2inc_interp) ./ 2;
incLo = incRangeDeg(1);
incHi = incRangeDeg(2);
lookmask_inc = isfinite(meanInc) & (meanInc >= incLo) & (meanInc <= incHi);

% Combine masks
lookmask_interp = lookmask_interp & lookmask_inc;

% Fill output grids
Bperp_out          = nan(Nout,1);
r_slant_out        = nan(Nout,1);
incident_out       = nan(Nout,1);
r2_slant_out       = nan(Nout,1);
r2_incident_out    = nan(Nout,1);
losBearingN_out    = nan(Nout,1);
flatInc_out        = nan(Nout,1);
r2_losBearingN_out = nan(Nout,1);
r2_flatInc_out     = nan(Nout,1);
lookmask_out       = false(Nout,1);

% Mask interpolated results with lookmask
Bperp_interp(~lookmask_interp)        = NaN;
rslant_interp(~lookmask_interp)       = NaN;
inc_interp(~lookmask_interp)          = NaN;
r2slant_interp(~lookmask_interp)      = NaN;
r2inc_interp(~lookmask_interp)        = NaN;
losBearN_interp(~lookmask_interp)     = NaN;
flatInc_interp(~lookmask_interp)      = NaN;
r2losBearN_interp(~lookmask_interp)   = NaN;
r2flatInc_interp(~lookmask_interp)    = NaN;

% In-fill the raster
Bperp_out(validIx)          = Bperp_interp;
r_slant_out(validIx)        = rslant_interp;
incident_out(validIx)       = inc_interp;
r2_slant_out(validIx)       = r2slant_interp;
r2_incident_out(validIx)    = r2inc_interp;
losBearingN_out(validIx)    = losBearN_interp;
flatInc_out(validIx)        = flatInc_interp;
r2_losBearingN_out(validIx) = r2losBearN_interp;
r2_flatInc_out(validIx)     = r2flatInc_interp;
lookmask_out(validIx)       = lookmask_interp;

% Reshape if grid
if n > 1
    Bperp_out          = reshape(Bperp_out,m,n);
    r_slant_out        = reshape(r_slant_out,m,n);
    incident_out       = reshape(incident_out,m,n);
    r2_slant_out       = reshape(r2_slant_out,m,n);
    r2_incident_out    = reshape(r2_incident_out,m,n);
    losBearingN_out    = reshape(losBearingN_out,m,n);
    flatInc_out        = reshape(flatInc_out,m,n);
    r2_losBearingN_out = reshape(r2_losBearingN_out,m,n);
    r2_flatInc_out     = reshape(r2_flatInc_out,m,n);
    lookmask_out       = reshape(lookmask_out,m,n);
end

% Diagnostic plots (optional)
% figure(); imagesc(reshape(beamAngle_interp,m,n)); colorbar; title('Beam Angle [deg]'); daspect([1 1 1]);
end

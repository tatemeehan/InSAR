%% ================================================================
%  Synthetic Trajectory + InSAR Geometry Planning
% ================================================================

%% ----------------------------------------------------------------
% 1. Paths and Processing Parameters
% ----------------------------------------------------------------
function InSAR_Planning(paths, params, AGL, verticalShift, shiftDistance, dx)
% paths.traj  = '/Users/rdcrlcjw/Desktop/InSAR-main/Trajectories/Alpine_Garden.txt';
% paths.dem   = '/Users/rdcrlcjw/Desktop/InSAR-main/DEMs/mtwash_DEM_3.tif';
% paths.cr    = '/Users/rdcrlcjw/Desktop/InSAR-main/CRs/MCS_CRs_2026.csv';
% paths.cmap  = '/Users/rdcrlcjw/Desktop/InSAR-main/colormaps/RdYlBu.csv';
% paths.out   = '/Users/rdcrlcjw/Desktop/InSAR-main/Export_Figures';
% 
% params.c = 0.3;      % wave speed (m/ns)
% params.f = 1.3;      % radar frequency (GHz)
% params.lambda = params.c / params.f;
% 
% AGL = 100;           % flight altitude above ground (m)
% verticalShift = 20;  % vertical separation between tracks
% shiftDistance = 0;   % horizontal track shift
% dx = 20;             % interpolation spacing (m)

params.lambda = params.c / params.f;

%% ----------------------------------------------------------------
% 2. Load and Build Trajectory
% ----------------------------------------------------------------

tmp = readtable(paths.traj);

MCStrajectory = [ ...
    tmp.Var9(2:3), ...
    tmp.Var10(2:3), ...
    [tmp.Var11(1); tmp.Var11(1)] ];

% Convert lat/lon → UTM
[X,Y] = utils.deg2utm(MCStrajectory(:,1),MCStrajectory(:,2));

testTraj = [X,Y,MCStrajectory(:,3)];
testTraj(:,3) = testTraj(:,3) + AGL;

% Heading
delta   = diff(testTraj);
heading = mod(atan2d(delta(:,1),delta(:,2)),360);

% Second trajectory (vertical offset)
testTraj2 = testTraj;
testTraj2(:,1) = testTraj(:,1) + shiftDistance*cosd(heading-90);
testTraj2(:,2) = testTraj(:,2) + shiftDistance*sind(heading-90);
testTraj2(:,3) = testTraj(:,3) + verticalShift;

%% ----------------------------------------------------------------
% 3. Interpolate Trajectories
% ----------------------------------------------------------------

[arclen,~] = utils.arclength(testTraj(:,1),testTraj(:,2),testTraj(:,3),'pchip');
trajAxis = 0:dx/arclen:1;
r1 = utils.interparc(trajAxis,testTraj(:,1),testTraj(:,2),testTraj(:,3),'pchip');

[arclen,~] = utils.arclength(testTraj2(:,1),testTraj2(:,2),testTraj2(:,3),'pchip');
trajAxis = 0:dx/arclen:1;
r2 = utils.interparc(trajAxis,testTraj2(:,1),testTraj2(:,2),testTraj2(:,3),'pchip');

%% ----------------------------------------------------------------
% 4. Load DEM
% ----------------------------------------------------------------

[demData.dem, demData.R, ~, ~, demData.lon, demData.lat, ...
 demData.X, demData.Y, demData.EPSG] = io.readLidarTif(paths.dem);

demData.dem(demData.dem < 0) = NaN;

% Surface normals
[demData.surfaceNormal,demData.aspect,demData.slope] = ...
    utils.compute_surface_normals(demData.dem,demData.lat,demData.lon,demData.EPSG);

% Flat-earth assumption for CR incidence
demData.surfaceNormal = repmat([0 0 1],length(demData.surfaceNormal),1);

%% ----------------------------------------------------------------
% 5. Compute InSAR Geometry
% ----------------------------------------------------------------

Xi = demData.X(:);
Yi = demData.Y(:);

tic
[baseline,slantRange,incidence,lookmask,slantRange2,incidence2,...
 losBearing,flatIncidence,losBearing2,flatIncidence2] = ...
 insar.InSARgeometry2_withFlatEarth(Xi,Yi,demData.dem,r1,r2,demData.surfaceNormal);
toc

CRtilt = flatIncidence2 - asind(1/sqrt(3));

meanSlantRange = median(slantRange(~isnan(slantRange)),'all');  % in meters
meanIncidence  = median(incidence(~isnan(incidence)),'all');    % in degrees

lambda = params.lambda;  % radar wavelength
hDesired = 0.1;          % example: desired vertical height resolution in meters

% Compute recommended Bperp
BperpRecommended = (lambda*meanSlantRange)/(4*pi*sin(deg2rad(meanIncidence))*hDesired);
fprintf('Recommended B⊥ based on data: %.2f m\n', BperpRecommended);

%% ----------------------------------------------------------------
% 6. Corner Reflector Alignment
% ----------------------------------------------------------------

CR = io.read_corner_reflectors(paths.cr);
crName = {'CR1','CR2'};

k = 5;
weightType = 'gaussian';
sigma = 1;

numCR = numel(crName);
CR_idxList = cell(numCR,1);
CR_flatIncidence = zeros(numCR,1);
CRbearing = zeros(numCR,1);

for c = 1:numCR
    
    [crIdxs,~,~,~,~] = ...
        utils.get_cr_dem_indices(CR,crName{c},demData,k,weightType,sigma);
    
    CR_idxList{c} = crIdxs;
    
    CR_flatIncidence(c) = mean(flatIncidence(crIdxs));
    CRbearing(c) = mean(losBearing(crIdxs));
    
end

CR.losElevationAngle = CR_flatIncidence;
CR.bearing = CRbearing;
CR.baseTiltDeg = CR_flatIncidence - asind(1/sqrt(3));
CR.ix = CR_idxList;

%% ----------------------------------------------------------------
% 7. Antenna Angle Optimization
% ----------------------------------------------------------------

out = insar.optimize_antenna_angle( ...
    flatIncidence2,slantRange2,lookmask,'WeightMode','none');%,'CR',CR);

%% ----------------------------------------------------------------
% 8. Coherence Model
% ----------------------------------------------------------------

kz = insar.compute_vertical_wavenumber(params.lambda,baseline,slantRange,incidence);

sigmaZ = 0.3;
coherenceMod = exp(-0.5*(sigmaZ*kz).^2);

%% ----------------------------------------------------------------
% 9. Load Colormap
% ----------------------------------------------------------------

cmap = csvread('/Users/rdcrlcjw/Desktop/InSAR-main/colormaps/RdYlBu.csv');
cmap = flipud(cmap);

X = demData.X(1,:);
Y = demData.Y(:,1);

%% ----------------------------------------------------------------
% 10. Geometry Figure
% ----------------------------------------------------------------

figure

% --- Slant Range
subplot(1,4,1)
plotTerrainBackground(demData)
hold on
imagesc(X./1000,Y./1000,slantRange,'AlphaData',0.625)
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
colormap([[1 1 1];cmap])
colorbar; clim([0 1200])
title('a) Slant Range')
formatMapAxes

% --- Incidence
subplot(1,4,2)
plotTerrainBackground(demData)
hold on
imagesc(X./1000,Y./1000,incidence,'AlphaData',0.625)
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
colormap([[1 1 1];cmap])
colorbar; clim([0 90])
title('b) Incidence Angle')
formatMapAxes

% --- Baseline
subplot(1,4,3)
plotTerrainBackground(demData)
hold on
imagesc(X./1000,Y./1000,baseline,'AlphaData',0.625)
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
colormap([[1 1 1];cmap])
colorbar; clim([min(baseline(~isnan(baseline)),[],'all'), max(baseline(~isnan(baseline)),[],'all')])
title('c) B_\perp')
formatMapAxes

% --- Coherence
subplot(1,4,4)
plotTerrainBackground(demData)
hold on
imagesc(X./1000,Y./1000,coherenceMod,'AlphaData',0.625)
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
colormap([[1 1 1];cmap])
colorbar; clim([0 1])
title('d) Coherence')
formatMapAxes

% Antenna annotation
antennaStr = sprintf('Optimal Antenna Angle: %.2f°',out.tiltOptDeg_discrete);

annotation('textbox',[0.32 0.94 0.36 0.05], ...
    'String',antennaStr,...
    'EdgeColor','k',...
    'BackgroundColor','white',...
    'HorizontalAlignment','center',...
    'FontSize',14,'FontWeight','bold');

%% ----------------------------------------------------------------
% 11. Corner Reflector Figure
% ----------------------------------------------------------------

figure

subplot(1,2,1)
plotTerrainBackground(demData)
hold on
imagesc(X./1000,Y./1000,CRtilt,'AlphaData',0.625)
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
colormap([[1 1 1];cmap])
colorbar; clim([-45 45])
title('CR Tilt')
formatMapAxes

subplot(1,2,2)
plotTerrainBackground(demData)
hold on
imagesc(X./1000,Y./1000,losBearing2,'AlphaData',0.625)
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
colormap([[1 1 1];cmap])
colorbar; clim([0 359])
title('CR Bearing')
formatMapAxes

%% ----------------------------------------------------------------
% 12. Export Outputs
% ----------------------------------------------------------------

% KMZ CR tilt
io.georaster2kmz(CRtilt,demData.R, ...
    fullfile(paths.out,'Camas2026CRtiltAngle240AGL.kmz'),...
    'CRS',demData.EPSG,'Name','Tilt angle (deg)','AddLegend',true,'SampleStep',10);

% Export trajectory
crs = projcrs(demData.EPSG);
[lat,lon] = projinv(crs,r1(:,1),r1(:,2));

kmlwriteline(fullfile(paths.out,'SAR_trajectory.kml'),lat,lon,r1(:,3), ...
    'AltitudeMode','relativeToSeaLevel','Color','yellow','LineWidth',3);

% Export incidence
io.georaster2kmz(incidence2,demData.R,...
    fullfile(paths.out,'IncidenceAngle.kmz'),...
    'CRS',demData.EPSG,'Name','Incidence angle','AddLegend',true,'SampleStep',10);

geotiffwrite(fullfile(paths.out,'IncidenceAngle.tif'), ...
             incidence2,demData.R,'CoordRefSysCode',demData.EPSG);

%% ================================================================
% Helper Functions
% ================================================================

function plotTerrainBackground(demData)

imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000, ...
((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5*demData.slope))

colormap(bone)
utils.freezeColors

end


function formatMapAxes

daspect([1 1 1])
set(gca,'YDir','normal','FontName','serif','FontWeight','bold','FontSize',14)

end

end
% Create Synthetic Trajectory
% WCP Entire Length
% testTraj = [715282, 4901860, 1665; 714100,4903660, 1675];
% WCP North Transect 2024 100 m AGL
% testTraj = [7.152854372026424e+05,4.903041641240263e+06,1.751923531197726e+03;...
    % 7.143103406620359e+05,4.903743909067830e+06,1.749685156920199e+03];
GCSnorth = [44.250062,-114.313299];
GCSsouth = [ 44.242464,-114.307302];
[GCSnorthX,GCSnorthY] = utils.deg2utm(GCSnorth(1),GCSnorth(2));
[GCSsouthX,GCSsouthY] = utils.deg2utm(GCSsouth(1),GCSsouth(2));

% WCP North Transect 2025
testLatLon = [44.2473663 -114.3045473 1656.75;44.2549514 -114.3152641 1656.75];
% WCP South Transect 2025
testLatLon = [44.2397457 -114.3046707 1664.630;44.2487267 -114.3130258 1664.63];
% MCS SnoTel 2026
tmp = readtable('E:\MCS\MCS_CR_SM_25\CarSar_GLSAR-120.txt');
tmp([3,5,15,17,27,29],:) = [];
MCStrajectory = [tmp.Var9(3:4),tmp.Var10(3:4),[tmp.Var11(1);tmp.Var11(1)]];
testLatLon = MCStrajectory;
[X,Y]=utils.deg2utm(testLatLon(:,1),testLatLon(:,2));
testTraj = [X,Y,testLatLon(:,3)];
% Set Flight AGL
AGL = 240;%-7.5;
testTraj(:,3) = testTraj(:,3) + AGL;
delta = diff(testTraj);
% Compute the heading angle
heading = mod(atan2d(delta(:,1), delta(:,2)),360);
% Shift Traj Coordinates
shiftDistance = 0; 
testTraj2(:,1) = testTraj(:,1)+shiftDistance.*cosd(heading-90);
testTraj2(:,2) = testTraj(:,2)+shiftDistance.*sind(heading-90);
testTraj2(:,3) = testTraj(:,3)+0;
% heading = atan2d(sind(testLL(2,2)-testLL(1,2)) * cosd(testLL(1,2)), cosd(testLL(1,1)) * sind(testLL(1,2)) - sind(testLL(1,1)) * cosd(testLL(1,2)) * cosd(testLL(2,2)-testLL(1,2)))
% Interpolate Trajectory
[arclen,seglen] = utils.arclength(testTraj(:,1),testTraj(:,2),testTraj(:,3),'pchip');
dx = 10;
trajAxis  = [0:dx./arclen:1];
r1 = utils.interparc(trajAxis,testTraj(:,1),testTraj(:,2),testTraj(:,3),'pchip');
% Interpolate Trajectory 2
[arclen,seglen] = utils.arclength(testTraj2(:,1),testTraj2(:,2),testTraj2(:,3),'pchip');
dx = 10;
trajAxis  = [0:dx./arclen:1];
r2 = utils.interparc(trajAxis,testTraj2(:,1),testTraj2(:,2),testTraj2(:,3),'pchip');
%% Load DEM Data
% DEM Paths
% demData.dir = 'E:\WCP\0924Campaign\LiDAR';
% demData.demFn = '240911_White_Cloud_Preserve_1m_DTM.tif';
% demData.vegFn  = 'WCPvegetationHeight.tif';

demData.dir = 'E:\MCS\MCS032024\lidar\tif';
demData.demFn = '20240320_MCS_REFDEM_WGS84_CarSAR_UAS_50cm.tif';
demData.vegFn  = '20240320_MCS_canopyHeight_CarSAR_UAS_50cm.tif';

% Processing Parameters
params.c          = 0.3;      % Wave speed sonstant (m/ns)
params.f          = 1.3;      % Frequency (GHz)
% params.f          = 3.2;      % Frequency (GHz)

params.lambda     = params.c / params.f;    % Radar wavelength (m)

[demData.dem, demData.R, ~, ~, demData.lon, demData.lat,demData.X, demData.Y, demData.EPSG] =...
    io.readLidarTif(fullfile(demData.dir, demData.demFn));
[demData.veg, ~, ~, ~, ~, ~,~, ~, ~] =...
    io.readLidarTif(fullfile(demData.dir, demData.vegFn));

% Compute Surface Normals
% Surface Normals
[demData.surfaceNormal,demData.aspect,demData.slope] = ...
    utils.compute_surface_normals(demData.dem, demData.lat, demData.lon, demData.EPSG);
% For Corner Reflector Incidence assume flat Earth n_hat = [0,0,1];
demData.surfaceNormal = repmat([0,0,1],length(demData.surfaceNormal),1);
%% Compute Synthetic InSAR Geometry
dem = demData.dem; Xi = demData.X(:);Yi = demData.Y(:); surfaceNormal = demData.surfaceNormal;
lambda = params.lambda;
tic;
% [baseline, slantRange, incidence, lookmask] = insar.InSARgeometry(Xi, Yi, dem, r1, r2, surfaceNormal);
% [baseline, slantRange, incidence, lookmask, slantRange2, incidence2] = ...
%     insar.InSARgeometry2(Xi, Yi, dem, ...
%     r1, r2, surfaceNormal);

[baseline, slantRange, incidence, lookmask, slantRange2, incidence2, ...
          losBearing, flatIncidence, losBearing2,flatIncidence2] = ...
    insar.InSARgeometry2_withFlatEarth(Xi, Yi, dem, r1, r2, surfaceNormal);
toc;

%% Corner Reflector Alignment
% Corner Reflectors
% MCS
crDir = 'E:\MCS\MCS_CR_SM_25';
crFn = 'MCS_CRs_2025.csv';
crPath = fullfile(crDir,crFn);
crName = {'CR1','CR2'};

k = 5; % Nearest Neighborscr
weightType = 'gaussian'; % Use Gaussian Weighting Kernel
sigma = 1; % meters

CR = io.read_corner_reflectors(crPath);

numCRs = numel(crName);
CR_idxList = cell(numCRs,1);
CR_flatIncidence = zeros(numCRs,1);
CRbearing = zeros(numCRs,1);
for c = 1:numCRs
    [crIdxs, ~, xCR0, yCR0, zCR0] = utils.get_cr_dem_indices(CR, crName{c}, demData, k, weightType, sigma);
    CR_idxList{c} = crIdxs;
    CR_flatIncidence(c) = mean(flatIncidence(crIdxs));
    CRbearing(c) = mean(losBearing(crIdxs));
end
CR.losElevationAngle = CR_flatIncidence;
CR.bearing = CRbearing;
CR.baseTiltDeg = CR_flatIncidence - asind(1./sqrt(3));
CR.ix = CR_idxList;

%% Antenna Angle Optimization

out = insar.optimize_antenna_angle(flatIncidence, slantRange, lookmask, 'CR',CR);

%% Coherence Modeling for Baseline Flights
% Compute Vertical Wave Number
kz = insar.compute_vertical_wavenumber(lambda, baseline, slantRange, incidence);
% Model the Coherence
sigmaZ = 0.35;
coherenceMod = exp(-0.5.*(sigmaZ.*kz).^2);

L_az = 3.5;
% B_crit = compute_critical_baseline(lambda, slantRange, L_az);
% Bcrit = compute_critical_baseline(lambda, slantRange, incidence, 1);
%% Figure
% Colorbars
cmap = csvread('./colormaps/RdYlBu.csv');
cmap = flipud(cmap);
X = demData.X(1,:);
Y = demData.Y(:,1);
figure();
subplot(1,3,3)
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; 
hold on;
% hI = imagesc(X./1000,Y./1000,baseline,'AlphaData',lookmask);hc = colorbar;daspect([1,1,1]);%title('B\perp')
hI = imagesc(X./1000,Y./1000,baseline,'AlphaData',0.625);hc = colorbar;daspect([1,1,1]);%title('B\perp')
hold on;
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
% colormap('bone')
ylabel(hc,'B\perp (m)','fontname','serif','fontweight','bold','fontsize',14)
xlabel('Easting (km)');%ylabel('Northing (km)');
title('c)                                                    ','FontSize',16)
% colormap(cmap)
colormap([[1 1 1];cmap])
clim([0 5])
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
set(gca,'YTickLabel',[]);

subplot(1,3,2)
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; hold on;
% hI = imagesc(X./1000,Y./1000,incidence,'AlphaData',lookmask);hc = colorbar;daspect([1,1,1]);%title('Incidence')
hI = imagesc(X./1000,Y./1000,incidence,'AlphaData',0.625);hc = colorbar;daspect([1,1,1]);%title('Incidence')
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
ylabel(hc,'Incidence Angle (\circ)','fontname','serif','fontweight','bold','fontsize',14)
xlabel('Easting (km)');%ylabel('Northing (km)');
title('b)                                                    ','FontSize',16)


% colormap('bone')
% colormap(cmap)
colormap([[1 1 1];cmap])
clim([0 90])
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
set(gca,'YTickLabel',[]);

subplot(1,3,1)
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; hold on;
% hI = imagesc(X./1000,Y./1000,slantRange,'AlphaData',lookmask);
hI = imagesc(X./1000,Y./1000,slantRange,'AlphaData',0.625);
plot(r1(:,1)./1000,r1(:,2)./1000,'k','LineWidth',2)
hc = colorbar;daspect([1,1,1]);%title('Slant Range')
% colormap('bone')
ylabel(hc,'Slant Range (m)','fontname','serif','fontweight','bold','fontsize',14)
xlabel('Easting (km)');ylabel('Northing (km)');
title('a)                                                   ','FontSize',16)
% colormap(cmap)
colormap([[1 1 1];cmap])

clim([0 1200])
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
% subplot(1,4,4)
% hI = imagesc(X./1000,Y./1000,lookmask);colorbar;daspect([1,1,1]);title('Look Mask')
% colormap('bone')
% ax = ancestor(hI, 'axes');
% ax.XAxis.Exponent = 0;
% xtickformat('%.1f')
% ax.YAxis.Exponent = 0;
% ytickformat('%.1f')
% set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
exportgraphics(gcf, ['E:\WCP\0925campaign\InSAR','\','geometrysmallerTerrainfontbiggeryetline.png'],'ContentType','image','Resolution',300);

%% Export Trajectory
% lat = trajLat(:);       % degrees
% lon = trajLon(:);       % degrees
% alt = trajAlt(:);       % meters (MSL if using 'relativeToSeaLevel')
% 
% fn = fullfile(sarData.dir,'sar_trajectory.kml');
% kmlwriteline(fn, lat, lon, alt, ...
%     'Name','SAR Trajectory', ...
%     'Color','yellow', 'Alpha',0.9, 'LineWidth',3, ...
%     'AltitudeMode','relativeToSeaLevel');   % or 'relativeToGround'

crs = projcrs(demData.EPSG);         % e.g., 32612 for UTM 12N
[lat, lon] = projinv(crs, r1(:,1), r1(:,2));
alt = r1(:,3);                           % meters
kmlwriteline(fullfile('E:\WCP\0925campaign\InSAR','WCP2025StestTrajectory120AGL.kml'), lat, lon, alt, 'AltitudeMode','relativeToSeaLevel','Name','SAR Trajectory', ...
    'Color','yellow', 'Alpha',0.9, 'LineWidth',3);
% Write kml raster
io.georaster2kmz(incidence2,demData.R,fullfile('E:\WCP\0925campaign\InSAR','WCP2025SfineTrajectoryIncidenceAngle120AGL.kmz'),'CRS', demData.EPSG,'Name', 'Incidence angle (deg)','AddLegend',true,'SampleStep',10 );
% Write Geotiff
geotiffwrite(fullfile('E:\WCP\0925campaign\InSAR','WCP2025NorthIncidenceAngle120AGL.tif'),incidence2,demData.R,'CoordRefSysCode',demData.EPSG);
%%
tic;
% [baseline, slantRange, incidence, lookmask] = insar.InSARgeometry(Xi, Yi, dem, r1, r2, surfaceNormal);
[baseline10, slantRange10, incidence10, lookmask10, slantRange210, incidence210] = ...
    insar.InSARgeometry2(Xi, Yi, dem, ...
    r1, r2, surfaceNormal);
toc;
L_az = 3.5;
% B_crit = compute_critical_baseline(lambda, slantRange, L_az);
Bcrit10 = compute_critical_baseline(lambda, slantRange10, incidence10, 1);
%%
tic;
% [baseline, slantRange, incidence, lookmask] = insar.InSARgeometry(Xi, Yi, dem, r1, r2, surfaceNormal);
[baseline15, slantRange15, incidence15, lookmask15, slantRange215, incidence215] = ...
    insar.InSARgeometry2(Xi, Yi, dem, ...
    r1, r2, surfaceNormal);
toc;
L_az = 3.5;
% B_crit = compute_critical_baseline(lambda, slantRange, L_az);
Bcrit15 = compute_critical_baseline(lambda, slantRange15, incidence15, 1);
%% 12.5 m baseline
tic;
% [baseline, slantRange, incidence, lookmask] = insar.InSARgeometry(Xi, Yi, dem, r1, r2, surfaceNormal);
[baseline12, slantRange12, incidence12, lookmask12, slantRange212, incidence212] = ...
    insar.InSARgeometry2(Xi, Yi, dem, ...
    r1, r2, surfaceNormal);
toc;
L_az = 3.5;
% B_crit = compute_critical_baseline(lambda, slantRange, L_az);
Bcrit12 = compute_critical_baseline(lambda, slantRange12, incidence12, 1);

%% 20 m vertical baseline (WCP2024)
tic;
% [baseline20, slantRange20, incidence20, lookmask20] = insar.InSARgeometry(Xi, Yi, dem, r1, r2, surfaceNormal);
[baseline20, slantRange20, incidence20, lookmask20, slantRange220, incidence220] = ...
    insar.InSARgeometry2(Xi, Yi, dem, ...
    r2, r1, surfaceNormal);
toc;
L_az = 3.5;
% Bcrit20 = compute_critical_baseline(lambda, slantRange, L_az);
Bcrit20 = compute_critical_baseline(lambda, slantRange20, incidence20, 1);
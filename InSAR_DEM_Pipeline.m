%% InSAR DEM Generation Pipeline
clear; close all; clc;

%% --------------------- USER SETTINGS ---------------------
% Paths
slcDir  = 'E:\WCP\0924Campaign\GLSAR_Data\20240913_north\';
mliDir  = slcDir;
demDir  = 'E:\WCP\0924Campaign\LiDAR\';
trajDir = 'E:\WCP\0924Campaign\GLSAR_TRAJ\09132024\';
figOut  = [slcDir, 'figures\'];

% Filenames
slcFiles = { ...
    'CRREL_TMA_1ms_200MHz_20240913T223727_TMA_TX2_RX2_V_V_V.slc', ...
    'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX2_RX2_V_V_V.slc', ...
    'CRREL_TMA_1ms_200MHz_20240913T223727_TMA_TX1_RX1_H_H_H.slc', ...
    'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX1_RX1_H_H_H.slc' };
mliFiles = { ...
    'CRREL_TMA_1ms_200MHz_20240913T223727_TMA_TX2_RX2_V_V_V.mli_geo.tif', ...
    'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX2_RX2_V_V_V.mli_geo.tif' };
demFn  = '240911_White_Cloud_Preserve_1m_DTM.tif';
vegFn  = 'WCPvegetationHeight.tif';
trajFn = '20240913_Lift4_fixed.txt';

% Processing Parameters
c          = 0.3;      % Wave speed sonstant (m/ns)
f          = 1.3;      % Frequency (GHz)
lambda     = c / f;    % Radar wavelength (m)
quality    = 0.3;      % Coherence threshold for unwrapping
dx         = 1;        % Trajectory resampling resolution (m)
filterSize = 5;        % Interferogram filter window size (pixels)

% Data Export
isWriteGeoTiff = 0;

% Colorbars
cmap = csvread('./colormaps/RdYlBu.csv');
cmap = flipud(cmap);

%% --------------------- LOAD MLI & DEM ---------------------
[mli1, R, ~, ~, lon, lat, utmX, utmY, EPSG] = io.readLidarTif(fullfile(mliDir, mliFiles{1}));
[mli2, ~, ~, ~, ~, ~, ~, ~, ~] = io.readLidarTif(fullfile(mliDir, mliFiles{2}));
[dem, Rdem, ~, ~, ~, ~, ~, ~, ~] = io.readLidarTif(fullfile(demDir, demFn));
[veg, Rveg, ~, ~, ~, ~, ~, ~, ~] = io.readLidarTif(fullfile(demDir, vegFn));

[Xq, Yq] = worldGrid(R);
dem = utils.interp_to_radar_grid(dem, Rdem, Xq, Yq);
veg = utils.interp_to_radar_grid(veg, Rveg, Xq, Yq);
nanMask = dem < 0 | veg > 100 | isnan(dem);
dem(nanMask) = NaN; veg(nanMask) = NaN;
Xi = utmX(:); Yi = utmY(:); DEMi = dem(:);

%% --------------------- LOAD SLC DATA ---------------------
frameSize = size(mli1);
slcs = cell(1, 4);
for i = 1:4
    slcs{i} = io.read_slc(fullfile(slcDir, slcFiles{i}), frameSize, nanMask);
end

%% --------------------- TRAJECTORY PROCESSING ---------------------
[r1, r2] = io.load_trajectory_and_indices(trajFn, trajDir, ...
    slcFiles{1}(1:end-5), slcFiles{2}(1:end-5), dx);

%% --------------------- SURFACE NORMAL ---------------------
[surfaceNormal,aspect,slope] = utils.compute_surface_normals(dem, lat, lon, EPSG);

%% --------------------- GEOMETRY COMPUTATION ---------------------
tic;
% [baseline, slantRange, incidence, lookmask] = insar.InSARgeometry(Xi, Yi, dem, r1, r2, surfaceNormal);
[baseline, slantRange, incidence, lookmask, slantRange2, incidence2] = ...
    insar.InSARgeometry2(Xi, Yi, dem, ...
    r1, r2, surfaceNormal);
toc;
%% --------------------- SINGLE LOOK PROCESSING ------------------------
amp = cell(1, length(slcFiles)); pow = amp; db = amp;
for ii = 1:length(slcFiles)
    if ii == 1 || ii == 3
        tmpIncidence = incidence;
    elseif ii == 2 || ii == 4
        tmpIncidence = incidence2;
    end
    % Calibration Parameters
    % CR power is ~25 db P = 10^2.5
    % Trihedral RCS
    a = 1;
    P = 10^2.533;
    pixelArea = R.SampleSpacingInWorldX.*R.SampleSpacingInWorldY;
    % Generate Single-look Images & Apply True Geometric Terrain Correction
    [amp{ii}, pow{ii}, db{ii}, meta] = insar.compute_sli(slcs{ii},P,a,lambda,pixelArea,tmpIncidence,[],{'pixelarea','gamma0'});
    % Multi-Look
    pow{ii} = medfilt2(pow{ii},[filterSize,filterSize]);
    db{ii} = medfilt2(db{ii},[filterSize,filterSize]);
    % Apply Look Direcection Mask
    amp{ii}(~lookmask) = NaN; pow{ii}(~lookmask) = NaN; db{ii}(~lookmask) = NaN;
end

%% --------------------- INTERFEROGRAMS & COHERENCE ---------------------
[phzVV, corVV] = insar.compute_interferogram(slcs{1}, slcs{2}, quality, filterSize, true);
[phzHH, corHH] = insar.compute_interferogram(slcs{3}, slcs{4}, quality, filterSize, true);

phzVV(~lookmask) = NaN; corVV(~lookmask) = NaN;
phzHH(~lookmask) = NaN; corHH(~lookmask) = NaN;

% Unwrapping
tic;
phzUnwrappedVV = insar.region_grow_unwrap_multiseed(phzVV, corVV, 0.9);
phzUnwrappedHH = insar.region_grow_unwrap_multiseed(phzHH, corHH, 0.9);
toc;

% Simulate InSAR Phase
% [simPhzUnwrapped, simPhzWrapped] = insar.simulate_insar_phase_from_dem(dem,baseline,slantRange,lambda,incidence);
[simPhzUnwrapped,simPhzWrapped] = insar.simulate_insar_phase_from_range(lambda, slantRange, slantRange2);
[simCor] = insar.phase_quality(simPhzWrapped,5);
%% --------------------- DEM COMPUTATION ---------------------
demVV = insar.insar_phase_to_dem(phzUnwrappedVV, lambda, slantRange, baseline, incidence);
demHH = insar.insar_phase_to_dem(phzUnwrappedHH, lambda, slantRange, baseline, incidence);
demVVHH = demVV-demHH;
%% --------------------- Penetration COMPUTATION ---------------------
% Note Same as DEM Computation Due to Backprojection Focusing (Sign
% Reversed)
penetrationVV = insar.phase_to_penetration(phzUnwrappedVV, lambda, baseline, slantRange, incidence);
penetrationHH = insar.phase_to_penetration(phzUnwrappedHH, lambda, baseline, slantRange, incidence);
%% --------------------- VISUALIZATION ---------------------
figure;
imagesc(utmX(1,:)/1000, utmY(:,1)/1000, ...
    ((cosd(aspect+45)+sind(aspect+45)).*sind(2.5*slope))); 
colormap(bone); utils.freezeColors; hold on;

hI = imagesc(utmX(1,:)/1000, utmY(:,1)/1000, demVV, 'AlphaData', 0.625);
daspect([1,1,1]); colormap([[1 1 1]; cmap]);  hc = colorbar;
ylabel(hc,'InSAR Signal Penetration VV - HH (m)','fontname','serif','fontweight','bold','fontsize',14)
title('White Clouds Preserve: 09/13/24');
xlabel('Easting (km)'); ylabel('Northing (km)');
set(gca, 'YDir','normal', 'FontName','serif', 'FontSize', 14, 'FontWeight','bold');
clim([-2.5 2.5]);

% exportgraphics(gcf, fullfile(figOut, 'glsarWCPinsarDEMdifferenceVVHH.png'), 'Resolution', 300);

%% Synthetic Phase Images
figure();
subplot(1,2,1)
imagesc(utmX(1,:)/1000, utmY(:,1)/1000, ...
    ((cosd(aspect+45)+sind(aspect+45)).*sind(2.5*slope))); 
colormap(bone); utils.freezeColors; hold on;

hI = imagesc(utmX(1,:)/1000, utmY(:,1)/1000, simPhzUnwrapped, 'AlphaData', 0.625);
daspect([1,1,1]); colormap([[1 1 1]; cmap]);  hc = colorbar;
ylabel(hc,'Synthetic Unwrapped Phase (rad)','fontname','serif','fontweight','bold','fontsize',14)
title('White Clouds Preserve: 09/13/24');
xlabel('Easting (km)'); ylabel('Northing (km)');
set(gca, 'YDir','normal', 'FontName','serif', 'FontSize', 14, 'FontWeight','bold');
clim([-500 0]);

subplot(1,2,2)
imagesc(utmX(1,:)/1000, utmY(:,1)/1000, ...
    ((cosd(aspect+45)+sind(aspect+45)).*sind(2.5*slope))); 
colormap(bone); utils.freezeColors; hold on;

hI = imagesc(utmX(1,:)/1000, utmY(:,1)/1000, simPhzWrapped, 'AlphaData', 0.625);
daspect([1,1,1]); colormap([[1 1 1]; cmap]);  hc = colorbar;
ylabel(hc,'Synthetic Wrapped Phase (rad)','fontname','serif','fontweight','bold','fontsize',14)
title('White Clouds Preserve: 09/13/24');
xlabel('Easting (km)'); ylabel('Northing (km)');
set(gca, 'YDir','normal', 'FontName','serif', 'FontSize', 14, 'FontWeight','bold');
clim([-pi pi]);

% exportgraphics(gcf, fullfile(figOut, 'syntheticInSARphaseWCP.png'), 'Resolution', 300);
%% Export GeoTiff Data
if isWriteGeoTiff
    geotiffwrite([slcDir,'WCPsignalPenetrationVV20240913T223727.tif'],demVV,R,'CoordRefSysCode',EPSG);
    geotiffwrite([slcDir,'WCPsignalPenetrationHH20240913T223727.tif'],demHH,R,'CoordRefSysCode',EPSG);
    geotiffwrite([slcDir,'WCPsignalPenetrationVV-HH20240913T223727.tif'],demVVHH,R,'CoordRefSysCode',EPSG);
end
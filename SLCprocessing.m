%% SLC Processing and Visualization
clear; close all; clc;

%% --------------------- USER SETTINGS ---------------------
% 2026 MCS GLSAR
% sarData(1).date = '20260223';
% sarData(1).dir = 'E:\MCS\MCS022326\GLSAR\lift1\slc';
% sarData(1).trajDir = 'E:\MCS\MCS022326\GLSAR\lift1\traj';
% sarData(1).slcFiles = {'lband-mcs_2026_TMA_2ms_200MHz_20260223T224040_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'lband-mcs_2026_TMA_2ms_200MHz_20260223T224040_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(2).date = '20260226';
% sarData(2).dir = 'E:\MCS\MCS022626\GLSAR\lift2\slc';
% sarData(2).trajDir = 'E:\MCS\MCS022626\GLSAR\lift2\traj';
% sarData(2).slcFiles = {'lband-mcs_2026_TMA_2ms_200MHz_20260226T184834_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'lband-mcs_2026_TMA_2ms_200MHz_20260226T184834_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% 2026 Camas GLSAR
% sarData(1).date = '20260219';
% sarData(1).dir = 'E:\CAMAS\GLSAR\20260219\lift4';
% sarData(1).trajDir = 'E:\CAMAS\GLSAR\20260219\lift4';
% sarData(1).slcFiles = {'camas_2026_02_19_GLSAR_240m_lift4_TMA_2ms_200MHz_20260219T223457_TMA_TX2_RX2_V_V_V.slc',...
%     'camas_2026_02_19_GLSAR_240m_shifty_lift4_TMA_2ms_200MHz_20260219T224106_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'camas_2026_02_19_GLSAR_240m_lift4_TMA_2ms_200MHz_20260219T223457_TMA_TX2_RX2_V_V_V.mli_geo.tif',...
%     'camas_2026_02_19_GLSAR_240m_shifty_lift4_TMA_2ms_200MHz_20260219T224106_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
%
% % % 240-240-240 Triplet
% sarData(1).date = '20260219';
% sarData(1).dir = 'E:\CAMAS\GLSAR\20260219\lift4';
% sarData(1).trajDir = 'E:\CAMAS\GLSAR\20260219\lift4';
% sarData(1).slcFiles = {'camas_2026_02_19_GLSAR_240m_lift4_TMA_2ms_200MHz_20260219T223457_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'camas_2026_02_19_GLSAR_240m_lift4_TMA_2ms_200MHz_20260219T223457_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(2).date = '20260219';
% sarData(2).dir = 'E:\CAMAS\GLSAR\20260219\lift3\slc';
% sarData(2).trajDir = 'E:\CAMAS\GLSAR\20260219\lift3\traj';
% sarData(2).slcFiles = {'camas_2026_02_19_GLSAR_240m_lift3_TMA_2ms_200MHz_20260219T220308_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'camas_2026_02_19_GLSAR_240m_lift3_TMA_2ms_200MHz_20260219T220308_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(3).date = '20260219';
% sarData(3).dir = 'E:\CAMAS\GLSAR\20260219\lift2\slc';
% sarData(3).trajDir = 'E:\CAMAS\GLSAR\20260219\lift2\traj';
% sarData(3).slcFiles = {'camas_2026_02_19_GLSAR_240m_lift2_TMA_2ms_200MHz_20260219T205715_TMA_TX2_RX2_V_V_V.slc'};
% sarData(3).mliFiles = {'camas_2026_02_19_GLSAR_240m_lift2_TMA_2ms_200MHz_20260219T205715_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% GLSAR 120 AGL Triplet
sarData(1).date = '20260219';
sarData(1).dir = 'E:\CAMAS\GLSAR\20260219\lift1\slc';
sarData(1).trajDir = 'E:\CAMAS\GLSAR\20260219\lift1\traj';
sarData(1).slcFiles = {'camas_2026_02_19_GLSAR_120m_lift1_TMA_2ms_200MHz_20260219T203126_TMA_TX2_RX2_V_V_V.slc'};
sarData(1).mliFiles = {'camas_2026_02_19_GLSAR_120m_lift1_TMA_2ms_200MHz_20260219T203126_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

sarData(2).date = '20260219';
sarData(2).dir = 'E:\CAMAS\GLSAR\20260219\lift5\slc';
sarData(2).trajDir = 'E:\CAMAS\GLSAR\20260219\lift5\traj';
sarData(2).slcFiles = {'camas_2026_02_19_GLSAR_120m_lift5_TMA_2ms_200MHz_20260219T231031_TMA_TX1_RX1_H_H_H.slc',...
    'camas_2026_02_19_GLSAR_120m_lift5_TMA_2ms_200MHz_20260219T231637_TMA_TX1_RX1_H_H_H.slc'};
sarData(2).mliFiles = {'camas_2026_02_19_GLSAR_120m_lift5_TMA_2ms_200MHz_20260219T231031_TMA_TX1_RX1_H_H_H.mli_geo.tif',...
    'camas_2026_02_19_GLSAR_120m_lift5_TMA_2ms_200MHz_20260219T231637_TMA_TX1_RX1_H_H_H.mli_geo.tif'};

% WCP GLSAR 091024
% sarData(1).date = '20250910';
% sarData(1).dir = 'E:\WCP\0925campaign\GLSAR\South\20250910\20250910T184409';
% sarData(1).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250910\20250910T184409';
% sarData(1).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250910T184409_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250910T184409_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(2).date = '20250910';
% sarData(2).dir = 'E:\WCP\0925campaign\GLSAR\South\20250910\20250910T190936';
% sarData(2).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250910\20250910T190936';
% sarData(2).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250910T190936_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250910T190936_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% sarData(2).date = '20250912';
% sarData(2).dir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\20250912T205604';
% sarData(2).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\20250912T205604';
% sarData(2).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T205604_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T205604_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(2).date = '20250912';
% sarData(2).dir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_2_120-5m-acrosstrackoffset\20250912T212016';
% sarData(2).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_2_120-5m-acrosstrackoffset\20250912T212016';
% sarData(2).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T212016_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T212016_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% Vertical Baseline 120-110 
% sarData(1).date = '20250912';
% sarData(1).dir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\20250912T210102';
% sarData(1).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\20250912T210102';
% sarData(1).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T210102_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T210102_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(2).date = '20250912';
% sarData(2).dir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\20250912T205604';
% sarData(2).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\20250912T205604';
% sarData(2).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T205604_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T205604_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% sarData(2).date = '20250912';
% sarData(2).dir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_2_120-5m-acrosstrackoffset\20250912T212016';
% sarData(2).trajDir = 'E:\WCP\0925campaign\GLSAR\South\20250912\lift_2_120-5m-acrosstrackoffset\20250912T212016';
% sarData(2).slcFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T212016_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'WCP_2025_TMA_2ms_200MHz_20250912T212016_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% MCS CarSAR 031824
% sarData(1).date = '20240319';
% sarData(1).dir = 'E:\MCS\MCS031924\CarSAR\processed\redux';
% sarData(1).trajDir = 'E:\MCS\MCS031924\CarSAR\processed\redux';
% sarData(1).slcFiles = {'SP_200MHz_V_6_20240319_161920_V_R3_3_V.slc','SP_200MHz_V_6_20240319_210757_V_R3_3_V.slc'};
% sarData(1).mliFiles = {'SP_200MHz_V_6_20240319_161920_V_R3_3_V.mli_geo.tif','SP_200MHz_V_6_20240319_210757_V_R3_3_V.mli_geo.tif'};

% MCS 021225 - 022025
% sarData(1).date = '20250212';
% sarData(1).dir = 'E:\MCS\MCS021225\GLSAR';
% sarData(1).trajDir = 'E:\MCS\MCS021225\GLSAR';
% sarData(1).slcFiles = {'CRREL_TMA_1ms_200MHz_20250212T192953_TMA_TX2_RX2_V_V_V.slc','CRREL_TMA_1ms_200MHz_20250212T200646_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'CRREL_TMA_1ms_200MHz_20250212T192953_TMA_TX2_RX2_V_V_V.mli_geo.tif','CRREL_TMA_1ms_200MHz_20250212T200646_TMA_TX2_RX2_V_V_V.mli_geo.tif'};
% 
% sarData(2).date = '20250220';
% sarData(2).dir = 'E:\MCS\MCS022225\GLSAR';
% sarData(2).trajDir = 'E:\MCS\MCS022225\GLSAR';
% sarData(2).slcFiles = {'CRREL_TMA_1ms_200MHz_20250220T175134_TMA_TX2_RX2_V_V_V.slc','CRREL_TMA_1ms_200MHz_20250220T183403_TMA_TX2_RX2_V_V_V.slc'};
% sarData(2).mliFiles = {'CRREL_TMA_1ms_200MHz_20250220T175134_TMA_TX2_RX2_V_V_V.mli_geo.tif','CRREL_TMA_1ms_200MHz_20250220T183403_TMA_TX2_RX2_V_V_V.mli_geo.tif'};

% sarData(1).date = '20240913';
% sarData(1).dir = 'E:\WCP\0924Campaign\GLSAR_Data\20240913_north';
% sarData(1).trajDir = 'E:\WCP\0924Campaign\GLSAR_TRAJ\09132024';
% sarData(1).slcFiles = {'CRREL_TMA_1ms_200MHz_20240913T223727_TMA_TX2_RX2_V_V_V.slc', ...
%     'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX2_RX2_V_V_V.slc'};
% sarData(1).mliFiles = {'CRREL_TMA_1ms_200MHz_20240913T223727_TMA_TX2_RX2_V_V_V.mli_geo.tif', ...
%     'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX2_RX2_V_V_V.mli_geo.tif' };

% WCP HH, HH, HV, triplet
% sarData(1).date = '20240913';
% sarData(1).dir = 'E:\WCP\0924Campaign\GLSAR_Data\20240913_north';
% sarData(1).trajDir = 'E:\WCP\0924Campaign\GLSAR_TRAJ\09132024';
% sarData(1).slcFiles = {'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX1_RX1_H_H_H.slc', ...
%     'CRREL_TMA_1ms_200MHz_20240913T232028_TMA_TX1_RX1_H_H_H.slc','CRREL_TMA_1ms_200MHz_20240913T232028_TMA_TX1_RX2_H_V_H.slc'};
% sarData(1).mliFiles = {'CRREL_TMA_1ms_200MHz_20240913T223129_TMA_TX1_RX1_H_H_H.mli_geo.tif', ...
%     'CRREL_TMA_1ms_200MHz_20240913T232028_TMA_TX1_RX1_H_H_H.mli_geo.tif','CRREL_TMA_1ms_200MHz_20240913T232028_TMA_TX1_RX2_H_V_H.mli_geo.tif'};

% Define input SLC files and directories
% % sarData(1).date     = '20240909';
% % sarData(1).dir      = 'E:\WCP\0924Campaign\carSAR\20240909\processed';
% % sarData(1).slcFiles = {'WCP_100MHz_V_6_20240909_215709_V_R3_3_H.slc'};
% % sarData(1).mliFiles = {'WCP_100MHz_V_6_20240909_215709_V_R3_3_H.mli_geo.tif'};
% % 
% % sarData(2).date     = '20240913';
% % sarData(2).dir      = 'E:\WCP\0924Campaign\carSAR\20240913\processed';
% % sarData(2).slcFiles = {'WCP_100MHz_V_6_20240913_173044_V_R3_3_H.slc'};
% % sarData(2).mliFiles = {'WCP_100MHz_V_6_20240913_173044_V_R3_3_H.mli_geo.tif'};
% 
% sarData(1).date     = '20240909';
% sarData(1).dir      = 'E:\WCP\0924Campaign\carSAR\20240909\processed';
% sarData(1).slcFiles = {'WCP_100MHz_V_6_20240909_220844_V_R3_3_H.slc'};
% sarData(1).mliFiles = {'WCP_100MHz_V_6_20240909_220844_V_R3_3_H.mli_geo.tif'};
% 
% sarData(2).date     = '20240911';
% sarData(2).dir      = 'E:\WCP\0924Campaign\carSAR\20240911\processed';
% sarData(2).slcFiles = {'WCP_100MHz_V_6_20240911_165621_V_R3_3_H.slc'};
% sarData(2).mliFiles = {'WCP_100MHz_V_6_20240911_165621_V_R3_3_H.mli_geo.tif'};
% 
% sarData(3).date     = '20240912';
% sarData(3).dir      = 'E:\WCP\0924Campaign\carSAR\20240912\processed';
% sarData(3).slcFiles = {'WCP_100MHz_V_6_20240912_183732_V_R3_3_H.slc'};
% sarData(3).mliFiles = {'WCP_100MHz_V_6_20240912_183732_V_R3_3_H.mli_geo.tif'};
% 
% sarData(4).date     = '20240913';
% sarData(4).dir      = 'E:\WCP\0924Campaign\carSAR\20240913\processed';
% sarData(4).slcFiles = {'WCP_100MHz_V_6_20240913_173044_V_R3_3_H.slc'};
% sarData(4).mliFiles = {'WCP_100MHz_V_6_20240913_173044_V_R3_3_H.mli_geo.tif'};

% DEM Paths
% demData.dir = 'E:\WCP\0925campaign\DEM';
% demData.demFn = '20250915_WCP_1m_DTM.tif';
% demData.vegFn  = 'WCPvegetationHeight.tif';

% demData.dir = 'E:\WCP\0924Campaign\LiDAR';
% demData.demFn = '240911_White_Cloud_Preserve_1m_DTM.tif';
% demData.vegFn  = 'WCPvegetationHeight.tif';

% demData.dir = 'E:\MCS\MCS032024\lidar\tif';
% demData.demFn = '20240320_MCS_REFDEM_WGS84_CarSAR_UAS_50cm.tif';
% demData.vegFn  = '20240320_MCS_canopyHeight_CarSAR_UAS_50cm.tif';

% CAMAS
demData.dir = 'E:\CAMAS\DEM';
demData.demFn = '2026_Camas_South.tif';
demData.vegFn  = 'WCPvegetationHeight.tif';

% MCS 2026
% demData.dir = 'E:\MCS\MCS022326\dem';
% demData.demFn = '2026_MCS_Trim.tif';
% demData.vegFn  = 'WCPvegetationHeight.tif';

% Processing Parameters
params.c          = 0.3;      % Wave speed sonstant (m/ns)
params.f          = 1.3;      % Frequency (GHz)
params.lambda     = params.c / params.f;    % Radar wavelength (m)
params.quality    = 0.75;     % Coherence threshold for unwrapping
params.dx         = 5;        % Trajectory resampling resolution (m)
params.filterSize = 5;        % Multilook filter window size (pixels)

% Data Export (Under Construction)
% Out Direcrtory
outDir = 'E:\CAMAS\GLSAR\20260219\Export120';

% Colorbars
cmap = csvread('./colormaps/RdYlBu.csv');
cmap = flipud(cmap);
storeCmap = cmap;

%% ---------- Load DEM and Vegetation Rasters ----------
% Read once, interpolate to MLI grid later
[demRaw, Rdem, ~, ~, ~, ~, ~, ~, ~] = io.readLidarTif(fullfile(demData.dir, demData.demFn));
% [vegRaw, Rveg, ~, ~, ~, ~, ~, ~, ~] = io.readLidarTif(fullfile(demData.dir, demData.vegFn));

%% ---------- Load Each MLI Frame & Store DEM Data ----------
fprintf('Loading MLI - SLC - DEM - Data\n')
tic
for ii = 1:numel(sarData)
    mliDir = sarData(ii).dir;

    for jj = 1:numel(sarData(ii).mliFiles)
        % Read MLI
        mliPath = fullfile(mliDir, sarData(ii).mliFiles{jj});
        if jj == 1
        [sarData(ii).mli{jj}, R, ~, ~, lon, lat, utmX, utmY, epsg] = io.readLidarTif(mliPath);
        if ii == 1
        % Store projection + coords
        demData.R       = R;
        demData.lon     = lon;
        demData.lat     = lat;
        demData.X       = utmX;
        demData.Y       = utmY;
        demData.EPSG    = epsg;
        demData.frameSize = size(sarData(ii).mli{jj});
        demData.pixelArea = R.SampleSpacingInWorldX.*R.SampleSpacingInWorldY;

        % Interpolate DEM + vegetation to MLI raster grid
        [Xq, Yq] = worldGrid(R);
        demInterp = utils.interp_to_radar_grid(demRaw, Rdem, Xq, Yq);
        % vegInterp = utils.interp_to_radar_grid(vegRaw, Rveg, Xq, Yq);

        % NaN masking
        % nanMask = demInterp < 0 | vegInterp > 100 | isnan(demInterp);
        nanMask = demInterp < 0  | isnan(demInterp);
        demInterp(nanMask) = NaN;
        % vegInterp(nanMask) = NaN;

        % Store DEM + veg + mask
        demData.dem = demInterp;
        % demData.veg = vegInterp;
        demData.nanMask = nanMask;
        end
        else
            [sarData(ii).mli{jj}, ~, ~, ~, ~, ~, ~, ~, ~] = io.readLidarTif(mliPath);
        end
    end
end
clear("demInterp","vegInterp","demRaw","vegRaw","Rdem","Rveg","Xq","Yq","nanMask","mliDir","lat","lon","utmX","utmY","R","epsg","mliPath")
%% ---------- Load SLC Data ----------
for ii = 1:numel(sarData)
    mliDir = sarData(ii).dir;

    for jj = 1:numel(sarData(ii).slcFiles)
        slcPath = fullfile(mliDir, sarData(ii).slcFiles{jj});
        sarData(ii).slc{jj} = io.read_slc(slcPath, demData.frameSize, demData.nanMask);
        sarData(ii).polarization{jj} = [slcPath(end-8),slcPath(end-6)];
    end
end
clear("mliDir","slcPath")
fprintf('Data Loaded. Elapsed time: %.2f seconds\n', toc);
%% --------------------- TRAJECTORY PROCESSING ---------------------
fprintf('Interpolating - Aligning - Truncating - Trajectories\n')
tic
trajDir = sarData(ii).dir;
for ii = 1:numel(sarData)
    if isfield(sarData, 'trajDir')
        trajDir = sarData(ii).trajDir;
    else
        trajDir = sarData(ii).dir;
    end
    for jj = 1:numel(sarData(ii).slcFiles)
        slcBase = sarData(ii).slcFiles{jj}(1:end-5);  % remove ".slc"
        % Trajectory Paths
        sarData(ii).trajPath{jj} = fullfile(trajDir, [slcBase, 'pos.dat.geod']);
        % Load Geodetic Trajectory and Resample
        [sarData(ii).traj{jj}, sarData(ii).velocity{jj}, sarData(ii).slowTime{jj}, sarData(ii).azimuthAxis{jj}, sarData(ii).trajLength{jj}] = io.load_geodetic_trajectory(sarData(ii).trajPath{jj}, params.dx);
    end
end
% ---------- Align and Truncate Trajectories ----------
[sarData] = utils.align_and_truncate_trajectories(sarData, params.dx);
clear("slcBase")
fprintf('Trajectory processing complete. Elapsed time: %.2f seconds\n', toc);

%% --------------------- GEOMETRY COMPUTATION ---------------------
% Must fix for instance of one SLC
fprintf('Computing SAR Geometry\n')
tic;
% Surface Normals
[demData.surfaceNormal,demData.aspect,demData.slope] = ...
    utils.compute_surface_normals(demData.dem, demData.lat, demData.lon, demData.EPSG);
% SAR Geometry
% [geomData] = insar.computeSARgeometry(sarData, demData);
% 1) Build all pairs (default)
geomData = insar.computeSARgeometry(sarData, demData);
% 
% % 2) Build only a short list of pairs
% pairs = [1 1 2 1;   3 1 4 1];        % [i1 b1 i2 b2] rows
% geomData = insar.computeSARgeometry(sarData, demData, struct('mode','selected','pairs',pairs));
% 
% % 3) Build by trajectory tolerance (meters between trajectories)
% geomData = insar.computeSARgeometry(sarData, demData, struct('mode','trajtol','trajTol',0.6));

dopplerOpts.downsample = [1,1];
geomData = insar.augment_geometry_with_icmap(geomData, sarData, demData, dopplerOpts);
dopplerOpts2 = struct('windowMeters',50, 'useR2Weight',true, 'smoothVel',51);
% % Doppler Frequency Centroids
% sarData = insar.compute_doppler_centroids_from_icmap( ...
%               sarData, geomData, demData, params.lambda, ...
%               'windowMeters', 25, 'useR2Weight', true, 'smoothVel', 25, 'force', false);
fprintf('SAR Geometry Computed. Elapsed time: %.2f seconds\n', toc);

% Geometry Plot
isGeometryPlot= 0;
if isGeometryPlot
figure();
subplot(1,4,1)
imagesc(geomData.sarGeometry{geomData.pairResults.geomIndex}.Bperp);colorbar;daspect([1,1,1]);title('B\perp')
subplot(1,4,2)
imagesc(geomData.sarGeometry{geomData.pairResults.geomIndex}.incidence);colorbar;daspect([1,1,1]);title('Incidence')
subplot(1,4,3)
imagesc(geomData.sarGeometry{geomData.pairResults.geomIndex}.slant);colorbar;daspect([1,1,1]);title('Slant Range')
subplot(1,4,4)
imagesc(geomData.sarGeometry{geomData.pairResults.geomIndex}.lookMask);colorbar;daspect([1,1,1]);title('Look Mask')
end
%% -------------------- Doppler Frequency Corrections ---------------------
% fprintf('Doppler Centroid Computation\n')
% tic
% for ii = 1:numel(sarData)
%     for jj = 1:numel(sarData(ii).slcFiles)
%         dopOpts.t = sarData(ii).slowTime{jj};
%         dopOpts.downsample = [1,1];
%         dopOpts.mask = geomData.sarGeometry{geomData.slcGeomIndex(ii).idx}.lookMask;
%         sarData(ii).dopplerCentroid{jj} = insar.doppler_from_trajectory_map(...
%             demData.X, demData.Y, demData.dem, sarData(ii).traj{jj}, sarData(ii).velocity{jj}, params.lambda, dopOpts);
%     end
% end
% fprintf('Doppler Centroids Computed. Elapsed time: %.2f seconds\n', toc);

%% -------------------- Corner Reflector RCS Calibration ------------------
fprintf('RCS Calibration\n')
tic
% Load CR info WCP
% % 2025
% crDir = 'E:\WCP\0925campaign\CR';
% crFn = 'WCP-2025-09-CR-Positions.csv';
% crPath = fullfile(crDir,crFn);
% crName = {'CRS1','CRS2'};
% 2024
% crDir = 'E:\WCP\0924Campaign\CR';
% crFn = 'WCP090924_CR_positions-t8.csv';
% crPath = fullfile(crDir,crFn);
% % Choose CR
% crName = {'CRC1','CRC2'};
% % MCS
% crDir = 'E:\MCS\MCS_CR_SM_25';
% % crFn = 'MCS_CRs_2025.csv';
% % crFn = 'MCS_CRs_2024.csv';
% crFn = 'MCS_CRs_2026_Corrected.csv';
% crPath = fullfile(crDir,crFn);
% crName = {'CR1','CR2'};

% Camas 2026
crDir = 'E:\CAMAS\CRs';
crFn = 'Camas_2026_CRs_field.csv';
crPath = fullfile(crDir,crFn);
crName = {'CR1','CR2','CR3','CR4'};

% crDir = [];
% crFn = [];
% crPath = fullfile(crDir,crFn);
% crName = {[]};

% crPath = fullfile(crDir,crFn);
% Choose CR
% crName = {'CR1','CR2'};
% crName = {'CR1'};


params.k = 5; % Nearest Neighborscr
params.weightType = 'gaussian'; % Use Gaussian Weighting Kernel
params.sigma = 1; % meters
% Check for CR existance and if no CR make empty CR table
CR = utils.CRpower(crPath, crName, demData, sarData, geomData, params.lambda, params.k, params.weightType, params.sigma);
fprintf('RCS Calibration Computed. Elapsed time: %.2f seconds\n', toc);

% Optional CR Figure
isCRfig = 0;
if isCRfig
tmpCR = readtable(crPath);
figure();
hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
hold on;
plot(tmpCR.Easting./1000,tmpCR.Northing./1000,'^k','MarkerSize',5,'MarkerFaceColor','k')
text(tmpCR.Easting./1000+0.01,tmpCR.Northing./1000+0.01, tmpCR.Name, 'FontSize', 12, 'Color', 'm','FontName','Serif','FontWeight','bold');
daspect([1,1,1])
set(gca,'Ydir','normal', 'fontname','serif','fontweight','bold','fontsize',14)
title('Camas: 02/26 Corner Reflectors')
% title('White Clouds Preserve: 09/25 Corner Reflectors')
% title('Mores Creek Summit: 2025 Corner Reflectors')
xlabel('Easting (km)');ylabel('Northing (km)');
ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
end

%% ------------------- SINGLE LOOK IMAGE PROCESSING ----------------------
fprintf('SLC Calibration - Correction \n');
tic
params.correction = {'gamma0','pixelarea','range2'};
sarData = insar.process_single_look_images(sarData, demData, geomData, CR, params.lambda, params.correction, params.filterSize);
fprintf('SLC processing complete. Elapsed time: %.2f seconds\n', toc);

%% --------------------- INTERFEROGRAMS ---------------------
fprintf('Interferometric Processing \n');
tic
params.filterSize = 9;
unwrapOpts.method = 'multiseed'; %'multiseed', 'none'
unwrapOpts.qualityThresh = 0.3;
unwrapOpts.rampRemoval.enable = false;
unwrapOpts.sigma = 5;
unwrapOpts.pairingMode = 'allCombos';
[insarData, unwrapOpts] = insar.process_interferometric_phase(sarData, geomData, CR, params, unwrapOpts);
insarClosureData = insar.compute_insar_closure_bias(insarData);
fprintf('InSAR processing complete. Elapsed time: %.2f seconds\n', toc);

%% -------------------------- EXPORT DATA ---------------------------------
exportOpts = struct();

% Save mode
exportOpts.saveMode = 'standard';
exportOpts.useCompression = false;

% Export behavior
exportOpts.exportGeotiff = true;
exportOpts.overwrite = true;

% Data inclusion
exportOpts.includeComplex = true;      % helps verify complex → mag/phase exports
exportOpts.includePredictors = false;  % skip heavy predictor fields for now
exportOpts.includeLonLat = false;      % keep file size reasonable

% Naming linkage (IMPORTANT)
exportOpts.sarDataForNaming  = sarData;
exportOpts.geomDataForNaming = geomData;
exportOpts.runLabel = 'CAMAS';
% exportOpts.exportSubdir = 'CAMAS_validation_run01'; % Manual Dir Naming

% Figure Export
cmap = [[1,1,1];storeCmap];
exportOpts.exportFigures = true;
exportOpts.figureFormat = 'png';
exportOpts.figureDPI = 300;
exportOpts.figureVisible = 'on';   % 'on' if you want to inspect interactively
exportOpts.useKmAxes = true;
exportOpts.cmap = cmap;
exportOpts.showHillshade = true;

exportOpts.plotOpts.common = struct( ...
    'visible', 'on', ...
    'useKm', true, ...
    'showHillshade', true,...
    'alpha', 0.625);

exportOpts.plotOpts.sar.mli = struct( ...
    'cmap', cmap, ...
    'climVals', [],...
    'showHillshade', true, ...
    'titleText', 'Multi-Looked Intensity',...
    'cbarLabel', 'Power');

exportOpts.plotOpts.sar.db = struct( ...
    'cmap', cmap, ...
    'climVals', [-5 20],...
    'showHillshade', true,...
    'titleText', 'Multi-Looked Intensity',...
    'cbarLabel', 'Power (dB)');

exportOpts.plotOpts.insar.phzWrapped = struct( ...
    'cmap', cmap, ...
    'climVals', [-pi pi], ...
    'showHillshade', true,...
    'titleText', 'Wrapped Phase',...
    'cbarLabel', '\Delta \phi(rad)');

exportOpts.plotOpts.insar.phzReferenced = struct( ...
    'cmap', cmap, ...
    'climVals', [], ...
    'showHillshade', true,...
    'phasePiLocked', true,...
    'titleText', 'Referenced Phase',...
    'cbarLabel', '\Delta \phi(rad)');

exportOpts.plotOpts.insar.coherence = struct( ...
    'cmap', cmap, ...
    'climVals', [0 1],...
    'showHillshade', true,...
    'titleText', 'Coherence');

exportOpts.plotOpts.insar.penetration = struct( ...
    'cmap', cmap, ...
    'climVals', [-1 1], ...
    'showHillshade', true,...
    'titleText', 'Penetraion Depth',...
    'cbarLabel', '(m)');

exportOpts.plotOpts.insar.phaseNoiseStd = struct( ...
    'cmap', cmap, ...
    'climVals', [],...
    'showHillshade', true,...
    'titleText', 'Phase Noise',...
    'cbarLabel', '\sigma_{\Delta \phi}(rad)');

exportOpts.plotOpts.closure.closureMap = struct( ...
    'cmap', cmap, ...
    'climVals', [-pi pi], ...
    'showHillshade', true,...
    'titleText', 'Phase Closure',...
    'cbarLabel', '\Delta \phi Bias (rad)');

% Big Smoke
fprintf('Data Exporting \n');
tic
exportInfo = io.export_insar_project(outDir, demData, sarData, geomData, insarData, insarClosureData, exportOpts);
fprintf('Export complete. Elapsed time: %.2f seconds\n', toc);
%% Make Figures
% Subplot
db = sarData(1).db{1};
dbMask = db<-30;
dbIx = find(dbMask);
% coh = medfilt2(insarData(1).coherence,[11,11]);
% coh = medfilt2(cor,[11,11]);
coh = insarData(1).coherence;
cohMask = coh < 0.5;
cohIx = find(cohMask);
% coh = cor;
% cor = coh;
% Power
figure();
subplot(1,2,1)
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; hold on;
hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,db,'AlphaData',0.625);daspect([1,1,1]);colormap([[1 1 1];cmap]);hc=colorbar;
ylabel(hc,'Power (dB)','fontname','serif','fontweight','bold','fontsize',14)
% xlabel('Longitude');ylabel('Latitude');
xlabel('Easting (km)');ylabel('Northing (km)');
clim([quantile(db(:),[0.1,0.975])])
clim([-30 10])
% clim([-30 10])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
title('a) White Clouds Preserve: 09/10/25')

% title('White Clouds Preserve: 09/13/24')
% title('a)   Mores Creek Summit: 02/20/25')
% title('a)   Mores Creek Summit: 03/19/24')


ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')

% Interferogram
% phz = phz_equiv;
% phz = insarData(1).phzUnwrapped;
phz = insarData(1).phzReferenced;
% phz = insarData(1).phzWrapped;
% phz = phz_corr;
% phz = angle(insarData(1).complexCoherence);
phz(dbIx) = NaN;
phz(cohIx) = NaN;
% phz = insar.vector_median_filter(phz,5);

subplot(1,2,2)
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; hold on;
hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,phz,'AlphaData',0.625);daspect([1,1,1]);colormap([[1 1 1];cmap]);hc=colorbar;
ylabel(hc,'\Delta \phi (rad)','fontname','serif','fontweight','bold','fontsize',14)
xlabel('Easting (km)');ylabel('Northing (km)');
% xlabel('Longitude');ylabel('Latitude');
clim([-2.*pi 2.*pi])
% clim([quantile(deltaLWC(deltaLWC>0),[0.05,0.95])])
% clim([-3.*pi 3.*pi])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
title('b) White Clouds Preserve: 09/10/25')

% title('White Clouds Preserve: 09/13/24')
% title('b)   Mores Creek Summit: 02/20/25')
% title('b)   Mores Creek Summit: 03/19/24')


ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% set(gcf,'Position',[165	351.400000000000	1369.60000000000	420.000000000000])
set(gcf,'Position',[416.200000000000	185.800000000000	1099.20000000000	420.000000000000])
% exportgraphics(gcf,[sarData(1).dir,'\figures\WCP091324GLSARHHHV.png'],Resolution=300)
% exportgraphics(gcf,[sarData(2).dir,'\figures\MCS022025GLSARbistatic.png'],Resolution=300)
% exportgraphics(gcf,[sarData(1).dir,'\figures\WCP001025south120-120-rampRemoved.png'],Resolution=300)

% Coherence
figure();
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; hold on;
hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,coh,'AlphaData',0.625);daspect([1,1,1]);colormap([[1 1 1];cmap]);hc=colorbar;
ylabel(hc,'Coherence','fontname','serif','fontweight','bold','fontsize',14)
% xlabel('Longitude');ylabel('Latitude');
xlabel('Easting (km)');ylabel('Northing (km)');
% clim([quantile(coh(:),[0.1,0.975])])
clim([0 1])
% clim([-30 10])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
title('White Clouds Preserve: 09/10/25')

% title('White Clouds Preserve: 09/13/24')
% title('a)   Mores Creek Summit: 02/20/25')
% title('a)   Mores Creek Summit: 03/19/24')


ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% exportgraphics(gcf,[sarData(1).dir,'\figures\WCP001025south120-120-coherence.png'],Resolution=300)

% Signal Penetration
if insarData.penetrationValid
    penetration = insarData.penetration;
    % penetration(coh<0.3) = NaN;
    % penetration(penetration > 0) = NaN;
    % penetration
    % penetration = -exp(penetration);
    figure();
    imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
    utils.freezeColors; hold on;
    hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,penetration,'AlphaData',0.625);daspect([1,1,1]);colormap([[1 1 1];cmap]);hc=colorbar;
    ylabel(hc,'Signal Penetration (m)','fontname','serif','fontweight','bold','fontsize',14)
    % xlabel('Longitude');ylabel('Latitude');
    xlabel('Easting (km)');ylabel('Northing (km)');
    % clim([quantile(coh(:),[0.1,0.975])])
    clim([-1 0])
    % clim([-30 10])
    set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
    title('White Clouds Preserve: 09/12/25')

    % title('White Clouds Preserve: 09/13/24')
    % title('a)   Mores Creek Summit: 02/20/25')
    % title('a)   Mores Creek Summit: 03/19/24')


    ax = ancestor(hI, 'axes');
    ax.XAxis.Exponent = 0;
    xtickformat('%.1f')
    ax.YAxis.Exponent = 0;
    ytickformat('%.1f')
    % exportgraphics(gcf,['E:\WCP\0925campaign\GLSAR\South\20250912\lift_1_120_110\figures\WCP091225south120-110-penetration.png'],Resolution=300)
end

% Phase Screen
figure();
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(demData.aspect+45)+sind(demData.aspect+45))).*sind(2.5.*demData.slope));colormap(bone)
utils.freezeColors; hold on;
hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,insarData(1).phzReferenced - insarData(1).refScreen,'AlphaData',0.625);daspect([1,1,1]);colormap([[1 1 1];cmap]);hc=colorbar;
% ylabel(hc,'Phase Reference Screen (\Delta \phi)','fontname','serif','fontweight','bold','fontsize',14)
ylabel(hc,'Residual Phase (\Delta \phi)','fontname','serif','fontweight','bold','fontsize',14)

% xlabel('Longitude');ylabel('Latitude');
xlabel('Easting (km)');ylabel('Northing (km)');
% clim([quantile(db(:),[0.1,0.975])])
clim([-2.*pi 2.*pi])
% clim([-30 10])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)
title('a) White Clouds Preserve: 09/10/25')

% title('White Clouds Preserve: 09/13/24')
title('Mores Creek Summit: 02/23/26 - 02/26/26')
% title('a)   Mores Creek Summit: 03/19/24')


ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')
% 
% figure();
% for kk = 1:6
% subplot(2,3,kk)
% alpha_mask = ones(size(insarData(kk).phzUnwrapped)); % Initialize with opaque (1)
% alpha_mask(isnan(insarData(kk).phzUnwrapped)) = 0;
% % tmp = floor(insarData(kk).phzReferenced./(2.*pi));
% imagesc(insarData(kk).phzWrapped,'AlphaData',alpha_mask);
% daspect([1,1,1])
% title(['Pair '],num2str(insarData(kk).pair))
% colorbar; colormap(cmap)
% % clim([-2.*pi,2.*pi])
% % clim([-5,5])
% end

% Trajectory Figure
% figure();
% subplot(3,1,1)
% plot(sarData(1).traj{1}(:,1),'k','LineWidth',2);hold on; plot(sarData(2).traj{1}(:,1),'m','LineWidth',2)
% ylabel('Easting (m)')
% grid on; grid minor;
% subplot(3,1,2)
% plot(sarData(1).traj{1}(:,2),'k','LineWidth',2);hold on; plot(sarData(2).traj{1}(:,2),'m','LineWidth',2)
% ylabel('Northing (m)')
% grid on; grid minor;
% subplot(3,1,3)
% plot(sarData(1).traj{1}(:,3),'k','LineWidth',2);hold on; plot(sarData(2).traj{1}(:,3),'m','LineWidth',2)
% ylabel('Altitude (m)')
% xlabel('Distance (m)')
% grid on; grid minor;

% Trajectory Figure
% rmsd = sqrt((sarData(1).traj{1}(:,1)- sarData(2).traj{1}(:,1)).^2+(sarData(1).traj{1}(:,2)- sarData(2).traj{1}(:,2)).^2+(sarData(1).traj{1}(:,3)- sarData(2).traj{1}(:,3)).^2);
rmsd = sqrt((sarData(geomData.pairResults.ii1).traj{geomData.pairResults.jj1}(:,1)- sarData(geomData.pairResults.ii2).traj{geomData.pairResults.jj2}(:,1)).^2+(sarData(geomData.pairResults.ii1).traj{geomData.pairResults.jj1}(:,2)- sarData(geomData.pairResults.ii2).traj{geomData.pairResults.jj2}(:,2)).^2+(sarData(geomData.pairResults.ii1).traj{geomData.pairResults.jj1}(:,3)- sarData(geomData.pairResults.ii2).traj{geomData.pairResults.jj2}(:,3)).^2);
azAxis = sarData(geomData.pairResults.ii1).azimuthAxis{1};
figure();
subplot(4,1,1)
plot(azAxis,sarData(geomData.pairResults.ii1).traj{geomData.pairResults.jj1}(:,1)-sarData(geomData.pairResults.ii2).traj{geomData.pairResults.jj2}(:,1),'k','LineWidth',2);
ylabel('Easting Difference (m)')
% ylim([-1 1])
% ylim([mean(sarData(1).traj{1}(:,1)-sarData(1).traj{2}(:,1))-1,mean(sarData(1).traj{1}(:,1)-sarData(1).traj{2}(:,1))+1])
grid on; grid minor;
subplot(4,1,2)
plot(azAxis,sarData(geomData.pairResults.ii1).traj{geomData.pairResults.jj1}(:,2)-sarData(geomData.pairResults.ii2).traj{geomData.pairResults.jj2}(:,2),'k','LineWidth',2);
ylabel('Northing Difference (m)')
% ylim([-1 1])
% ylim([mean(sarData(1).traj{1}(:,2)-sarData(1).traj{2}(:,2))-1, mean(sarData(1).traj{1}(:,2)-sarData(1).traj{2}(:,2))+1])
grid on; grid minor;
subplot(4,1,3)
plot(azAxis,sarData(geomData.pairResults.ii1).traj{geomData.pairResults.jj1}(:,3)- sarData(geomData.pairResults.ii2).traj{geomData.pairResults.jj2}(:,3),'k','LineWidth',2);
ylabel('Altitude Difference (m)')
xlabel('Distance (m)')
grid on; grid minor;
subplot(4,1,4)
plot(azAxis,rmsd,'k','LineWidth',2);
ylabel('Mean Difference (m)')
xlabel('Distance (m)')
grid on; grid minor;


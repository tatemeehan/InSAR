function InSAR_Planning_GUI()

%% =====================================================
% Default values
%% =====================================================

defaults.traj = '/Users/rdcrlcjw/Desktop/InSAR-main/Trajectories/Alpine_Garden.txt';
defaults.dem  = '/Users/rdcrlcjw/Desktop/InSAR-main/DEMs/mtwash_DEM_3.tif';
defaults.cr   = '/Users/rdcrlcjw/Desktop/InSAR-main/CRs/MCS_CRs_2026.csv';
defaults.cmap = '/Users/rdcrlcjw/Desktop/InSAR-main/colormaps/RdYlBu.csv';
defaults.out  = '/Users/rdcrlcjw/Desktop/InSAR-main/Export_Figures';

defaults.c = 0.3;
defaults.f = 1.3;

defaults.AGL = 100;
defaults.verticalShift = 20;
defaults.shiftDistance = 0;
defaults.dx = 20;

%% =====================================================
% GUI Window
%% =====================================================

fig = uifigure( ...
    'Name','InSAR Trajectory Planning Tool', ...
    'Position',[100 100 700 550]);

%% Title
uilabel(fig,...
    'Text','InSAR Geometry & Trajectory Planner',...
    'FontSize',20,...
    'FontWeight','bold',...
    'Position',[180 510 400 30]);

%% =====================================================
% FILE PATH PANEL
%% =====================================================

pathPanel = uipanel(fig,...
    'Title','Input / Output Paths',...
    'Position',[20 280 660 210]);

trajField = addPathField(pathPanel,'Trajectory File',defaults.traj,160);
demField  = addPathField(pathPanel,'DEM File',defaults.dem,120);
crField   = addPathField(pathPanel,'Corner Reflector File',defaults.cr,80);
cmapField = addPathField(pathPanel,'Colormap File',defaults.cmap,40);
outField  = addPathField(pathPanel,'Output Directory',defaults.out,0,true);

%% =====================================================
% RADAR PANEL
%% =====================================================

radarPanel = uipanel(fig,...
    'Title','Radar Parameters',...
    'Position',[20 180 320 90]);

cField = addNumField(radarPanel,'Wave Speed (m/ns)',defaults.c,35);
fField = addNumField(radarPanel,'Frequency (GHz)',defaults.f,5);

%% =====================================================
% FLIGHT GEOMETRY PANEL
%% =====================================================

flightPanel = uipanel(fig,...
    'Title','Flight Geometry',...
    'Position',[360 180 320 90]);

AGLField    = addNumField(flightPanel,'Flight AGL (m)',defaults.AGL,35);
vShiftField = addNumField(flightPanel,'Vertical Track Shift (m)',defaults.verticalShift,5);

shiftPanel = uipanel(fig,...
    'Title','Trajectory Controls',...
    'Position',[20 80 660 90]);

hShiftField = addNumField(shiftPanel,'Horizontal Track Shift (m)',defaults.shiftDistance,30);
dxField     = addNumField(shiftPanel,'Interpolation Spacing (m)',defaults.dx,0);

%% =====================================================
% RUN BUTTON
%% =====================================================

runBtn = uibutton(fig,...
    'Text','Run InSAR Planning',...
    'FontSize',16,...
    'FontWeight','bold',...
    'Position',[260 25 180 45],...
    'ButtonPushedFcn',@runPlanning);

%% =====================================================
% RUN CALLBACK
%% =====================================================

function runPlanning(~,~)

paths.traj = trajField.Value;
paths.dem  = demField.Value;
paths.cr   = crField.Value;
paths.cmap = cmapField.Value;
paths.out  = outField.Value;

params.c = cField.Value;
params.f = fField.Value;

AGL = AGLField.Value;
verticalShift = vShiftField.Value;
shiftDistance = hShiftField.Value;
dx = dxField.Value;

disp("Running InSAR Planning...")

InSAR_Planning(paths, params, AGL, verticalShift, shiftDistance, dx);

end

end

function field = addPathField(parent,label,default,y,isFolder)

if nargin < 5
    isFolder = false;
end

uilabel(parent,'Text',label,'Position',[10 y+5 160 22]);

field = uieditfield(parent,'text',...
    'Position',[170 y+5 380 22],...
    'Value',default);

uibutton(parent,'Text','Browse',...
    'Position',[560 y+5 80 22],...
    'ButtonPushedFcn',@(btn,event)browse);

    function browse
        if isFolder
            p = uigetdir;
        else
            [f,pth] = uigetfile('*.*');
            p = fullfile(pth,f);
        end
        
        if p~=0
            field.Value = p;
        end
    end

end

function field = addNumField(parent,label,default,y)

uilabel(parent,'Text',label,'Position',[10 y+5 180 22]);

field = uieditfield(parent,'numeric',...
    'Position',[200 y+5 80 22],...
    'Value',default);

end
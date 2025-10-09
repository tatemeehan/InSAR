function setup_insar_toolbox(basePath)
%SETUP_INSAR_TOOLBOX Add InSAR DEM Toolbox folders to MATLAB path
%
%   setup_insar_toolbox()         % Adds current folder
%   setup_insar_toolbox(path)     % Adds specified toolbox path

    if nargin < 1
        basePath = fileparts(mfilename('fullpath'));
    end

    addpath(genpath(basePath));
    disp(['[InSAR Toolbox] Path added: ', basePath]);
end

function [traj, v, t, s, L] = load_geodetic_trajectory(geodPath,dx)
%LOAD_GEODETIC_TRAJECTORY Load geodetic trajectory from pos.dat.geod file
%
%   traj = load_geodetic_trajectory(geodPath)
%
%   Input:
%       geodPath - full path to the pos.dat.geod file (lat, lon, height)
%
%   Output:
%       traj - Nx3 matrix of [X_utm, Y_utm, Z_elev] in meters

    % Check if file exists
    if ~isfile(geodPath)
        error('Trajectory file not found: %s', geodPath);
    end
    aziTimePath = [geodPath(1:end-4),'azi_times'];
    if ~isfile(aziTimePath)
        warning('Azimuth times not supplied')
    else
        time = load(aziTimePath);
    end


    % Load [lat (deg), lon (deg), height (m)]
    geod = load(geodPath);
    if size(geod,2) < 3
        error('Trajectory file must have at least 3 columns: lat, lon, height');
    end

    % Convert to UTM (Zone inferred from lat/lon)
    lat = geod(:,1);
    lon = geod(:,2);
    alt = geod(:,3);

    [x, y, zone] = utils.deg2utm(lat, lon);

    traj = [x, y, alt];
    
    % Resample Trajectory
    if ~isfile(aziTimePath)
        [traj, L] = utils.resample_trajectory(traj,dx);
        t = []; s = []; v = [];
    else
        % Compute Velocity
        [traj, t, s, L, v] = utils.resample_trajectory_with_time(traj, time, dx);
    end

    % Optional: store metadata
    % trajInfo.zone = zone;
    % trajInfo.hemisphere = hemisphere;
    % assignin('base', 'trajInfo', trajInfo); % optional debug visibility
end

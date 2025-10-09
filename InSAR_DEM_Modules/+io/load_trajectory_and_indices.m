function [r1, r2] = load_trajectory_and_indices(trajFile, dirPath, fn1_base, fn2_base, dx, isAlignTrajectories)
%LOAD_TRAJECTORY_AND_INDICES Match and align SAR trajectory segments
%
%   [r1, r2] = load_trajectory_and_indices(trajFile, dirPath, fn1_base, fn2_base, dx)
%
%   Inputs:
%       trajFile - path to trajectory .txt file
%       dirPath  - folder where pos.dat.azi_times / az_ind are stored
%       fn1_base, fn2_base - base filenames (without extension) of SLCs
%       dx       - desired spacing (m)
% isAlignTrajectories - (Optional) Default = true;
%
%   Outputs:
%       r1, r2 - aligned trajectory points (UTM)

if nargin < 6
    isAlignTrajectories = true;
end

    traj = csvread(fullfile(dirPath, trajFile));
    traj(:,2:3) = traj(:,2:3) * 180/pi;  % radians to degrees

    azi1 = csvread(fullfile(dirPath, [fn1_base, 'pos.dat.azi_times']));
    ix1  = csvread(fullfile(dirPath, [fn1_base, 'pos.az_ind']));
    azi2 = csvread(fullfile(dirPath, [fn2_base, 'pos.dat.azi_times']));
    ix2  = csvread(fullfile(dirPath, [fn2_base, 'pos.az_ind']));

    [~, i1] = min(abs(azi1(ix1)' - traj(:,1))); i1 = i1(1):i1(2);
    [~, i2] = min(abs(azi2(ix2)' - traj(:,1))); i2 = i2(1):i2(2);

    r1 = traj(i1,:);
    r2 = traj(i2,:);
    if isAlignTrajectories
        % Align Trajetories and Maintain Vertical Offset
        [r1, r2] = utils.alignTrajectoriesXY(r1, r2, dx);
    end
end
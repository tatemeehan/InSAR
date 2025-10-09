function [r1, r2] = alignTrajectoriesXY(traj1, traj2, dx)
% alignTrajectoriesXY Align two trajectories in XY and resample at 1 m spacing
% 
% Inputs:
%   traj1 - Nx4 array [time, lat, lon, alt] for reference trajectory (e.g. lower altitude)
%   traj2 - Mx4 array [time, lat, lon, alt] for comparison trajectory
%   dx    - desired spacing along-track in meters (default: 1)
%
% Outputs:
%   r1 - Px3 array [x y z] for reference trajectory (uniform spacing)
%   r2 - Px3 array [x y z] for comparison trajectory, aligned in XY

if nargin < 3
    dx = 1;
end

% Convert to UTM
[trajX1, trajY1] = utils.deg2utm(traj1(:,2), traj1(:,3));
trajZ1 = traj1(:,4);

[trajX2, trajY2] = utils.deg2utm(traj2(:,2), traj2(:,3));
trajZ2 = traj2(:,4);

% Interpolate traj1 (reference) at uniform spacing
[arc1,~] = utils.arclength(trajX1, trajY1, trajZ1, 'pchip');
s_uniform = [0:dx:arc1]';  % new utils.arclength steps
r1 = utils.interparc(s_uniform / arc1, trajX1, trajY1, trajZ1, 'pchip');  % reference

% Interpolate traj2 to match XY of r1
[~,seg2] = utils.arclength(trajX2, trajY2, trajZ2, 'pchip');
s2 = [0;cumsum(seg2)];
% Remove duplicates in s2 to prevent interp1 failure
[s2_unique, uniqueIx2] = unique(s2);
z2_unique = trajZ2(uniqueIx2);

% Interpolate z2 onto s1
z2_interp = interp1(s2_unique, z2_unique, s_uniform, 'pchip', 'extrap');

r2 = [r1(:,1:2), z2_interp];

end

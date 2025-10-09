function [traj_out, L] = resample_trajectory(traj, dx)
    trajX = traj(:,1); trajY = traj(:,2); trajZ = traj(:,3);
    [~,xIx] = unique(trajX);[~,yIx] = unique(trajY);
    ix = intersect(xIx,yIx);
    trajX = trajX(ix); trajY = trajY(ix); trajZ = trajZ(ix);
    [L, ~] = utils.arclength(trajX, trajY, trajZ, 'pchip');
    s = (0:dx:L)';
    traj_out = utils.interparc(s / L, trajX, trajY, trajZ, 'pchip');
end
function [traj_out, t_out, s_out, L, V_out] = resample_trajectory_with_time(trajXYZ, t_in, dx)
% RESAMPLE_TRAJECTORY_WITH_TIME
%   Resample a trajectory to uniform along-track spacing dx (m), carrying time.
%   Handles seconds-of-week wrap and duplicates; returns velocity dX/dt too.
%
% Inputs:
%   trajXYZ : [N x 3] positions (x,y,z) in meters (UTM/ENU)
%   t_in    : [N x 1] time stamps (e.g., seconds-of-week, possibly wrapping)
%   dx      : scalar along-track spacing in meters
%
% Outputs:
%   traj_out : [M x 3] resampled (x,y,z) at s = 0:dx:L
%   t_out    : [M x 1] time interpolated at the same s samples
%   s_out    : [M x 1] along-track distances (0...L)
%   L        : total arclength (m)
%   V_out    : [M x 3] velocity (m/s) computed as dX/dt on (traj_out,t_out)

    arguments
        trajXYZ (:,3) double
        t_in    (:,1) double
        dx      (1,1) double {mustBePositive}
    end

    % --- 0) Clean duplicates (keep order), drop NaNs
    good = all(isfinite(trajXYZ),2) & isfinite(t_in);
    trajXYZ = trajXYZ(good,:);  t_in = t_in(good);

    % Unique consecutive points (avoid zero arclength steps)
    dxyz = [true; any(abs(diff(trajXYZ))>0,2)];
    trajXYZ = trajXYZ(dxyz,:);  t_in = t_in(dxyz);

    % --- 1) Unwrap seconds-of-week (detect wrap ~ 604800 s)
    sow = 7*24*3600;  % 604800
    t = t_in;
    for k = 2:numel(t)
        if t(k) < t(k-1) - sow/2
            % wrapped forward (e.g., 604799 -> 0)
            t(k:end) = t(k:end) + sow;
        elseif t(k) > t(k-1) + sow/2
            % extremely rare: wrapped backward
            t(k:end) = t(k:end) - sow;
        end
    end

    % --- 2) Compute arclength s along the original track
    % If you already have utils.arclength, feel free to use it; here is inline:
    % ds = sqrt(sum(diff(trajXYZ,1,1).^2, 2));
    % s  = [0; cumsum(ds)];
    % L  = s(end);
    [L, s] = utils.arclength(trajXYZ(:,1), trajXYZ(:,2), trajXYZ(:,3), 'pchip');
    s  = [0; cumsum(s)];
    % s = (0:dx:L)';

    if L == 0
        % Degenerate trajectory
        traj_out = trajXYZ(1,:); t_out = t(1); s_out = 0; V_out = [0 0 0];
        return;
    end

    % Some users’ tracks can stall briefly; ensure s is strictly increasing for interp1:
    [s, uniqIdx] = unique(s, 'stable');
    trajXYZ = trajXYZ(uniqIdx,:); t = t(uniqIdx);

    % --- 3) Target arclength samples
    s_out = (0:dx:L).';  % column vector
    u = s_out / L;       % normalized arclength in [0,1]

    % --- 4) Interpolate positions on uniform s (uses your interparc-like helper)
    % If you have utils.interparc(u, x, y, z, 'pchip'):
    traj_out = utils.interparc(u, trajXYZ(:,1), trajXYZ(:,2), trajXYZ(:,3), 'pchip');

    % --- 5) Interpolate time vs arclength (monotonic in s; pchip avoids overshoot)
    t_out = interp1(s, t, s_out, 'pchip', 'extrap');

    % --- 6) Velocity via time-derivative on the resampled track
    % Use gradient with non-uniform t_out (robust and simple)
    V_out = zeros(size(traj_out));
    dt = diff(t_out);
    if any(dt <= 0)
        % If something weird happened with time, fall back to uniform finite differencing
        warning('Non-monotonic t_out detected; velocities may be unreliable.');
        dt_mean = median(diff(t_out));
        V_out(:,1) = gradient(traj_out(:,1), dt_mean);
        V_out(:,2) = gradient(traj_out(:,2), dt_mean);
        V_out(:,3) = gradient(traj_out(:,3), dt_mean);
    else
        V_out(:,1) = gradient(traj_out(:,1), t_out);
        V_out(:,2) = gradient(traj_out(:,2), t_out);
        V_out(:,3) = gradient(traj_out(:,3), t_out);
    end
end

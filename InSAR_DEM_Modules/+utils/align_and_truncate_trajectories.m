% utils/align_and_truncate_trajectories.m
function [sarData] = align_and_truncate_trajectories(sarData, dx)
%ALIGN_AND_TRUNCATE_TRAJECTORIES Align and truncate all trajectories in sarData
%
%   [sarData, startCoord, endCoord] = align_and_truncate_trajectories(sarData, dx)
%
%   Updates sarData in place with truncated and aligned trajectories.

    trajIn = {};
    trajMap = [];  % Mapping from trajIn index to (ii,jj)

    % Collect all resampled trajectories
    for ii = 1:numel(sarData)
        for jj = 1:numel(sarData(ii).slcFiles)
            trajIn{end+1} = sarData(ii).traj{jj};
            trajMap(end+1,:) = [ii, jj];
        end
    end

    % Truncate to shared segment
    [trajOut, startCoord, endCoord, indices] = utils.truncate_trajectories_global(trajIn);
    N = (size(trajOut{1}, 1)-1) * dx;

    % Reassign to sarData structure
    for kk = 1:numel(trajOut)
        ii = trajMap(kk,1); jj = trajMap(kk,2);
        sarData(ii).traj{jj} = trajOut{kk};
        if ~isempty(sarData(ii).velocity{jj})%isfield(sarData,'velocity')
        sarData(ii).velocity{jj} = sarData(ii).velocity{jj}([indices(kk,1):indices(kk,2)],:);
        sarData(ii).slowTime{jj} = sarData(ii).slowTime{jj}([indices(kk,1):indices(kk,2)],:);
        sarData(ii).azimuthAxis{jj} = sarData(ii).azimuthAxis{jj}([indices(kk,1):indices(kk,2)],:)...
            -sarData(ii).azimuthAxis{jj}([indices(kk,1)],:);
        end
        sarData(ii).trajLength{jj}  = N;
        sarData(ii).trajStart{jj} = startCoord;
        sarData(ii).trajEnd{jj} = endCoord;
        % sarData(ii).truncateStart{jj} = indices(k,1);
        % sarData(ii).truncateEnd{jj}   = indices(k,2);
    end
end

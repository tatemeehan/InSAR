function sarData = compute_doppler_centroids_from_icmap(sarData, geomData, demData, lambda, opts)
% Compute per-SLC Doppler centroid maps using geometry closestIndex maps.
%
% Inputs:
%   sarData(ii).traj{jj}     [N x 3] UTM/ENU positions
%   sarData(ii).velocity{jj} [N x 3] UTM/ENU velocities
%   geomData.sarGeometry{idx}.closestIndex_master / _slave
%   geomData.sarGeometry{idx}.lookMask
%   demData.X, demData.Y, demData.dem grids (same size)
%   lambda  radar wavelength [m]
%   opts:
%     .force         (false) recompute even if cached
%     .windowSamples (51)    half-window in samples (±W)
%     .windowMeters  []      half-window in meters (overrides samples)
%     .windowTime    []      half-window in seconds (needs opts.t)
%     .useR2Weight   (true)  apply 1/R^2 weighting inside window
%     .smoothVel     (51)    Gaussian length for along-track V smoothing (0 = off)
%     .t             {}      cell array of slow time per (ii).slowTime{jj} if using windowTime
%
% Output:
%   sarData(ii).doppler.fdc{jj}  (Hz) same size as DEM

    arguments
        sarData
        geomData
        demData
        lambda (1,1) double
        opts.force logical = false
        opts.windowSamples (1,1) double = 51
        opts.windowMeters  = []
        opts.windowTime    = []
        opts.useR2Weight logical = true
        opts.smoothVel (1,1) double = 51
        opts.t = []
    end

    X = demData.X; Y = demData.Y; Z = demData.dem;
    nDir = numel(sarData);

    for ii = 1:nDir
        nBurst = numel(sarData(ii).slc);
        if ~isfield(sarData(ii),'doppler') || ~isfield(sarData(ii).doppler,'fdc')
            sarData(ii).doppler.fdc = cell(1,nBurst);
        end
        for jj = 1:nBurst
            % Skip if cached
            if ~opts.force && ~isempty(sarData(ii).doppler.fdc{jj}) ...
                    && all(isfinite(sarData(ii).doppler.fdc{jj}(:)))
                continue
            end

            % Geometry index & pass
            gi = geomData.slcGeomIndex(ii,jj);
            if isnan(gi.idx)
                warning('Missing geometry idx for (%d,%d); skipping Doppler.', ii, jj);
                continue
            end
            G = geomData.sarGeometry{gi.idx};
            if isfield(G,'lookMask') && ~isempty(G.lookMask)
                lookMask = logical(G.lookMask);
            else
                lookMask = true(size(X));
            end

            % Choose the correct closestIndex map for this pass
            switch lower(gi.pass)
                case 'master'
                    if ~isfield(G,'closestIndex_master')
                        warning('No closestIndex_master in geom idx %d; skipping (%d,%d).', gi.idx, ii, jj);
                        continue
                    end
                    ic = G.closestIndex_master;
                case 'slave'
                    if ~isfield(G,'closestIndex_slave')
                        warning('No closestIndex_slave in geom idx %d; skipping (%d,%d).', gi.idx, ii, jj);
                        continue
                    end
                    ic = G.closestIndex_slave;
                otherwise
                    ic = G.closestIndex_master;
            end

            % Traj & velocity
            P = sarData(ii).traj{jj};
            V = sarData(ii).velocity{jj};
            if size(P,1) ~= size(V,1)
                error('P and V length mismatch for (%d,%d).', ii, jj);
            end

            % Optional: pre-smooth velocity along track
            if opts.smoothVel > 0
                V = smoothdata(V, 1, 'gaussian', max(3, round(opts.smoothVel)));
            end

            % Window options
            dopOpts = struct('windowSamples', opts.windowSamples, ...
                             'windowMeters',  opts.windowMeters, ...
                             'windowTime',    opts.windowTime, ...
                             'useR2Weight',   opts.useR2Weight);
            if ~isempty(opts.t) && numel(sarData(ii).slowTime) >= jj && ~isempty(sarData(ii).slowTime{jj})
                dopOpts.t = sarData(ii).slowTime{jj};
            end

            % Compute FDC (Hz)
            fdc = insar.doppler_from_icmap_weighted(P, V, lambda, X, Y, Z, ic, lookMask, dopOpts);
            sarData(ii).doppler.fdc{jj} = fdc;
        end
    end
end

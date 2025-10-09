function out = doppler_delta_vs_incidence_utm(X,Y,Z, incDeg, lambda, ...
                                              trajA_pos, trajA_vel, ...
                                              trajB_pos, trajB_vel, opts)
% Returns struct with tile-wise Δf_dc and incidence bins

if nargin<10, opts=struct(); end
if ~isfield(opts,'downsample'), opts.downsample=[4 4]; end
if ~isfield(opts,'bins'), opts.bins = 0:2:90; end

fdA = doppler_from_trajectory_map(X,Y,Z, trajA_pos, trajA_vel, lambda, opts);
fdB = doppler_from_trajectory_map(X,Y,Z, trajB_pos, trajB_vel, lambda, opts);
df  = fdB - fdA;

% Collect valid samples (where we computed fd)
mask = isfinite(df) & isfinite(incDeg);
inc = incDeg(mask); dfdc = df(mask);

% Bin statistics
[~,~,binID] = histcounts(inc, opts.bins);
binCtr = 0.5*(opts.bins(1:end-1)+opts.bins(2:end));
binMean = accumarray(max(binID,1), dfdc, [], @mean, NaN);
binMed  = accumarray(max(binID,1), dfdc, [], @median, NaN);
binCnt  = accumarray(max(binID,1), 1, [], @sum, 0);

out.fdA = fdA; out.fdB = fdB; out.dfdc = df;
out.tbl = table(inc, dfdc, 'VariableNames', {'inc_deg','dfdc_Hz'});
out.bin = table(binCtr(:), binMean, binMed, binCnt, ...
    'VariableNames', {'inc_deg','dfdc_mean','dfdc_median','count'});
end

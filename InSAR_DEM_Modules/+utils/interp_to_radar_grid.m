function dataInterp = interp_to_radar_grid(data, Rdata, targetX, targetY)
%INTERP_TO_RADAR_GRID Interpolate georeferenced raster to radar grid
%
%   dataInterp = interp_to_radar_grid(data, Rdata, targetX, targetY)

    dataInterp = mapinterp(data, Rdata, targetX, targetY);
end
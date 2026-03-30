function [surfaceNormal, aspect, slope] = compute_surface_normals(dem, lat, lon, EPSG)
%COMPUTE_SURFACE_NORMALS Compute terrain normal vectors from DEM
%
%   [surfaceNormal, aspect, slope] = compute_surface_normals(dem, lat, lon, EPSG)

    % Work in double for numerical stability during gradient computation
    demD = double(dem);
    latD = double(lat);
    lonD = double(lon);

    latlim = [min(latD(:)), max(latD(:))];
    lonlim = [min(lonD(:)), max(lonD(:))];
    sizeDEM = size(demD);

    gref = georefpostings(latlim, lonlim, sizeDEM, ...
        'RowsStartFrom', 'west', ...
        'ColumnsStartFrom', 'north');

    gref.GeographicCRS = projcrs(EPSG).GeographicCRS;

    [aspectD, slopeD, gy, gx] = gradientm(demD, gref);

    normZ = ones(size(gx), 'double');
    surfaceNormalD = normalize([-gx(:), -gy(:), normZ(:)], 2, 'norm');

    % Store compactly
    surfaceNormal = single(surfaceNormalD);
    aspect = single(aspectD);
    slope = single(slopeD);
end
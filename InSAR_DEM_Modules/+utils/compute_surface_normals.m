function [surfaceNormal,aspect,slope] = compute_surface_normals(dem, lat, lon, EPSG)
%COMPUTE_SURFACE_NORMALS Compute terrain normal vectors from DEM
%
%   surfaceNormal = compute_surface_normals(dem, lat, lon, EPSG)

    latlim = [min(lat(:)), max(lat(:))];
    lonlim = [min(lon(:)), max(lon(:))];
    sizeDEM = size(dem);

    gref = georefpostings(latlim, lonlim, sizeDEM, ...
        'RowsStartFrom','west', 'ColumnsStartFrom','north');
    gref.GeographicCRS = projcrs(EPSG).GeographicCRS;

    [aspect,slope,gy,gx] = gradientm(dem, gref);
    normZ = ones(size(gx));
    surfaceNormal = normalize([-gx(:), -gy(:), normZ(:)], 2, 'norm');
end
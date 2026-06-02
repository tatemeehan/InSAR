function write_tif(A, fn, Coords)
%WRITE_TIF Write projected GeoTIFF with consistent coordinate metadata.

if isfield(Coords, 'EPSG') && ~isempty(Coords.EPSG)
    geotiffwrite(fn, A, Coords.R, 'CoordRefSysCode', Coords.EPSG);
else
    geotiffwrite(fn, A, Coords.R);
end

fprintf('  wrote %s\n', fn);

end
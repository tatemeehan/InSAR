function [A, R, info] = readLidarTif_big(fullfilename, debug)
%READLIDARTIF_BIG Robust reader for tiled BigTIFF/GeoTIFF DEM rasters.
%
% This function:
%   1) Tries READGEORASTER first
%   2) Falls back to low-level TIFF tile reading for BigTIFF/tiled files
%   3) Returns raster A as SINGLE
%   4) Reconstructs a MapPostingsReference when possible
%
% Inputs
%   fullfilename : full path to GeoTIFF
%   debug        : optional logical flag for diagnostic printing (default false)
%
% Outputs
%   A    : raster data, class single
%   R    : map.rasterref.MapPostingsReference when possible
%   info : struct of diagnostic information

    if nargin < 2 || isempty(debug)
        debug = false;
    end

    info = struct();
    info.filename = fullfilename;
    info.debug = debug;

    % -------------------------------------------------------------
    % Metadata
    % -------------------------------------------------------------
    tifInfo = imfinfo(fullfilename);
    info.imfinfo = tifInfo;

    try
        gt = geotiffinfo(fullfilename);
        info.geotiffinfo = gt;
    catch ME
        gt = [];
        info.geotiffinfo_error = ME.message;
    end

    % -------------------------------------------------------------
    % Fast path
    % -------------------------------------------------------------
    try
        [A, R0] = readgeoraster(fullfilename, 'OutputType', 'single');
        A = applyNoData(A, tifInfo, gt);
        R = forcePostingsReference(R0, size(A), tifInfo, gt);
        info.method = 'readgeoraster';
        return
    catch ME
        info.readgeoraster_error = ME.message;
    end

    % -------------------------------------------------------------
    % Low-level tiled TIFF fallback
    % -------------------------------------------------------------
    t = Tiff(fullfilename, 'r');
    cleaner = onCleanup(@() close(t)); %#ok<NASGU>

    width  = tifInfo.Width;
    height = tifInfo.Height;

    if ~isfield(tifInfo, 'TileWidth') || ~isfield(tifInfo, 'TileLength')
        error('This fallback expects a tiled TIFF, but TileWidth/TileLength were not found.');
    end

    tileWidth  = tifInfo.TileWidth;
    tileLength = tifInfo.TileLength;

    numTilesAcross = ceil(width  / tileWidth);
    numTilesDown   = ceil(height / tileLength);

    info.tileWidth = tileWidth;
    info.tileLength = tileLength;
    info.numTilesAcross = numTilesAcross;
    info.numTilesDown = numTilesDown;

    if debug
        fprintf('Reading tiled TIFF fallback:\n');
        fprintf('  size(A) = [%d %d]\n', height, width);
        fprintf('  tile size = [%d %d]\n', tileLength, tileWidth);
        fprintf('  tiles = [%d down x %d across]\n', numTilesDown, numTilesAcross);
    end

    A = zeros(height, width, 'single');

    % Suppress warning backtraces only
    warnState = warning('off', 'backtrace');
    warnCleaner = onCleanup(@() warning(warnState)); %#ok<NASGU>

    warningCount = 0;
    transposedTileCount = 0;

    for rowTile = 1:numTilesDown
        row0 = (rowTile - 1) * tileLength + 1;
        rr = min(tileLength, height - row0 + 1);

        for colTile = 1:numTilesAcross
            col0 = (colTile - 1) * tileWidth + 1;
            cc = min(tileWidth, width - col0 + 1);

            % MATLAB tile numbering for this workflow: 1-based row-major
            tileNumber = (rowTile - 1) * numTilesAcross + colTile;

            try
                txt = evalc('tileData = readEncodedTile(t, tileNumber);');
            catch ME
                error(['Failed reading tile.\n' ...
                       '  rowTile=%d colTile=%d tileNumber=%d\n' ...
                       '  row0=%d col0=%d rr=%d cc=%d\n' ...
                       '  Error: %s'], ...
                       rowTile, colTile, tileNumber, row0, col0, rr, cc, ME.message);
            end

            if ~isempty(strtrim(txt))
                if contains(txt, 'Using code not yet in table')
                    warningCount = warningCount + 1;
                elseif debug
                    fprintf('Tile %d emitted TIFF warning text:\n%s\n', tileNumber, txt);
                end
            end

            tileData = single(tileData);
            [tr, tc] = size(tileData);

            if debug && rowTile <= 2 && colTile <= 3
                fprintf(['Tile (%d,%d): tileNumber=%d, tileData size=[%d %d], ' ...
                         'target size=[%d %d]\n'], ...
                         rowTile, colTile, tileNumber, tr, tc, rr, cc);
            end

            % Most common case
            if tr == rr && tc == cc
                A(row0:row0+rr-1, col0:col0+cc-1) = tileData;

            % Transposed case
            elseif tr == cc && tc == rr
                A(row0:row0+rr-1, col0:col0+cc-1) = tileData.';
                transposedTileCount = transposedTileCount + 1;

            % Full interior tile that needs cropping
            elseif tr >= rr && tc >= cc
                A(row0:row0+rr-1, col0:col0+cc-1) = tileData(1:rr, 1:cc);

            % Cropped transposed edge tile
            elseif tr >= cc && tc >= rr
                tmp = tileData.';
                A(row0:row0+rr-1, col0:col0+cc-1) = tmp(1:rr, 1:cc);
                transposedTileCount = transposedTileCount + 1;

            else
                error(['Unexpected tile shape.\n' ...
                       '  rowTile=%d colTile=%d tileNumber=%d\n' ...
                       '  tileData size=[%d %d]\n' ...
                       '  target size=[%d %d]'], ...
                       rowTile, colTile, tileNumber, tr, tc, rr, cc);
            end
        end
    end

    A = applyNoData(A, tifInfo, gt);
    R = buildPostingsReference(tifInfo, gt, size(A));

    info.method = 'Tiff.readEncodedTile';
    info.tiffWarningCount = warningCount;
    info.transposedTileCount = transposedTileCount;
end

% =====================================================================
function A = applyNoData(A, tifInfo, gt)
%APPLYNODATA Convert NoData value to NaN if present.

    ndv = [];

    if isfield(tifInfo, 'GDAL_NODATA') && ~isempty(tifInfo.GDAL_NODATA)
        ndv = str2double(tifInfo.GDAL_NODATA);
    elseif ~isempty(gt) && isfield(gt, 'GeoTIFFTags') ...
            && isfield(gt.GeoTIFFTags, 'GDAL_NODATA')
        ndv = str2double(gt.GeoTIFFTags.GDAL_NODATA);
    end

    if ~isempty(ndv) && ~isnan(ndv)
        A(A == ndv) = NaN;
    end
end

% =====================================================================
function R = forcePostingsReference(R0, rasterSize, tifInfo, gt)
%FORCEPOSTINGSREFERENCE Return MapPostingsReference when possible.

    if isa(R0, 'map.rasterref.MapPostingsReference')
        R = R0;
    else
        R = buildPostingsReference(tifInfo, gt, rasterSize);
    end
end

% =====================================================================
function R = buildPostingsReference(tifInfo, gt, rasterSize)
%BUILDPOSTINGSREFERENCE Construct MapPostingsReference.

    % Preferred path: geotiffinfo RefMatrix
    if ~isempty(gt) && isfield(gt, 'RefMatrix') && ~isempty(gt.RefMatrix)
        R = refmatToMapRasterReference(gt.RefMatrix, rasterSize, 'postings');
        return
    end

    % Fallback from TIFF tags
    if ~isfield(tifInfo, 'ModelPixelScaleTag') || ~isfield(tifInfo, 'ModelTiepointTag')
        error('Unable to build map reference: required GeoTIFF tags are missing.');
    end

    s   = tifInfo.ModelPixelScaleTag;
    tie = tifInfo.ModelTiepointTag;

    dx = s(1);
    dy = s(2);

    xul = tie(4);
    yul = tie(5);

    nRows = rasterSize(1);
    nCols = rasterSize(2);

    xlimits = [xul, xul + dx*(nCols - 1)];
    ylimits = [yul - dy*(nRows - 1), yul];

    R = maprefpostings(xlimits, ylimits, rasterSize);
end
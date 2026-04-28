function [A,R,X,Y,lon,lat,utmX,utmY,epsgCode] = readLidarTif(fullfilename)
%READLIDARTIF Robust LiDAR GeoTIFF reader with BigTIFF fallback support.
%
% Outputs:
%   A        - raster, single
%   R        - spatial referencing object
%   X,Y      - projected/world grid if projected CRS, otherwise lon/lat grid
%   lon,lat  - geographic coordinates if available
%   utmX,utmY- projected coordinates in EPSG zone if available
%   epsgCode - EPSG code if available, otherwise 0

    % Output suppression logic
    isTilde  = utils.detectOutputSuppression(nargout);
    isOutput = ~isTilde;

    % Defaults
    X = 0; Y = 0;
    lon = 0; lat = 0;
    utmX = 0; utmY = 0;
    epsgCode = 0;

    % -------------------------------------------------------------
    % Read raster
    % -------------------------------------------------------------
    usedFallback = false;
    try
        [A,R] = readgeoraster(fullfilename, 'OutputType', 'single');
    catch
        [A,R,readInfo] = io.readLidarTif_big(fullfilename,false);
        usedFallback = true; %#ok<NASGU>
    end

    % -------------------------------------------------------------
    % Read GeoTIFF metadata
    % -------------------------------------------------------------
    try
        info = geotiffinfo(fullfilename);
    catch
        info = io.geotiffinfoT8(fullfilename);
    end

    % -------------------------------------------------------------
    % Resolve CRS / EPSG robustly
    % -------------------------------------------------------------
    [isProjected, p, epsgCode] = resolveProjectedCRS(info, R);

    % -------------------------------------------------------------
    % Geographic raster branch
    % -------------------------------------------------------------
    if ~isProjected
        if isOutput(3) || isOutput(4)
            warning('Geographic coordinate system detected. X and Y will be output as longitude and latitude.');
            X = single(ones(R.RasterSize(1),1) * ...
                linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2)));
            Y = single(linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1))' * ...
                ones(1, R.RasterSize(2)));
        end

        if isOutput(5) || isOutput(6)
            lon = single(ones(R.RasterSize(1),1) * ...
                linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2)));
            lat = single(linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1))' * ...
                ones(1, R.RasterSize(2)));
        end

        if isOutput(7) || isOutput(8)
            warning('Geographic coordinate system detected. UTM outputs not computed.');
        end

        return
    end

    % -------------------------------------------------------------
    % Projected raster branch
    % -------------------------------------------------------------
    needXY   = isOutput(3) || isOutput(4);
    needLL   = isOutput(5) || isOutput(6);
    needUTM  = isOutput(7) || isOutput(8);
    needEPSG = isOutput(9);

    if needXY || needLL || needUTM || needEPSG
        nRows = R.RasterSize(1);
        nCols = R.RasterSize(2);

        xvec = single(linspace(R.XWorldLimits(1), R.XWorldLimits(2), nCols));
        yvec = single(linspace(R.YWorldLimits(1), R.YWorldLimits(2), nRows));

        Xfull = single(ones(nRows,1) * xvec);
        Yfull = single(yvec(:) * ones(1,nCols));
        % Equivalent to meshgrid-style layout with rows aligned to raster
    else
        Xfull = [];
        Yfull = [];
    end

    if isOutput(3), X = Xfull; end
    if isOutput(4), Y = Yfull; end

    if needLL || needUTM
        [latD, lonD] = projinv(p, double(Xfull), double(Yfull));

        if isOutput(5), lon = single(lonD); end
        if isOutput(6), lat = single(latD); end

        if needUTM
            if epsgCode > 0
                utmprojection = projcrs(epsgCode);
                [utmXD, utmYD] = projfwd(utmprojection, latD, lonD);

                % Preserve your original behavior
                utmYD = flipud(utmYD);

                if isOutput(7), utmX = single(utmXD); end
                if isOutput(8), utmY = single(utmYD); end
            end
        end
    end

    if ~needEPSG
        epsgCode = 0;
    end
end

% =====================================================================
function [isProjected, p, epsgCode] = resolveProjectedCRS(info, R)
% Resolve projected CRS robustly from geotiffinfo and/or raster reference.

    isProjected = false;
    p = [];
    epsgCode = 0;

    % 1) Best case: geotiffinfo SpatialRef.ProjectedCRS
    if isstruct(info) && isfield(info,'SpatialRef') && ~isempty(info.SpatialRef)
        try
            if isprop(info.SpatialRef, 'ProjectedCRS') && isa(info.SpatialRef.ProjectedCRS, 'projcrs')
                p = info.SpatialRef.ProjectedCRS;
                isProjected = true;
            end
        catch
        end
    end

    % 2) Try raster reference ProjectedCRS
    if isempty(p)
        try
            if isprop(R, 'ProjectedCRS') && isa(R.ProjectedCRS, 'projcrs')
                p = R.ProjectedCRS;
                isProjected = true;
            end
        catch
        end
    end

    % 3) EPSG from GeoTIFFCodes.PCS if present
    if isstruct(info) && isfield(info,'GeoTIFFCodes') && ~isempty(info.GeoTIFFCodes)
        if isfield(info.GeoTIFFCodes,'PCS') && ~isempty(info.GeoTIFFCodes.PCS)
            epsgCode = info.GeoTIFFCodes.PCS;
        elseif isfield(info.GeoTIFFCodes,'GCS') && ~isempty(info.GeoTIFFCodes.GCS)
            epsgCode = info.GeoTIFFCodes.GCS;
        end
    end

    % 4) If we have EPSG but not p, construct projcrs
    if isempty(p) && epsgCode > 0
        try
            p = projcrs(epsgCode);
            isProjected = true;
        catch
        end
    end

    % 5) Final fallback for your current DEM family
    % Metadata previously showed WGS84 / UTM zone 11N.
    if isempty(p) && isfield(info, 'GeoTIFFTags')
        try
            p = projcrs(32611);
            epsgCode = 32611;
            isProjected = true;
        catch
        end
    end
end
% function [A,R,X,Y,lon,lat,utmX,utmY,epsgCode] = readLidarTif(fullfilename)
% 
% % Check outputs
% isTilde  = utils.detectOutputSuppression(nargout);
% isOutput = ~isTilde;
% 
% if nargin < 2
%     debug = false;
% end
% 
% % Read raster as single
% try
%     [A,R] = readgeoraster(fullfilename, 'OutputType', 'single');
% catch
%     [A, R, info] = io.readLidarTif_big(fullfilename,true);
% end
% 
% try
%     info = geotiffinfo(fullfilename);
% catch
%     info = io.geotiffinfoT8(fullfilename);
% end
% 
% 
% 
% if isempty(info.GeoTIFFCodes.PCS)
%     % Geographic coordinate system
%     if isOutput(3) || isOutput(4)
%         warning('Global Coordinate System Detected. X and Y coordinates will be output as Longitude and Latitude')
%         X = single(ones(R.RasterSize(1),1) * linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2)));
%         Y = single(linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1))' * ones(1, R.RasterSize(2)));
%     else
%         X = 0; Y = 0;
%     end
% 
%     if isOutput(5) || isOutput(6)
%         lon = single(ones(R.RasterSize(1),1) * linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2)));
%         lat = single(linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), R.RasterSize(1))' * ones(1, R.RasterSize(2)));
%     else
%         lon = 0; lat = 0;
%     end
% 
%     if isOutput(9)
%         epsgCode = info.GeoTIFFCodes.GCS;
%     else
%         epsgCode = 0;
%     end
% 
%     if isOutput(7) || isOutput(8)
%         warning('Global Coordinate System Detected.')
%         utmX = 0; utmY = 0;
%     else
%         utmX = 0; utmY = 0;
%     end
% 
% else
%     % Projected coordinate system
%     X = single(ones(R.RasterSize(1),1) * linspace(R.XWorldLimits(1), R.XWorldLimits(2), R.RasterSize(2)));
%     Y = single(linspace(R.YWorldLimits(1), R.YWorldLimits(2), R.RasterSize(1))' * ones(1, R.RasterSize(2)));
% 
%     if (all(isOutput(5:6)) || isOutput(9))
%         % Use double internally for projinv/projfwd, then cast outputs to single
%         [latD, lonD] = projinv(R.ProjectedCRS, double(X), double(Y));
% 
%         if isOutput(5)
%             lon = single(lonD);
%         else
%             lon = 0;
%         end
% 
%         if isOutput(6)
%             lat = single(latD);
%         else
%             lat = 0;
%         end
% 
%         if all(isOutput(7:8)) || isOutput(9)
%             epsgCode = info.GeoTIFFCodes.PCS;
%             utmprojection = projcrs(epsgCode);
% 
%             [utmXD, utmYD] = projfwd(utmprojection, latD, lonD);
% 
%             % Preserve original behavior exactly
%             utmYD = flipud(utmYD);
% 
%             if isOutput(7)
%                 utmX = single(utmXD);
%             else
%                 utmX = 0;
%             end
% 
%             if isOutput(8)
%                 utmY = single(utmYD);
%             else
%                 utmY = 0;
%             end
%         else
%             utmX = 0; utmY = 0; epsgCode = 0;
%         end
%     else
%         lat = 0; lon = 0;
%         utmX = 0; utmY = 0;
%         epsgCode = 0;
%     end
% end
% 
% end
% % % function [A,R,X,Y,lon,lat,utmX,utmY,epsgCode] = readLidarTif(fullfilename)
% % % 
% % % % Check outputs
% % % isTilde  = utils.detectOutputSuppression(nargout);
% % % isOutput = ~isTilde;
% % % 
% % % % Defaults for optional outputs
% % % X = [];
% % % Y = [];
% % % lon = [];
% % % lat = [];
% % % utmX = [];
% % % utmY = [];
% % % epsgCode = 0;
% % % 
% % % % Read raster directly as single
% % % [A,R] = readgeoraster(fullfilename, 'OutputType', 'single');
% % % 
% % % try
% % %     info = geotiffinfo(fullfilename);
% % % catch
% % %     info = io.geotiffinfoT8(fullfilename);
% % % end
% % % 
% % % hasProjectedCRS = ~isempty(info.GeoTIFFCodes.PCS);
% % % 
% % % if ~hasProjectedCRS
% % %     % Geographic raster
% % %     if (nargout >= 3 && isOutput(3)) || (nargout >= 4 && isOutput(4))
% % %         warning('Global Coordinate System Detected. X and Y coordinates will be output as Longitude and Latitude')
% % % 
% % %         xvec = single(linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2)));
% % %         yvec = single(linspace(R.LatitudeLimits(1),  R.LatitudeLimits(2),  R.RasterSize(1)));
% % % 
% % %         X = repmat(xvec, R.RasterSize(1), 1);
% % %         Y = repmat(yvec(:), 1, R.RasterSize(2));
% % %     end
% % % 
% % %     if (nargout >= 5 && isOutput(5)) || (nargout >= 6 && isOutput(6))
% % %         if isempty(X) || isempty(Y)
% % %             lonvec = single(linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), R.RasterSize(2)));
% % %             latvec = single(linspace(R.LatitudeLimits(1),  R.LatitudeLimits(2),  R.RasterSize(1)));
% % % 
% % %             lon = repmat(lonvec, R.RasterSize(1), 1);
% % %             lat = repmat(latvec(:), 1, R.RasterSize(2));
% % %         else
% % %             lon = X;
% % %             lat = Y;
% % %         end
% % %     end
% % % 
% % %     if nargout >= 9 && isOutput(9)
% % %         epsgCode = info.GeoTIFFCodes.GCS;
% % %     end
% % % 
% % %     if (nargout >= 7 && isOutput(7)) || (nargout >= 8 && isOutput(8))
% % %         warning('Global Coordinate System Detected.')
% % %         utmX = single(0);
% % %         utmY = single(0);
% % %     end
% % % 
% % % else
% % %     % Projected raster: build projected X/Y only if needed
% % %     if (nargout >= 3 && isOutput(3)) || (nargout >= 4 && isOutput(4)) || ...
% % %        (nargout >= 5 && isOutput(5)) || (nargout >= 6 && isOutput(6)) || ...
% % %        (nargout >= 7 && isOutput(7)) || (nargout >= 8 && isOutput(8)) || ...
% % %        (nargout >= 9 && isOutput(9))
% % % 
% % %         xvec = single(linspace(R.XWorldLimits(1), R.XWorldLimits(2), R.RasterSize(2)));
% % %         yvec = single(linspace(R.YWorldLimits(1), R.YWorldLimits(2), R.RasterSize(1)));
% % % 
% % %         Xfull = repmat(xvec, R.RasterSize(1), 1);
% % %         Yfull = repmat(yvec(:), 1, R.RasterSize(2));
% % %     else
% % %         Xfull = [];
% % %         Yfull = [];
% % %     end
% % % 
% % %     if (nargout >= 3 && isOutput(3))
% % %         X = Xfull;
% % %     end
% % %     if (nargout >= 4 && isOutput(4))
% % %         Y = Yfull;
% % %     end
% % % 
% % %     % Compute lat/lon only if requested
% % %     if ((nargout >= 5 && isOutput(5)) || (nargout >= 6 && isOutput(6)) || ...
% % %         (nargout >= 7 && isOutput(7)) || (nargout >= 8 && isOutput(8)) || ...
% % %         (nargout >= 9 && isOutput(9)))
% % % 
% % %         [latD, lonD] = projinv(R.ProjectedCRS, double(Xfull), double(Yfull));
% % %         latD = flipud(latD);
% % %         lonD = flipud(lonD);
% % % 
% % %         if nargout >= 5 && isOutput(5)
% % %             lon = single(lonD);
% % %         end
% % %         if nargout >= 6 && isOutput(6)
% % %             lat = single(latD);
% % %         end
% % % 
% % %         % EPSG
% % %         if nargout >= 9 && isOutput(9)
% % %             epsgCode = info.GeoTIFFCodes.PCS;
% % %         else
% % %             epsgCode = info.GeoTIFFCodes.PCS;
% % %         end
% % % 
% % %         % Compute UTM only if requested
% % %         if (nargout >= 7 && isOutput(7)) || (nargout >= 8 && isOutput(8))
% % %             utmprojection = projcrs(epsgCode);
% % %             [utmXD, utmYD] = projfwd(utmprojection, latD, lonD);
% % %             utmYD = flipud(utmYD);
% % % 
% % %             if nargout >= 7 && isOutput(7)
% % %                 utmX = single(utmXD);
% % %             end
% % %             if nargout >= 8 && isOutput(8)
% % %                 utmY = single(utmYD);
% % %             end
% % %         end
% % % 
% % %     else
% % %         if nargout >= 9 && isOutput(9)
% % %             epsgCode = info.GeoTIFFCodes.PCS;
% % %         end
% % %     end
% % % end
% % % 
% % % end
% % function [A,R,X,Y,lon,lat,utmX,utmY,epsgCode] = readLidarTif(fullfilename)
% % 
% % % Check outPuts
% % % https://www.mathworks.com/matlabcentral/fileexchange/79218-detectoutputsuppression
% % isTilde = utils.detectOutputSuppression(nargout);isOutput = ~isTilde;
% % 
% % % Extract the Geotiff data
% % % https://www.mathworks.com/matlabcentral/answers/517443-how-i-can-view-this-tif-file
% % [A,R] = readgeoraster(fullfilename, 'OutputType', 'single');
% % try
% %     info = geotiffinfo(fullfilename);
% % catch
% %     % https://www.mathworks.com/help/map/ref/geotiffinfo.html
% %     info = geotiffinfoT8(fullfilename);
% % end
% % if isempty(info.GeoTIFFCodes.PCS)
% %     if isOutput(3) | isOutput(4)
% %         warning('Global Coordinate System Detected. X and Y coordinates will be output as Longitude and Latitude')
% %         X = ones(R.RasterSize(1),1)*linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
% %         Y = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1))'*ones(1,R.RasterSize(2));
% %     else
% %         X = 0; Y = 0;
% %     end
% %     if isOutput(5)|isOutput(6)
% %         lon = ones(R.RasterSize(1),1)*linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
% %         lat = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1))'*ones(1,R.RasterSize(2));
% %     else
% %         lon = 0; lat = 0;
% %     end
% %     if isOutput(9)
% %         epsgCode = info.GeoTIFFCodes.GCS;
% %     else
% %         epsgCode = 0;
% %     end
% %     if isOutput(7)|isOutput(8)
% %         warning('Global Coordinate System Detected.')
% %         utmX = 0; utmY = 0;
% %     else
% %         utmX = 0; utmY = 0;
% %     end
% % else
% %     % Get X,Y MeshGrid like Matrix
% %     X = ones(R.RasterSize(1),1)*linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
% %     Y = linspace(R.YWorldLimits(1),R.YWorldLimits(2),R.RasterSize(1))'*ones(1,R.RasterSize(2));
% % 
% %     % Extract Lat Lon
% %     if all(isOutput(5:6))|isOutput(9)
% %         [lat,lon] = projinv(R.ProjectedCRS,X,Y);
% %         % Convert to UTM
% %         if all(isOutput(7:8))|isOutput(9)
% %             % EPSG Code from info
% %             epsgCode = info.GeoTIFFCodes.PCS;
% %             utmprojection = projcrs(epsgCode);
% %             % EPSG Code for Grand Mesa, CO Zone 12
% %             % https://georepository.com/crs_32612/WGS-84-UTM-zone-12N.html
% %             %         utmprojection = projcrs(32612);
% %             % EPS Code for Grand Mesa, CO Zone 13
% %             % https://epsg.org/crs_32613/WGS-84-UTM-zone-13N.html?sessionkey=of4suotgv0
% %             %         utmprojection = projcrs(32613);
% %             [utmX,utmY] = projfwd(utmprojection,lat,lon);
% %             utmY = flipud(utmY);
% %         else
% %             utmX = 0; utmY = 0; epsgCode = 0;
% %         end
% %     else
% %         lat = 0; lon = 0;
% %         utmX = 0; utmY = 0;
% %         epsgCode = 0;
% %     end
% % end
% % 
% % end
% % 

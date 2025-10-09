function [A,R,X,Y,lon,lat,utmX,utmY,epsgCode] = readLidarTif(fullfilename)

% Check outPuts
% https://www.mathworks.com/matlabcentral/fileexchange/79218-detectoutputsuppression
isTilde = utils.detectOutputSuppression(nargout);isOutput = ~isTilde;

% Extract the Geotiff data
% https://www.mathworks.com/matlabcentral/answers/517443-how-i-can-view-this-tif-file
[A,R] = readgeoraster(fullfilename, 'OutputType', 'double');
try
    info = geotiffinfo(fullfilename);
catch
    % https://www.mathworks.com/help/map/ref/geotiffinfo.html
    info = geotiffinfoT8(fullfilename);
end
if isempty(info.GeoTIFFCodes.PCS)
    if isOutput(3) | isOutput(4)
        warning('Global Coordinate System Detected. X and Y coordinates will be output as Longitude and Latitude')
        X = ones(R.RasterSize(1),1)*linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
        Y = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1))'*ones(1,R.RasterSize(2));
    else
        X = 0; Y = 0;
    end
    if isOutput(5)|isOutput(6)
        lon = ones(R.RasterSize(1),1)*linspace(R.LongitudeLimits(1),R.LongitudeLimits(2),R.RasterSize(2));
        lat = linspace(R.LatitudeLimits(1),R.LatitudeLimits(2),R.RasterSize(1))'*ones(1,R.RasterSize(2));
    else
        lon = 0; lat = 0;
    end
    if isOutput(9)
        epsgCode = info.GeoTIFFCodes.GCS;
    else
        epsgCode = 0;
    end
    if isOutput(7)|isOutput(8)
        warning('Global Coordinate System Detected.')
        utmX = 0; utmY = 0;
    else
        utmX = 0; utmY = 0;
    end
else
    % Get X,Y MeshGrid like Matrix
    X = ones(R.RasterSize(1),1)*linspace(R.XWorldLimits(1),R.XWorldLimits(2),R.RasterSize(2));
    Y = linspace(R.YWorldLimits(1),R.YWorldLimits(2),R.RasterSize(1))'*ones(1,R.RasterSize(2));

    % Extract Lat Lon
    if all(isOutput(5:6))|isOutput(9)
        [lat,lon] = projinv(R.ProjectedCRS,X,Y);
        % Convert to UTM
        if all(isOutput(7:8))|isOutput(9)
            % EPSG Code from info
            epsgCode = info.GeoTIFFCodes.PCS;
            utmprojection = projcrs(epsgCode);
            % EPSG Code for Grand Mesa, CO Zone 12
            % https://georepository.com/crs_32612/WGS-84-UTM-zone-12N.html
            %         utmprojection = projcrs(32612);
            % EPS Code for Grand Mesa, CO Zone 13
            % https://epsg.org/crs_32613/WGS-84-UTM-zone-13N.html?sessionkey=of4suotgv0
            %         utmprojection = projcrs(32613);
            [utmX,utmY] = projfwd(utmprojection,lat,lon);
            utmY = flipud(utmY);
        else
            utmX = 0; utmY = 0; epsgCode = 0;
        end
    else
        lat = 0; lon = 0;
        utmX = 0; utmY = 0;
        epsgCode = 0;
    end
end

end


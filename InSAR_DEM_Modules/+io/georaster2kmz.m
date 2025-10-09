function georaster2kmz(A, R, outFile, varargin)
%GEORASTER2KMZ  Export a raster + spatial reference to KML or KMZ GroundOverlay,
%               with optional legend (ScreenOverlay) and clickable sampled values.
%
% Example (KMZ):
%   georaster2kmz(A,R,'overlay.kmz','CRS',32612,'Name','Incidence (deg)', ...
%       'CLim',[0 60],'AddLegend',true,'LegendLabel','Incidence (deg)', ...
%       'SampleStep',250,'Marker','dot','MarkerSize',8,'MarkerColor',[0 0 0]);
%
% Example (KML + sibling PNGs):
%   georaster2kmz(A,R,'overlay.kml', ... same options ...);
% Name-Value options
%   'CRS'           : REQUIRED if R is Map*Reference (projected). EPSG (numeric),
%                     a projcrs, or a GeoTIFF filename whose CRS matches R.
%   'Name'          : overlay title (default: basename of outFile)
%   'ImageFormat'   : 'png' (default) or 'jpg'
%   'CLim'          : [min max] for single-band stretch (default: auto 2–98%)
%   'Colormap'      : Nx3 colormap for single-band (default: parula(256))
%   'TransparentIf' : value to make transparent (NaNs always transparent; default NaN only)
%   'AddLegend'     : true/false (default false) – adds a ScreenOverlay legend
%   'LegendLabel'   : label string (default '')
%   'LegendSize'    : [width height] px (default [160 360])
%   'LegendCorner'  : 'NE'|'NW'|'SE'|'SW' (default 'NE')
%   'SampleStep'    : pixel spacing for clickable points (default Inf = off)
%   'MaxSamples'    : cap for number of placemarks (default 3000)
%   'ValueFormat'   : sprintf format for value text (default '%.3g')
%   'Units'         : units suffix for values (default '')
%   'Marker'        : Incidence Angle Marker (default 'dot')
%   'MarkerSize'    : default 8
%   'MarkerColor'   : default [0,0,0] black
%   'MarkserScale'  : default 1.0
%   'ShowLabels'    : default true

% ---------- Parse options ----------
defaults.Name          = '';
defaults.ImageFormat   = 'png';
defaults.CLim          = [NaN NaN];
defaults.Colormap      = parula(256);
defaults.TransparentIf = NaN;
defaults.CRS           = [];
defaults.AddLegend     = false;
defaults.LegendLabel   = '';
defaults.LegendSize    = [160 360];
defaults.LegendCorner  = 'NE';
defaults.SampleStep    = Inf;
defaults.MaxSamples    = 3000;
defaults.ValueFormat   = '%.3g';
defaults.Units         = '';
defaults.Marker        = 'dot';     % 'dot' | 'square'
defaults.MarkerSize    = 30;         % pixels
defaults.MarkerColor   = [0 0 0];   % black (RGB)
defaults.MarkerScale   = .5;       % KML scale factor
defaults.ShowLabels    = true;     % show point labels by default

opts = parseNV(varargin, defaults);
fmt  = lower(char(opts.ImageFormat));
if ~ismember(fmt, {'png','jpg'})
    error('ImageFormat must be ''png'' or ''jpg''.');
end

% Output mode by extension
[outDir, outBase, outExt] = fileparts(outFile);
isKMZ = strcmpi(outExt, '.kmz');
isKML = strcmpi(outExt, '.kml');
if ~(isKMZ || isKML)
    error('Output must end with .kmz or .kml');
end

% ---------- Build RGB + alpha (and get effective CLim for legend) ----------
[rgb, alpha, climEff] = toUint8RGB(A, opts.CLim, opts.Colormap, opts.TransparentIf);

% ---------- Compute LatLonBox (north/south/east/west) ----------
if isprop(R,'LatitudeLimits') % Geographic*Reference (lat/lon)
    latlim = R.LatitudeLimits;  lonlim = R.LongitudeLimits;
    if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
        dlat = diff(latlim) / max(R.RasterSize(1)-1, 1);
        dlon = diff(lonlim) / max(R.RasterSize(2)-1, 1);
        latlim = latlim + [-0.5 0.5]*dlat;
        lonlim = lonlim + [-0.5 0.5]*dlon;
    end
    north = latlim(2); south = latlim(1); east = lonlim(2); west = lonlim(1);
else % Map*Reference (projected) -> need CRS to projinv corners
    xw = R.XWorldLimits;  yw = R.YWorldLimits;
    if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
        dx = diff(xw) / max(R.RasterSize(2)-1, 1);
        dy = diff(yw) / max(R.RasterSize(1)-1, 1);
        xw = xw + [-0.5 0.5]*dx;
        yw = yw + [-0.5 0.5]*dy;
    end
    x = [xw(1) xw(2) xw(2) xw(1)];
    y = [yw(2) yw(2) yw(1) yw(1)];
    crs = normalizeCRS(opts.CRS);
    if isempty(crs)
        error(['Projected reference detected. Provide ''CRS'' as EPSG code, projcrs, ' ...
               'or a GeoTIFF filename with matching CRS.']);
    end
    [lat,lon] = projinv(crs, x, y);
    north = max(lat); south = min(lat); east = max(lon); west = min(lon);
end

% ---------- File staging ----------
tmp = tempname; mkdir(tmp);
if isKMZ
    imgDir  = fullfile(tmp,'files'); mkdir(imgDir);
else
    imgDir  = tmp;
end
imgName = [outBase '.' fmt];
if isKMZ
    imgRel = fullfile('files', imgName);
else
    imgRel = imgName;
end
imgPath = fullfile(imgDir, imgName);

% Write overlay image
if strcmp(fmt,'png')
    imwrite(rgb, imgPath, 'Alpha', alpha);
else
    imwrite(rgb, imgPath, 'Quality', 95);
end

% Optional legend (ScreenOverlay)
legendRel = ''; legendName = ''; legDir = '';
if opts.AddLegend
    legendName = [outBase '_legend.png'];
    if isKMZ
        legDir    = imgDir;                      % put legend in 'files' too
        legendRel = fullfile('files', legendName);
    else
        legDir    = tmp;
        legendRel = legendName;
    end
    legendPath = fullfile(legDir, legendName);
    makeLegendPNG(legendPath, opts.Colormap, climEff, opts.LegendLabel, opts.LegendSize);
end

% Optional marker icon for samples
markerRel = ''; markerName = ''; markDir = ''; kmlColor = '#ff000000';
if isfinite(opts.SampleStep) && opts.SampleStep>=1
    markerName = sprintf('%s_marker_%s_%dpx.png', outBase, lower(opts.Marker), round(opts.MarkerSize));
    if isKMZ
        markDir   = imgDir;
        markerRel = fullfile('files', markerName);
    else
        markDir   = tmp;
        markerRel = markerName;
    end
    markerPath = fullfile(markDir, markerName);
    makeMarkerPNG(markerPath, lower(opts.Marker), round(opts.MarkerSize));
    kmlColor = rgbToKmlColor(opts.MarkerColor, 255); % ABGR hex
end

% ---------- Write KML (doc.kml in temp; KMZ expects that as root) ----------
kmlTempPath = fullfile(tmp, 'doc.kml');
overlayName = opts.Name; if isempty(overlayName), overlayName = outBase; end
overlayName = char(overlayName);
imgHref = strrep(imgRel, filesep, '/');           % forward slashes
if ~isempty(legendRel), legendHref = strrep(legendRel, filesep, '/'); end
if ~isempty(markerRel), markerHref = strrep(markerRel, filesep, '/'); end

fid = fopen(kmlTempPath,'w'); assert(fid>0, 'Failed to create KML.');
fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n');
fprintf(fid, '  <name>%s</name>\n', overlayName);

% --- Raster GroundOverlay ---
fprintf(fid, '  <GroundOverlay>\n');
fprintf(fid, '    <name>%s</name>\n', overlayName);
fprintf(fid, '    <Icon><href>%s</href></Icon>\n', imgHref);
fprintf(fid, '    <LatLonBox>\n');
fprintf(fid, '      <north>%.10f</north>\n', north);
fprintf(fid, '      <south>%.10f</south>\n', south);
fprintf(fid, '      <east>%.10f</east>\n',  east );
fprintf(fid, '      <west>%.10f</west>\n',  west );
fprintf(fid, '      <rotation>0</rotation>\n');
fprintf(fid, '    </LatLonBox>\n');
fprintf(fid, '  </GroundOverlay>\n');

% --- Optional ScreenOverlay legend ---
if opts.AddLegend
    [ox, oy] = cornerXY(opts.LegendCorner);
    fprintf(fid, '  <ScreenOverlay>\n');
    fprintf(fid, '    <name>Legend</name>\n');
    fprintf(fid, '    <Icon><href>%s</href></Icon>\n', legendHref);
    fprintf(fid, '    <overlayXY x="%.3f" y="%.3f" xunits="fraction" yunits="fraction"/>\n', ox, oy);
    fprintf(fid, '    <screenXY  x="%.3f" y="%.3f" xunits="fraction" yunits="fraction"/>\n', ox, oy);
    fprintf(fid, '    <size x="%d" y="%d" xunits="pixels" yunits="pixels"/>\n', opts.LegendSize(1), opts.LegendSize(2));
    fprintf(fid, '  </ScreenOverlay>\n');
end

% --- Optional sampled value placemarks (using simple local icon) ---
if isfinite(opts.SampleStep) && opts.SampleStep>=1
    writeSamplePlacemarks(fid, A, R, opts, crsIfMapRef(R, opts.CRS), ...
        markerHref, kmlColor);
end

% Close KML
fprintf(fid, '</Document>\n</kml>\n');
fclose(fid);

% ---------- Package to destination ----------
cwd = pwd;
try
    if isKMZ
        % zip then rename to .kmz
        cd(tmp);
        zipName = fullfile(tmp,'package.zip');
        zip(zipName, {'doc.kml','files'});
        cd(cwd);
        if ~isempty(outDir) && exist(outDir,'dir')==0, mkdir(outDir); end
        movefile(zipName, outFile, 'f');
    else
        if ~isempty(outDir) && exist(outDir,'dir')==0, mkdir(outDir); end
        destKML = fullfile(outDir, [outBase '.kml']);
        destIMG = fullfile(outDir, imgName);
        copyfile(imgPath, destIMG, 'f');
        movefile(kmlTempPath, destKML, 'f');
        if opts.AddLegend
            copyfile(fullfile(legDir, legendName), fullfile(outDir, legendName), 'f');
        end
        if isfinite(opts.SampleStep) && opts.SampleStep>=1
            copyfile(fullfile(markDir, markerName), fullfile(outDir, markerName), 'f');
        end
    end
catch ME
    cd(cwd);
    rethrow(ME);
end
end

% ================== helpers ==================
function opts = parseNV(nv, defaults)
    opts = defaults;
    if isempty(nv), return; end
    if mod(numel(nv),2)~=0, error('Options must be name-value pairs.'); end
    for k = 1:2:numel(nv)
        name = nv{k};  val = nv{k+1};
        if ~ischar(name) && ~isstring(name), error('Option names must be strings.'); end
        name = char(name);
        if ~isfield(opts, name), error('Unknown option: %s', name); end
        opts.(name) = val;
    end
end

function [rgb, alpha, climUsed] = toUint8RGB(A, CLim, cmap, TransparentIf)
    if ndims(A)==2
        data = double(A);
        finiteMask = isfinite(data);

        if any(isnan(CLim))
            if any(finiteMask(:))
                p = prctile(data(finiteMask), [2 98]);
            else
                p = [0 1];
            end
        else
            p = CLim;
        end
        if ~all(isfinite(p)) || numel(p)~=2
            p = [0 1];
        end
        if p(1)==p(2)
            d = max(1,abs(p(1)))*1e-6;
            p = p + [-d d];
        end
        climUsed = p;

        Inorm = mat2gray(data, p);
        rgb   = ind2rgb(gray2ind(Inorm, size(cmap,1)), cmap);
        if isnan(TransparentIf)
            tmask = isnan(data);
        else
            tmask = isnan(data) | (data==TransparentIf);
        end
        alpha = uint8(255 * (~tmask));
    elseif size(A,3)==3
        if isa(A,'uint8')
            rgb = A;
        else
            mn = min(A(:)); mx = max(A(:));
            if isfloat(A) && isfinite(mn) && mn>=0 && isfinite(mx) && mx<=1
                rgb = im2uint8(A);
            else
                rgb = im2uint8(mat2gray(A));
            end
        end
        alpha = uint8(255*ones(size(A,1), size(A,2)));
        climUsed = [0 1];
    else
        error('A must be single-band or 3-band RGB.');
    end
end

function [ox, oy] = cornerXY(which)
    switch upper(which)
        case 'NE', ox = 1; oy = 1;
        case 'NW', ox = 0; oy = 1;
        case 'SE', ox = 1; oy = 0;
        case 'SW', ox = 0; oy = 0;
        otherwise, ox = 1; oy = 1;
    end
end

function crs = crsIfMapRef(R, CRSopt)
    if isprop(R,'LatitudeLimits')
        crs = [];
    else
        crs = normalizeCRS(CRSopt);
    end
end

function crs = normalizeCRS(CRSopt)
    crs = [];
    if isempty(CRSopt), return; end
    if isnumeric(CRSopt), crs = projcrs(CRSopt); return; end
    if isa(CRSopt,'projcrs'), crs = CRSopt; return; end
    if ischar(CRSopt) || isstring(CRSopt)
        f = char(CRSopt);
        if exist(f,'file')==2
            try
                info = georasterinfo(f);
                if isfield(info,'CoordinateReferenceSystem') && ~isempty(info.CoordinateReferenceSystem)
                    crs = info.CoordinateReferenceSystem; return;
                end
            catch, end
            try
                info = geotiffinfo(f);
                if isfield(info,'PCS') && ~isempty(info.PCS)
                    crs = projcrs(info.PCS); return;
                end
            catch, end
        end
    end
end

function makeLegendPNG(path, cmap, clim, label, sz)
    if numel(clim)~=2 || ~all(isfinite(clim)), clim = [0 1]; end
    if clim(1)==clim(2)
        d = max(1,abs(clim(1)))*1e-6;
        clim = clim + [-d d];
    end
    w = sz(1); h = sz(2);
    f = figure('Visible','off','Position',[100 100 w h], 'Color','w');
    axes('Position',[0.35 0.08 0.4 0.84]); %#ok<LAXES>
    g = linspace(clim(1), clim(2), 256)'; % data-space gradient
    imagesc([0 1],[clim(1) clim(2)], g);
    set(gca,'YDir','normal','XTick',[],'Box','on');
    colormap(cmap); caxis(clim);
    ylabel(label, 'Interpreter','none');
    try
        exportgraphics(f, path, 'BackgroundColor','white', 'Resolution', 150);
    catch
        set(f,'InvertHardcopy','off');
        print(f, path, '-dpng','-r150');
    end
    close(f);
end

function makeMarkerPNG(path, shape, sz)
    % Create a tiny PNG icon with transparent background and an opaque shape.
    % PNG must be MxN or MxNx3; pass alpha mask separately.

    sz = max(2, round(sz));

    % Base RGB (white so KML <color> tint works nicely)
    rgb   = uint8(255 * ones(sz, sz, 3));

    % Alpha mask: 0 = transparent, 255 = opaque
    alpha = zeros(sz, sz, 'uint8');

    switch lower(shape)
        case 'dot'
            [X,Y] = meshgrid(1:sz, 1:sz);
            cx = (sz+1)/2; cy = (sz+1)/2;
            r  = sz/2 * 0.9;
            M  = (X - cx).^2 + (Y - cy).^2 <= r^2;

        case 'square'
            pad = max(1, round(0.1*sz));
            M = false(sz);
            M(1+pad:sz-pad, 1+pad:sz-pad) = true;

        otherwise
            % fallback to small square
            pad = max(1, round(0.1*sz));
            M = false(sz);
            M(1+pad:sz-pad, 1+pad:sz-pad) = true;
    end

    alpha(M) = 255;

    % Write as RGB + separate alpha channel
    imwrite(rgb, path, 'png', 'Alpha', alpha);
end


function hex = rgbToKmlColor(rgb, alpha)
    % Convert [r g b] (0-1 or 0-255) + alpha (0-255) to KML AABBGGRR hex
    rgb = double(rgb(:)');
    if max(rgb)<=1, rgb = round(255*rgb); end
    r = rgb(1); g = rgb(2); b = rgb(3);
    a = max(0,min(255,round(alpha)));
    hex = sprintf('#%02x%02x%02x%02x', a, b, g, r); % AABBGGRR
end

function writeSamplePlacemarks(fid, A, R, opts, crs, markerHref, kmlColor)
    % Sampling grid from center; fallback to one valid point if needed.
    step = max(1, round(opts.SampleStep));
    [nrows, ncols, ~] = size(A);
    r0 = max(1, min(nrows, round(nrows/2)));
    c0 = max(1, min(ncols, round(ncols/2)));
    rr = unique([r0:-step:1, r0:step:nrows]);
    cc = unique([c0:-step:1, c0:step:ncols]);
    [RR, CC] = ndgrid(rr, cc);

    if ndims(A)==2
        vals = A(sub2ind([nrows ncols], RR(:), CC(:)));
    else
        vals = A(:,:,1);
        vals = vals(sub2ind([nrows ncols], RR(:), CC(:)));
    end
    good = isfinite(vals);
    RR = RR(good); CC = CC(good); vals = vals(good);

    if isempty(vals)
        if ndims(A)==2
            [ii,jj] = find(isfinite(A), 1, 'first');
        else
            [ii,jj] = find(isfinite(A(:,:,1)), 1, 'first');
        end
        if isempty(ii), return; end
        RR = ii; CC = jj; vals = A(ii,jj);
    end

    if numel(vals) > opts.MaxSamples
        idx = randperm(numel(vals), opts.MaxSamples);
        RR = RR(idx); CC = CC(idx); vals = vals(idx);
    end
    npts = numel(vals);

    if isprop(R,'LatitudeLimits')
        [lat, lon] = intrinsicToGeographic(R, CC, RR);
    else
        [x, y] = intrinsicToWorld(R, CC, RR);
        [lat, lon] = projinv(crs, x, y);
    end

    % Simple local style using tiny icon and color tint; label scale optional
    fprintf(fid, '  <Style id="ptStyle">\n');
    fprintf(fid, '    <IconStyle>\n');
    fprintf(fid, '      <scale>%.2f</scale>\n', opts.MarkerScale);
    fprintf(fid, '      <color>%s</color>\n', kmlColor);
    fprintf(fid, '      <Icon><href>%s</href></Icon>\n', markerHref);
    fprintf(fid, '    </IconStyle>\n');
    % Show/hide labels: 0 = hidden, 0.8 = readable
    labelScale = 0.0;
    if islogical(opts.ShowLabels) && opts.ShowLabels
        labelScale = 0.8;
    elseif isnumeric(opts.ShowLabels) && opts.ShowLabels~=0
        labelScale = 0.8;
    end
    fprintf(fid, '    <LabelStyle><scale>%.2f</scale></LabelStyle>\n', labelScale);
    fprintf(fid, '  </Style>\n');

    fprintf(fid, '  <Folder>\n');
    fprintf(fid, '    <name>Samples (%d)</name>\n', npts);
    fprintf(fid, '    <visibility>1</visibility>\n');

    fmtVal = opts.ValueFormat;
    hasUnits = ~isempty(opts.Units);

    for k = 1:npts
        vtxt = sprintf(fmtVal, vals(k));
        if hasUnits, vtxt = [vtxt ' ' opts.Units]; end
        fprintf(fid, '    <Placemark>\n');
        fprintf(fid, '      <styleUrl>#ptStyle</styleUrl>\n');
        fprintf(fid, '      <name>%s</name>\n', vtxt);
        fprintf(fid, '      <description><![CDATA[Value: <b>%s</b>]]></description>\n', vtxt);
        fprintf(fid, '      <Point><coordinates>%.10f,%.10f,0</coordinates></Point>\n', lon(k), lat(k));
        fprintf(fid, '      <ExtendedData><Data name="value"><value>%s</value></Data></ExtendedData>\n', vtxt);
        fprintf(fid, '    </Placemark>\n');
    end
    fprintf(fid, '  </Folder>\n');
end

% function georaster2kmz(A, R, outFile, varargin)
% %GEORASTER2KMZ  Export a raster + spatial reference to KML or KMZ GroundOverlay,
% %               with optional legend (ScreenOverlay) and clickable sampled values.
% %
% % KMZ:
% %   georaster2kmz(A,R,'overlay.kmz','CRS',32612,'Name','Incidence (deg)', ...
% %       'CLim',[0 60],'AddLegend',true,'LegendLabel','Incidence (deg)', ...
% %       'SampleStep',200,'ValueFormat','%.2f','Units','deg');
% %
% % KML (image placed next to KML):
% %   georaster2kmz(A,R,'overlay.kml', ... same options ...);
% %
% % Name-Value options
% %   'CRS'           : REQUIRED if R is Map*Reference (projected). EPSG (numeric),
% %                     a projcrs, or a GeoTIFF filename whose CRS matches R.
% %   'Name'          : overlay title (default: basename of outFile)
% %   'ImageFormat'   : 'png' (default) or 'jpg'
% %   'CLim'          : [min max] for single-band stretch (default: auto 2–98%)
% %   'Colormap'      : Nx3 colormap for single-band (default: parula(256))
% %   'TransparentIf' : value to make transparent (NaNs always transparent; default NaN only)
% %   'AddLegend'     : true/false (default false) – adds a ScreenOverlay legend
% %   'LegendLabel'   : label string (default '')
% %   'LegendSize'    : [width height] px (default [160 360])
% %   'LegendCorner'  : 'NE'|'NW'|'SE'|'SW' (default 'NE')
% %   'SampleStep'    : pixel spacing for clickable points (default Inf = off)
% %   'MaxSamples'    : cap for number of placemarks (default 3000)
% %   'ValueFormat'   : sprintf format for value text (default '%.3g')
% %   'Units'         : units suffix for values (default '')
% 
% % ---------- Parse options ----------
% defaults.Name          = '';
% defaults.ImageFormat   = 'png';
% defaults.CLim          = [NaN NaN];
% defaults.Colormap      = parula(256);
% defaults.TransparentIf = NaN;
% defaults.CRS           = [];
% defaults.AddLegend     = false;
% defaults.LegendLabel   = '';
% defaults.LegendSize    = [160 360];
% defaults.LegendCorner  = 'NE';
% defaults.SampleStep    = Inf;
% defaults.MaxSamples    = 3000;
% defaults.ValueFormat   = '%.3g';
% defaults.Units         = '';
% 
% opts = parseNV(varargin, defaults);
% fmt  = lower(char(opts.ImageFormat));
% if ~ismember(fmt, {'png','jpg'})
%     error('ImageFormat must be ''png'' or ''jpg''.');
% end
% 
% % Output mode by extension
% [outDir, outBase, outExt] = fileparts(outFile);
% isKMZ = strcmpi(outExt, '.kmz');
% isKML = strcmpi(outExt, '.kml');
% if ~(isKMZ || isKML)
%     error('Output must end with .kmz or .kml');
% end
% 
% % ---------- Build RGB + alpha (and get effective CLim for legend) ----------
% [rgb, alpha, climEff] = toUint8RGB(A, opts.CLim, opts.Colormap, opts.TransparentIf);
% 
% % ---------- Compute LatLonBox (north/south/east/west) ----------
% if isprop(R,'LatitudeLimits') % Geographic*Reference (lat/lon)
%     latlim = R.LatitudeLimits;  lonlim = R.LongitudeLimits;
%     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
%         dlat = diff(latlim) / max(R.RasterSize(1)-1, 1);
%         dlon = diff(lonlim) / max(R.RasterSize(2)-1, 1);
%         latlim = latlim + [-0.5 0.5]*dlat;
%         lonlim = lonlim + [-0.5 0.5]*dlon;
%     end
%     north = latlim(2); south = latlim(1); east = lonlim(2); west = lonlim(1);
% else % Map*Reference (projected) -> need CRS to projinv corners
%     xw = R.XWorldLimits;  yw = R.YWorldLimits;
%     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
%         dx = diff(xw) / max(R.RasterSize(2)-1, 1);
%         dy = diff(yw) / max(R.RasterSize(1)-1, 1);
%         xw = xw + [-0.5 0.5]*dx;
%         yw = yw + [-0.5 0.5]*dy;
%     end
%     x = [xw(1) xw(2) xw(2) xw(1)];
%     y = [yw(2) yw(2) yw(1) yw(1)];
%     crs = normalizeCRS(opts.CRS);
%     if isempty(crs)
%         error(['Projected reference detected. Provide ''CRS'' as EPSG code, projcrs, ' ...
%                'or a GeoTIFF filename with matching CRS.']);
%     end
%     [lat,lon] = projinv(crs, x, y);
%     north = max(lat); south = min(lat); east = max(lon); west = min(lon);
% end
% 
% % ---------- File staging ----------
% tmp = tempname; mkdir(tmp);
% if isKMZ
%     imgDir  = fullfile(tmp,'files'); mkdir(imgDir);
% else
%     imgDir  = tmp;
% end
% imgName = [outBase '.' fmt];
% if isKMZ
%     imgRel = fullfile('files', imgName);   % path inside KMZ
% else
%     imgRel = imgName;                      % sibling to .kml
% end
% 
% imgPath = fullfile(imgDir, imgName);
% 
% % Write overlay image
% if strcmp(fmt,'png')
%     imwrite(rgb, imgPath, 'Alpha', alpha);
% else
%     imwrite(rgb, imgPath, 'Quality', 95);
% end
% 
% % Optional legend (ScreenOverlay) next to overlay image
% legendRel = ''; legendName = ''; legDir = '';
% if opts.AddLegend
%     legendName = [outBase '_legend.png'];
%     if isKMZ
%         legDir    = imgDir;                      % put legend in 'files' too
%         legendRel = fullfile('files', legendName);
%     else
%         legDir    = tmp;
%         legendRel = legendName;
%     end
%     legendPath = fullfile(legDir, legendName);
%     makeLegendPNG(legendPath, opts.Colormap, climEff, opts.LegendLabel, opts.LegendSize);
% end
% 
% % ---------- Write KML (doc.kml in temp; KMZ expects that as root) ----------
% kmlTempPath = fullfile(tmp, 'doc.kml');
% overlayName = opts.Name; if isempty(overlayName), overlayName = outBase; end
% overlayName = char(overlayName);
% imgHref = strrep(imgRel, filesep, '/');           % forward slashes for KMZ
% if ~isempty(legendRel), legendHref = strrep(legendRel, filesep, '/'); end
% 
% fid = fopen(kmlTempPath,'w'); assert(fid>0, 'Failed to create KML.');
% fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
% fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n');
% fprintf(fid, '  <name>%s</name>\n', overlayName);
% 
% % --- Raster GroundOverlay ---
% fprintf(fid, '  <GroundOverlay>\n');
% fprintf(fid, '    <name>%s</name>\n', overlayName);
% fprintf(fid, '    <Icon><href>%s</href></Icon>\n', imgHref);
% fprintf(fid, '    <LatLonBox>\n');
% fprintf(fid, '      <north>%.10f</north>\n', north);
% fprintf(fid, '      <south>%.10f</south>\n', south);
% fprintf(fid, '      <east>%.10f</east>\n',  east );
% fprintf(fid, '      <west>%.10f</west>\n',  west );
% fprintf(fid, '      <rotation>0</rotation>\n');
% fprintf(fid, '    </LatLonBox>\n');
% fprintf(fid, '  </GroundOverlay>\n');
% 
% % --- Optional ScreenOverlay legend ---
% if opts.AddLegend
%     [ox, oy] = cornerXY(opts.LegendCorner);
%     fprintf(fid, '  <ScreenOverlay>\n');
%     fprintf(fid, '    <name>Legend</name>\n');
%     fprintf(fid, '    <Icon><href>%s</href></Icon>\n', legendHref);
%     fprintf(fid, '    <overlayXY x="%.3f" y="%.3f" xunits="fraction" yunits="fraction"/>\n', ox, oy);
%     fprintf(fid, '    <screenXY  x="%.3f" y="%.3f" xunits="fraction" yunits="fraction"/>\n', ox, oy);
%     fprintf(fid, '    <size x="%d" y="%d" xunits="pixels" yunits="pixels"/>\n', opts.LegendSize(1), opts.LegendSize(2));
%     fprintf(fid, '  </ScreenOverlay>\n');
% end
% 
% % --- Optional sampled value placemarks ---
% if isfinite(opts.SampleStep) && opts.SampleStep>=1
%     writeSamplePlacemarks(fid, A, R, opts, crsIfMapRef(R, opts.CRS));
% end
% 
% % Close KML
% fprintf(fid, '</Document>\n</kml>\n');
% fclose(fid);
% 
% % ---------- Package to destination ----------
% cwd = pwd;
% try
%     if isKMZ
%         % Create real .zip then rename to .kmz
%         cd(tmp);
%         zipName = fullfile(tmp,'package.zip');
%         zip(zipName, {'doc.kml','files'});
%         cd(cwd);
%         if ~isempty(outDir) && exist(outDir,'dir')==0, mkdir(outDir); end
%         movefile(zipName, outFile, 'f');   % rename .zip -> .kmz
%     else
%         if ~isempty(outDir) && exist(outDir,'dir')==0, mkdir(outDir); end
%         destKML = fullfile(outDir, [outBase '.kml']);
%         destIMG = fullfile(outDir, imgName);
%         copyfile(imgPath, destIMG, 'f');
%         movefile(kmlTempPath, destKML, 'f');
%         if opts.AddLegend
%             copyfile(fullfile(legDir, legendName), fullfile(outDir, legendName), 'f');
%         end
%     end
% catch ME
%     cd(cwd);
%     rethrow(ME);
% end
% end
% 
% % ================== helpers ==================
% function opts = parseNV(nv, defaults)
%     opts = defaults;
%     if isempty(nv), return; end
%     if mod(numel(nv),2)~=0, error('Options must be name-value pairs.'); end
%     for k = 1:2:numel(nv)
%         name = nv{k};  val = nv{k+1};
%         if ~ischar(name) && ~isstring(name), error('Option names must be strings.'); end
%         name = char(name);
%         if ~isfield(opts, name), error('Unknown option: %s', name); end
%         opts.(name) = val;
%     end
% end
% 
% function [rgb, alpha, climUsed] = toUint8RGB(A, CLim, cmap, TransparentIf)
%     if ndims(A)==2
%         data = double(A);
%         finiteMask = isfinite(data);
% 
%         % Decide the stretch
%         if any(isnan(CLim))
%             if any(finiteMask(:))
%                 p = prctile(data(finiteMask), [2 98]);
%             else
%                 p = [0 1];
%             end
%         else
%             p = CLim;
%         end
%         % sanitize
%         if ~all(isfinite(p)) || numel(p)~=2
%             p = [0 1];
%         end
%         if p(1)==p(2)
%             d = max(1,abs(p(1)))*1e-6;
%             p = p + [-d d];
%         end
%         climUsed = p;
% 
%         % Build RGB + alpha
%         Inorm = mat2gray(data, p);
%         rgb   = ind2rgb(gray2ind(Inorm, size(cmap,1)), cmap);
%         if isnan(TransparentIf)
%             tmask = isnan(data);
%         else
%             tmask = isnan(data) | (data==TransparentIf);
%         end
%         alpha = uint8(255 * (~tmask));
% 
%     elseif size(A,3)==3
%         if isa(A,'uint8')
%             rgb = A;
%         else
%             mn = min(A(:)); mx = max(A(:));
%             if isfloat(A) && isfinite(mn) && mn>=0 && isfinite(mx) && mx<=1
%                 rgb = im2uint8(A);
%             else
%                 rgb = im2uint8(mat2gray(A));
%             end
%         end
%         alpha = uint8(255*ones(size(A,1), size(A,2)));
%         climUsed = [0 1]; % not used for RGB, but defined
%     else
%         error('A must be single-band or 3-band RGB.');
%     end
% end
% 
% function [ox, oy] = cornerXY(which)
%     switch upper(which)
%         case 'NE', ox = 1; oy = 1;
%         case 'NW', ox = 0; oy = 1;
%         case 'SE', ox = 1; oy = 0;
%         case 'SW', ox = 0; oy = 0;
%         otherwise, ox = 1; oy = 1;
%     end
% end
% 
% function crs = crsIfMapRef(R, CRSopt)
%     if isprop(R,'LatitudeLimits')
%         crs = []; % geographic ref
%     else
%         crs = normalizeCRS(CRSopt);
%     end
% end
% 
% function crs = normalizeCRS(CRSopt)
%     crs = [];
%     if isempty(CRSopt), return; end
%     if isnumeric(CRSopt), crs = projcrs(CRSopt); return; end
%     if isa(CRSopt,'projcrs'), crs = CRSopt; return; end
%     if ischar(CRSopt) || isstring(CRSopt)
%         f = char(CRSopt);
%         if exist(f,'file')==2
%             try
%                 info = georasterinfo(f);
%                 if isfield(info,'CoordinateReferenceSystem') && ~isempty(info.CoordinateReferenceSystem)
%                     crs = info.CoordinateReferenceSystem; return;
%                 end
%             catch, end
%             try
%                 info = geotiffinfo(f); % older MATLAB fallback
%                 if isfield(info,'PCS') && ~isempty(info.PCS)
%                     crs = projcrs(info.PCS); return;
%                 end
%             catch, end
%         end
%     end
% end
% 
% function makeLegendPNG(path, cmap, clim, label, sz)
%     % Robust legend with white background; gradient spans clim directly.
%     if numel(clim)~=2 || ~all(isfinite(clim)), clim = [0 1]; end
%     if clim(1)==clim(2)
%         d = max(1,abs(clim(1)))*1e-6;
%         clim = clim + [-d d];
%     end
% 
%     w = sz(1); h = sz(2);
%     f = figure('Visible','off','Position',[100 100 w h], 'Color','w');
%     axes('Position',[0.35 0.08 0.4 0.84]); %#ok<LAXES>
% 
%     % Make a vertical strip of data that runs exactly from clim(1) to clim(2)
%     g = linspace(clim(1), clim(2), 256)';         % 256x1
%     imagesc([0 1], [clim(1) clim(2)], g);         % CData in data units
%     set(gca,'YDir','normal','XTick',[],'Box','on');
%     colormap(cmap); caxis(clim);
%     ylabel(label, 'Interpreter','none');
% 
%     try
%         exportgraphics(f, path, 'BackgroundColor','white', 'Resolution', 150);
%     catch
%         set(f,'InvertHardcopy','off');
%         print(f, path, '-dpng','-r150');
%     end
%     close(f);
% end
% 
% function writeSamplePlacemarks(fid, A, R, opts, crs)
%     % Sample pixel centers every SampleStep and write placemarks in a Folder.
%     % - Starts near the center (avoids NaN borders).
%     % - Falls back to at least 1 valid sample if the grid hits only NaNs.
% 
%     % 1) Build the sampling grid
%     step = max(1, round(opts.SampleStep));
%     [nrows, ncols, ~] = size(A);
% 
%     % Start from ~center to avoid NaN borders
%     r0 = max(1, min(nrows, round(nrows/2)));
%     c0 = max(1, min(ncols, round(ncols/2)));
% 
%     % Expand outwards in both directions at the given step
%     rr = unique([r0:-step:1, r0:step:nrows]);
%     cc = unique([c0:-step:1, c0:step:ncols]);
%     [RR, CC] = ndgrid(rr, cc);
% 
%     % 2) Pull scalar values (if RGB, take first channel)
%     if ndims(A)==2
%         vals = A(sub2ind([nrows ncols], RR(:), CC(:)));
%     else
%         vals = A(:,:,1);
%         vals = vals(sub2ind([nrows ncols], RR(:), CC(:)));
%     end
% 
%     % Keep finite only
%     good = isfinite(vals);
%     RR = RR(good); CC = CC(good); vals = vals(good);
% 
%     % 3) If none survived (e.g., all NaN), fallback to one valid pixel
%     if isempty(vals)
%         if ndims(A)==2
%             [ii,jj] = find(isfinite(A), 1, 'first');
%         else
%             [ii,jj] = find(isfinite(A(:,:,1)), 1, 'first');
%         end
%         if isempty(ii)
%             % Truly no finite data; nothing to write
%             return
%         end
%         RR = ii; CC = jj; vals = vals; %#ok<NASGU>
%         if ndims(A)==2
%             vals = A(ii,jj);
%         else
%             vals = A(ii,jj,1);
%         end
%     end
% 
%     % 4) Cap total samples
%     if numel(vals) > opts.MaxSamples
%         idx = randperm(numel(vals), opts.MaxSamples);
%         RR = RR(idx); CC = CC(idx); vals = vals(idx);
%     end
%     npts = numel(vals);
% 
%     % 5) Pixel-center coordinates
%     if isprop(R,'LatitudeLimits')
%         [lat, lon] = intrinsicToGeographic(R, CC, RR);
%     else
%         [x, y] = intrinsicToWorld(R, CC, RR);
%         [lat, lon] = projinv(crs, x, y);
%     end
% 
%     % 6) Styles & folder
%     fprintf(fid, '  <Style id="ptStyle">\n');
%     fprintf(fid, '    <IconStyle><scale>0.9</scale></IconStyle>\n'); % default pushpin
%     fprintf(fid, '    <LabelStyle><scale>0.8</scale></LabelStyle>\n');
%     fprintf(fid, '  </Style>\n');
% 
%     fprintf(fid, '  <Folder>\n');
%     fprintf(fid, '    <name>Samples (%d)</name>\n', npts);
%     fprintf(fid, '    <visibility>1</visibility>\n');
% 
%     fmtVal = opts.ValueFormat;
%     hasUnits = ~isempty(opts.Units);
% 
%     for k = 1:npts
%         vtxt = sprintf(fmtVal, vals(k));
%         if hasUnits, vtxt = [vtxt ' ' opts.Units]; end
%         fprintf(fid, '    <Placemark>\n');
%         fprintf(fid, '      <styleUrl>#ptStyle</styleUrl>\n');
%         fprintf(fid, '      <name>%s</name>\n', vtxt);
%         fprintf(fid, '      <description><![CDATA[Value: <b>%s</b>]]></description>\n', vtxt);
%         fprintf(fid, '      <Point><coordinates>%.10f,%.10f,0</coordinates></Point>\n', lon(k), lat(k));
%         fprintf(fid, '      <ExtendedData><Data name="value"><value>%s</value></Data></ExtendedData>\n', vtxt);
%         fprintf(fid, '    </Placemark>\n');
%     end
%     fprintf(fid, '  </Folder>\n');
% end
% 
% % function writeSamplePlacemarks(fid, A, R, opts, crs)
% %     % Sample pixel centers every SampleStep and write placemarks in a Folder
% %     step = max(1, round(opts.SampleStep));
% %     [nrows, ncols, ~] = size(A);
% %     rr = 1:step:nrows; cc = 1:step:ncols;
% %     [RR, CC] = ndgrid(rr, cc);
% % 
% %     % Scalar values (if RGB, take first channel)
% %     if ndims(A)==2
% %         vals = A(sub2ind([nrows ncols], RR(:), CC(:)));
% %     else
% %         vals = A(:,:,1);
% %         vals = vals(sub2ind([nrows ncols], RR(:), CC(:)));
% %     end
% % 
% %     good = isfinite(vals);
% %     RR = RR(good); CC = CC(good); vals = vals(good);
% % 
% %     % Cap total
% %     if numel(vals) > opts.MaxSamples
% %         idx = randperm(numel(vals), opts.MaxSamples);
% %         RR = RR(idx); CC = CC(idx); vals = vals(idx);
% %     end
% %     npts = numel(vals);
% % 
% %     if npts==0, return; end
% % 
% %     % Pixel-center coords
% %     if isprop(R,'LatitudeLimits')
% %         [lat, lon] = intrinsicToGeographic(R, CC, RR);
% %     else
% %         [x, y] = intrinsicToWorld(R, CC, RR);
% %         [lat, lon] = projinv(crs, x, y);
% %     end
% % 
% %     % Styles & folder
% %     fprintf(fid, '  <Style id="ptStyle">\n');
% %     fprintf(fid, '    <IconStyle><scale>0.9</scale></IconStyle>\n'); % default pushpin, a bit smaller
% %     fprintf(fid, '    <LabelStyle><scale>0.8</scale></LabelStyle>\n');
% %     fprintf(fid, '  </Style>\n');
% % 
% %     fprintf(fid, '  <Folder>\n');
% %     fprintf(fid, '    <name>Samples (%d)</name>\n', npts);
% %     fprintf(fid, '    <visibility>1</visibility>\n');
% % 
% %     fmtVal = opts.ValueFormat;
% %     hasUnits = ~isempty(opts.Units);
% % 
% %     for k = 1:npts
% %         vtxt = sprintf(fmtVal, vals(k));
% %         if hasUnits, vtxt = [vtxt ' ' opts.Units]; end
% %         fprintf(fid, '    <Placemark>\n');
% %         fprintf(fid, '      <styleUrl>#ptStyle</styleUrl>\n');
% %         fprintf(fid, '      <name>%s</name>\n', vtxt);
% %         fprintf(fid, '      <description><![CDATA[Value: <b>%s</b>]]></description>\n', vtxt);
% %         fprintf(fid, '      <Point><coordinates>%.10f,%.10f,0</coordinates></Point>\n', lon(k), lat(k));
% %         fprintf(fid, '      <ExtendedData><Data name="value"><value>%s</value></Data></ExtendedData>\n', vtxt);
% %         fprintf(fid, '    </Placemark>\n');
% %     end
% %     fprintf(fid, '  </Folder>\n');
% % end
% 
% % function writeSamplePlacemarks(fid, A, R, opts, crs)
% %     % Sample a grid of pixel centers every SampleStep and write KML placemarks
% %     step = max(1, round(opts.SampleStep));
% %     [nrows, ncols, ~] = size(A);
% %     rr = 1:step:nrows; cc = 1:step:ncols;
% %     [RR, CC] = ndgrid(rr, cc);
% % 
% %     % Pull scalar values (if RGB, sample the first channel)
% %     if ndims(A)==2
% %         vals = A(sub2ind([nrows ncols], RR(:), CC(:)));
% %     else
% %         ch = A(:,:,1);
% %         vals = ch(sub2ind([nrows ncols], RR(:), CC(:)));
% %     end
% %     good = isfinite(vals);
% %     RR = RR(good); CC = CC(good); vals = vals(good);
% % 
% %     % Cap number of samples
% %     if numel(vals) > opts.MaxSamples
% %         idx = randperm(numel(vals), opts.MaxSamples);
% %         RR = RR(idx); CC = CC(idx); vals = vals(idx);
% %     end
% % 
% %     % Pixel-center coordinates
% %     if isprop(R,'LatitudeLimits')
% %         [lat, lon] = intrinsicToGeographic(R, CC, RR);
% %     else
% %         [x, y] = intrinsicToWorld(R, CC, RR);
% %         [lat, lon] = projinv(crs, x, y);
% %     end
% % 
% %     % Simple style for small square icons
% %     fprintf(fid, '  <Style id="ptStyle"><IconStyle><scale>0.6</scale><Icon>\n');
% %     fprintf(fid, '    <href>http://maps.google.com/mapfiles/kml/shapes/placemark_square.png</href>\n');
% %     fprintf(fid, '  </Icon></IconStyle></Style>\n');
% % 
% %     fmtVal = opts.ValueFormat;
% %     hasUnits = ~isempty(opts.Units);
% %     for k = 1:numel(vals)
% %         vtxt = sprintf(fmtVal, vals(k));
% %         if hasUnits, vtxt = [vtxt ' ' opts.Units]; end
% %         fprintf(fid, '  <Placemark>\n');
% %         fprintf(fid, '    <styleUrl>#ptStyle</styleUrl>\n');
% %         fprintf(fid, '    <name>%s</name>\n', vtxt);
% %         fprintf(fid, '    <description><![CDATA[Value: <b>%s</b>]]></description>\n', vtxt);
% %         fprintf(fid, '    <Point><coordinates>%.10f,%.10f,0</coordinates></Point>\n', lon(k), lat(k));
% %         fprintf(fid, '    <ExtendedData><Data name="value"><value>%s</value></Data></ExtendedData>\n', vtxt);
% %         fprintf(fid, '  </Placemark>\n');
% %     end
% % end
% 
% % function georaster2kmz(A, R, outFile, varargin)
% % %GEORASTER2KMZ  Export a raster + spatial reference to KML or KMZ GroundOverlay,
% % %               with optional legend (ScreenOverlay) and clickable value sampling.
% % %
% % % KMZ:
% % %   georaster2kmz(A,R,'incidence.kmz','CRS',32612,'Name','Incidence (deg)',...
% % %       'CLim',[0 60],'AddLegend',true,'LegendLabel','Incidence (deg)',...
% % %       'SampleStep',200,'ValueFormat','%.2f','Units','deg');
% % %
% % % KML:
% % %   georaster2kmz(A,R,'incidence.kml', ... same options ...);
% % %
% % % Key Name-Value options
% % %   'CRS'           : REQUIRED if R is Map*Reference (projected). EPSG (numeric),
% % %                     a projcrs, or a GeoTIFF filename with matching CRS.
% % %   'Name'          : overlay title (default: basename of outFile)
% % %   'ImageFormat'   : 'png' (default) or 'jpg'
% % %   'CLim'          : [min max] for single-band stretch (default: auto 2–98%)
% % %   'Colormap'      : Nx3 colormap for single-band (default: parula(256))
% % %   'TransparentIf' : value to make transparent (NaNs always transparent; default NaN only)
% % %   'AddLegend'     : true/false (default false) – adds a ScreenOverlay legend
% % %   'LegendLabel'   : label string (default '')
% % %   'LegendSize'    : [width height] in pixels (default [160 360])
% % %   'LegendCorner'  : 'NE'|'NW'|'SE'|'SW' (default 'NE')
% % %   'SampleStep'    : integer pixel spacing for clickable points (default Inf=off)
% % %   'MaxSamples'    : cap number of placemarks (default 3000)
% % %   'ValueFormat'   : sprintf format for value text (default '%.3g')
% % %   'Units'         : appended to values in balloons (default '')
% % %
% % % Notes
% % % - Pixel “inspection” is emulated via placemarks sampled at a grid—you control density.
% % % - Legend is a static image anchored on the screen (doesn’t scale with zoom).
% % % - For postings vs cells, extents are handled; sampling uses pixel centers.
% % 
% % % ---------- Parse options ----------
% % defaults.Name          = '';
% % defaults.ImageFormat   = 'png';
% % defaults.CLim          = [NaN NaN];
% % defaults.Colormap      = parula(256);
% % defaults.TransparentIf = NaN;
% % defaults.CRS           = [];
% % defaults.AddLegend     = false;
% % defaults.LegendLabel   = '';
% % defaults.LegendSize    = [160 360];
% % defaults.LegendCorner  = 'NE';
% % defaults.SampleStep    = Inf;      % off by default
% % defaults.MaxSamples    = 3000;
% % defaults.ValueFormat   = '%.3g';
% % defaults.Units         = '';
% % 
% % opts = parseNV(varargin, defaults);
% % fmt  = lower(char(opts.ImageFormat));
% % if ~ismember(fmt, {'png','jpg'})
% %     error('ImageFormat must be ''png'' or ''jpg''.');
% % end
% % 
% % % Output mode by extension
% % [outDir, outBase, outExt] = fileparts(outFile);
% % isKMZ = strcmpi(outExt, '.kmz');
% % isKML = strcmpi(outExt, '.kml');
% % if ~(isKMZ || isKML)
% %     error('Output must end with .kmz or .kml');
% % end
% % 
% % % ---------- Build RGB + alpha ----------
% % [rgb, alpha, climEff] = toUint8RGB(A, opts.CLim, opts.Colormap, opts.TransparentIf);
% % 
% % % ---------- Compute LatLonBox (north/south/east/west) ----------
% % if isprop(R,'LatitudeLimits') % Geographic*Reference (lat/lon)
% %     latlim = R.LatitudeLimits;  lonlim = R.LongitudeLimits;
% %     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
% %         dlat = diff(latlim) / max(R.RasterSize(1)-1, 1);
% %         dlon = diff(lonlim) / max(R.RasterSize(2)-1, 1);
% %         latlim = latlim + [-0.5 0.5]*dlat;
% %         lonlim = lonlim + [-0.5 0.5]*dlon;
% %     end
% %     north = latlim(2); south = latlim(1); east = lonlim(2); west = lonlim(1);
% % else % Map*Reference (projected) -> need CRS to projinv corners
% %     xw = R.XWorldLimits;  yw = R.YWorldLimits;
% %     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
% %         dx = diff(xw) / max(R.RasterSize(2)-1, 1);
% %         dy = diff(yw) / max(R.RasterSize(1)-1, 1);
% %         xw = xw + [-0.5 0.5]*dx;
% %         yw = yw + [-0.5 0.5]*dy;
% %     end
% %     x = [xw(1) xw(2) xw(2) xw(1)];
% %     y = [yw(2) yw(2) yw(1) yw(1)];
% %     crs = normalizeCRS(opts.CRS);
% %     if isempty(crs)
% %         error(['Projected reference detected. Provide ''CRS'' as EPSG code, projcrs, ' ...
% %                'or a GeoTIFF filename with matching CRS.']);
% %     end
% %     [lat,lon] = projinv(crs, x, y);
% %     north = max(lat); south = min(lat); east = max(lon); west = min(lon);
% % end
% % 
% % % ---------- File staging ----------
% % tmp = tempname; mkdir(tmp);
% % if isKMZ
% %     imgDir  = fullfile(tmp,'files'); mkdir(imgDir);
% %     imgName = [outBase '.' fmt];
% %     imgRel  = fullfile('files', imgName);   % relative path inside KMZ
% % else
% %     imgDir  = tmp;
% %     imgName = [outBase '.' fmt];
% %     imgRel  = imgName;
% % end
% % imgPath = fullfile(imgDir, imgName);
% % 
% % % Write overlay image
% % if strcmp(fmt,'png')
% %     imwrite(rgb, imgPath, 'Alpha', alpha);
% % else
% %     imwrite(rgb, imgPath, 'Quality', 95);
% % end
% % 
% % % Optionally make legend image (ScreenOverlay)
% % legendRel = ''; legendSize = [];
% % if opts.AddLegend
% %     if isKMZ, legDir = imgDir; else, legDir = tmp; end
% %     legendName = [outBase '_legend.png'];
% %     legendPath = fullfile(legDir, legendName);
% %     makeLegendPNG(legendPath, opts.Colormap, climEff, opts.LegendLabel, opts.LegendSize);
% %     legendRel  = fullfile(isKMZ*'files' + ~isKMZ*"", legendName); % will fix slashes below
% %     legendRel  = char(legendRel); % helper trick above returns char([]) sometimes
% %     if isempty(legendRel)
% %         legendRel = legendName;
% %     elseif contains(legendRel, char(0)) % in case of odd char math, fallback
% %         legendRel = fullfile('files', legendName);
% %     end
% %     legendSize = opts.LegendSize;
% % end
% % 
% % % ---------- Write KML (doc.kml in temp) ----------
% % kmlTempPath = fullfile(tmp, 'doc.kml');
% % overlayName = opts.Name; if isempty(overlayName), overlayName = outBase; end
% % overlayName = char(overlayName);
% % imgHref = strrep(imgRel, filesep, '/');  % forward slashes for KMZ
% % if ~isempty(legendRel), legendHref = strrep(legendRel, filesep, '/'); end
% % 
% % fid = fopen(kmlTempPath,'w'); assert(fid>0, 'Failed to create KML.');
% % fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
% % fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n');
% % fprintf(fid, '  <name>%s</name>\n', overlayName);
% % 
% % % --- Raster GroundOverlay ---
% % fprintf(fid, '  <GroundOverlay>\n');
% % fprintf(fid, '    <name>%s</name>\n', overlayName);
% % fprintf(fid, '    <Icon><href>%s</href></Icon>\n', imgHref);
% % fprintf(fid, '    <LatLonBox>\n');
% % fprintf(fid, '      <north>%.10f</north>\n', north);
% % fprintf(fid, '      <south>%.10f</south>\n', south);
% % fprintf(fid, '      <east>%.10f</east>\n',  east );
% % fprintf(fid, '      <west>%.10f</west>\n',  west );
% % fprintf(fid, '      <rotation>0</rotation>\n');
% % fprintf(fid, '    </LatLonBox>\n');
% % fprintf(fid, '  </GroundOverlay>\n');
% % 
% % % --- Optional ScreenOverlay legend ---
% % if opts.AddLegend
% %     [ox, oy] = cornerXY(opts.LegendCorner);
% %     fprintf(fid, '  <ScreenOverlay>\n');
% %     fprintf(fid, '    <name>Legend</name>\n');
% %     fprintf(fid, '    <Icon><href>%s</href></Icon>\n', legendHref);
% %     fprintf(fid, '    <overlayXY x="%.3f" y="%.3f" xunits="fraction" yunits="fraction"/>\n', ox, oy);
% %     fprintf(fid, '    <screenXY  x="%.3f" y="%.3f" xunits="fraction" yunits="fraction"/>\n', ox, oy);
% %     fprintf(fid, '    <size x="%d" y="%d" xunits="pixels" yunits="pixels"/>\n', legendSize(1), legendSize(2));
% %     fprintf(fid, '  </ScreenOverlay>\n');
% % end
% % 
% % % --- Optional sampled value placemarks ---
% % if isfinite(opts.SampleStep) && opts.SampleStep>=1
% %     writeSamplePlacemarks(fid, A, R, opts, crsIfMapRef(R, opts.CRS));
% % end
% % 
% % % Close KML
% % fprintf(fid, '</Document>\n</kml>\n');
% % fclose(fid);
% % 
% % % ---------- Package to destination ----------
% % cwd = pwd;
% % try
% %     if isKMZ
% %         % zip then rename to .kmz
% %         cd(tmp);
% %         zipName = fullfile(tmp,'package.zip');
% %         if opts.AddLegend
% %             zip(zipName, {'doc.kml','files'});  % legend sits in files/ too
% %         else
% %             % Ensure files/ exists even if legend off (it does, for main overlay)
% %             zip(zipName, {'doc.kml','files'});
% %         end
% %         cd(cwd);
% %         if ~isempty(outDir) && exist(outDir,'dir')==0, mkdir(outDir); end
% %         movefile(zipName, outFile, 'f');
% %     else
% %         if ~isempty(outDir) && exist(outDir,'dir')==0, mkdir(outDir); end
% %         destKML = fullfile(outDir, [outBase '.kml']);
% %         destIMG = fullfile(outDir, imgName);
% %         copyfile(imgPath, destIMG, 'f');
% %         movefile(kmlTempPath, destKML, 'f');
% %         if opts.AddLegend
% %             copyfile(fullfile(legDir, legendName), fullfile(outDir, legendName), 'f');
% %         end
% %     end
% % catch ME
% %     cd(cwd);
% %     rethrow(ME);
% % end
% % end
% % 
% % % ================== helpers ==================
% % function opts = parseNV(nv, defaults)
% %     opts = defaults;
% %     if isempty(nv), return; end
% %     if mod(numel(nv),2)~=0, error('Options must be name-value pairs.'); end
% %     for k = 1:2:numel(nv)
% %         name = nv{k};  val = nv{k+1};
% %         if ~ischar(name) && ~isstring(name), error('Option names must be strings.'); end
% %         name = char(name);
% %         if ~isfield(opts, name), error('Unknown option: %s', name); end
% %         opts.(name) = val;
% %     end
% % end
% % 
% % function [rgb, alpha, climUsed] = toUint8RGB(A, CLim, cmap, TransparentIf)
% %     if ndims(A)==2
% %         data = double(A);
% %         finiteMask = isfinite(data);
% % 
% %         % Decide the stretch
% %         if any(isnan(CLim))
% %             if any(finiteMask(:))
% %                 p = prctile(data(finiteMask), [2 98]);
% %             else
% %                 p = [0 1];
% %             end
% %         else
% %             p = CLim;
% %         end
% %         % sanitize in case p collapses
% %         if ~all(isfinite(p)) || numel(p)~=2
% %             p = [0 1];
% %         end
% %         if p(1)==p(2)
% %             d = max(1,abs(p(1)))*1e-6;
% %             p = p + [-d d];
% %         end
% %         climUsed = p;
% % 
% %         % Build RGB + alpha
% %         Inorm = mat2gray(data, p);
% %         rgb   = ind2rgb(gray2ind(Inorm, size(cmap,1)), cmap);
% % 
% %         if isnan(TransparentIf)
% %             tmask = isnan(data);
% %         else
% %             tmask = isnan(data) | (data==TransparentIf);
% %         end
% %         alpha = uint8(255 * (~tmask));
% % 
% %     elseif size(A,3)==3
% %         % RGB input
% %         if isa(A,'uint8')
% %             rgb = A;
% %         else
% %             mn = min(A(:)); mx = max(A(:));
% %             if isfloat(A) && isfinite(mn) && mn>=0 && isfinite(mx) && mx<=1
% %                 rgb = im2uint8(A);
% %             else
% %                 rgb = im2uint8(mat2gray(A));
% %             end
% %         end
% %         alpha = uint8(255*ones(size(A,1), size(A,2)));
% %         climUsed = [0 1]; % not used, but well-defined
% %     else
% %         error('A must be single-band or 3-band RGB.');
% %     end
% % end
% % 
% % 
% % function [ox, oy] = cornerXY(which)
% %     switch upper(which)
% %         case 'NE', ox = 1; oy = 1;
% %         case 'NW', ox = 0; oy = 1;
% %         case 'SE', ox = 1; oy = 0;
% %         case 'SW', ox = 0; oy = 0;
% %         otherwise, ox = 1; oy = 1;
% %     end
% % end
% % 
% % function crs = crsIfMapRef(R, CRSopt)
% %     if isprop(R,'LatitudeLimits')
% %         crs = []; % geographic ref
% %     else
% %         crs = normalizeCRS(CRSopt);
% %     end
% % end
% % 
% % function crs = normalizeCRS(CRSopt)
% %     crs = [];
% %     if isempty(CRSopt), return; end
% %     if isnumeric(CRSopt), crs = projcrs(CRSopt); return; end
% %     if isa(CRSopt,'projcrs'), crs = CRSopt; return; end
% %     if ischar(CRSopt) || isstring(CRSopt)
% %         f = char(CRSopt);
% %         if exist(f,'file')==2
% %             try
% %                 info = georasterinfo(f);
% %                 if isfield(info,'CoordinateReferenceSystem') && ~isempty(info.CoordinateReferenceSystem)
% %                     crs = info.CoordinateReferenceSystem; return;
% %                 end
% %             catch, end
% %             try
% %                 info = geotiffinfo(f); % older MATLAB fallback
% %                 if isfield(info,'PCS') && ~isempty(info.PCS)
% %                     crs = projcrs(info.PCS); return;
% %                 end
% %             catch, end
% %         end
% %     end
% % end
% % 
% % function makeLegendPNG(path, cmap, clim, label, sz)
% %     % Render a simple vertical legend PNG (transparent background) using MATLAB plotting.
% %     w = sz(1); h = sz(2);
% %     f = figure('Visible','off','Position',[100 100 w h], 'Color', 'none');
% %     ax = axes(f,'Position',[0.35 0.08 0.4 0.84]); %#ok<NASGU>
% %     imagesc([0 1],[clim(1) clim(2)], (0:255)'/255);
% %     set(gca,'YDir','normal','XTick',[]);
% %     colormap(cmap); caxis(clim);
% %     cb = colorbar('Location','eastoutside'); %#ok<NASGU>
% %     ylabel(label, 'Interpreter','none');
% %     try
% %         exportgraphics(f, path, 'BackgroundColor','none', 'Resolution', 150);
% %     catch
% %         set(f,'InvertHardcopy','off');
% %         print(f, path, '-dpng','-r150');
% %     end
% %     close(f);
% % end
% % 
% % function writeSamplePlacemarks(fid, A, R, opts, crs)
% %     % Sample a grid of pixel centers every SampleStep and write KML placemarks
% %     step = max(1, round(opts.SampleStep));
% %     [nrows, ncols, ~] = size(A);
% %     rr = 1:step:nrows; cc = 1:step:ncols;
% %     [RR, CC] = ndgrid(rr, cc);
% %     vals = A(sub2ind([nrows ncols], RR(:), CC(:)));
% %     good = isfinite(vals);
% %     RR = RR(good); CC = CC(good); vals = vals(good);
% % 
% %     % Cap number of samples
% %     if numel(vals) > opts.MaxSamples
% %         idx = randperm(numel(vals), opts.MaxSamples);
% %         RR = RR(idx); CC = CC(idx); vals = vals(idx);
% %     end
% % 
% %     % Map to lat/lon at pixel centers via intrinsic coordinates (col=x, row=y)
% %     if isprop(R,'LatitudeLimits')
% %         [lat, lon] = intrinsicToGeographic(R, CC, RR);
% %     else
% %         [x, y] = intrinsicToWorld(R, CC, RR);
% %         [lat, lon] = projinv(crs, x, y);
% %     end
% % 
% %     % Write a simple style for small cross icons
% %     fprintf(fid, '  <Style id="ptStyle"><IconStyle><scale>0.6</scale><Icon>\n');
% %     fprintf(fid, '    <href>http://maps.google.com/mapfiles/kml/shapes/placemark_square.png</href>\n');
% %     fprintf(fid, '  </Icon></IconStyle></Style>\n');
% % 
% %     % Placemarks
% %     fmtVal = opts.ValueFormat;
% %     hasUnits = ~isempty(opts.Units);
% %     for k = 1:numel(vals)
% %         vtxt = sprintf(fmtVal, vals(k));
% %         if hasUnits, vtxt = [vtxt ' ' opts.Units]; end
% %         fprintf(fid, '  <Placemark>\n');
% %         fprintf(fid, '    <styleUrl>#ptStyle</styleUrl>\n');
% %         fprintf(fid, '    <name>%s</name>\n', vtxt);
% %         fprintf(fid, '    <description><![CDATA[Value: <b>%s</b>]]></description>\n', vtxt);
% %         fprintf(fid, '    <Point><coordinates>%.10f,%.10f,0</coordinates></Point>\n', lon(k), lat(k));
% %         fprintf(fid, '    <ExtendedData><Data name="value"><value>%s</value></Data></ExtendedData>\n', vtxt);
% %         fprintf(fid, '  </Placemark>\n');
% %     end
% % end
% % 
% % % function georaster2kmz(A, R, outFile, varargin)
% % % %GEORASTER2KMZ  Export a raster + spatial reference to KML or KMZ GroundOverlay.
% % % %
% % % % Usage (KMZ):
% % % %   georaster2kmz(A, R, 'overlay.kmz', 'CRS', 32612, 'Name','Incidence (deg)');
% % % %
% % % % Usage (KML + image next to it):
% % % %   georaster2kmz(A, R, 'overlay.kml', 'CRS', 32612, 'Name','Incidence (deg)');
% % % %
% % % % Inputs
% % % %   A        : MxN (single-band) or MxNx3 (RGB) raster
% % % %   R        : Geographic*Reference or Map*Reference
% % % %   outFile  : output path ending with .kmz or .kml
% % % %
% % % % Name-Value options
% % % %   'CRS'           : REQUIRED if R is Map*Reference (projected). One of:
% % % %                     numeric EPSG (e.g., 32612), projcrs object, or a
% % % %                     GeoTIFF filename whose CRS matches R.
% % % %   'Name'          : overlay title (default: 'overlay')
% % % %   'ImageFormat'   : 'png' (default) or 'jpg'
% % % %   'CLim'          : [min max] for single-band contrast stretch (default: auto 2–98%)
% % % %   'Colormap'      : Nx3 colormap for single-band (default: parula(256))
% % % %   'TransparentIf' : value to make transparent (NaNs always transparent; default NaN only)
% % % 
% % % % ---------- Parse options ----------
% % % defaults.Name          = '';
% % % defaults.ImageFormat   = 'png';     % 'png' or 'jpg'
% % % defaults.CLim          = [NaN NaN];
% % % defaults.Colormap      = parula(256);
% % % defaults.TransparentIf = NaN;
% % % defaults.CRS           = [];
% % % 
% % % opts = parseNV(varargin, defaults);
% % % fmt  = lower(char(opts.ImageFormat));
% % % if ~ismember(fmt, {'png','jpg'}), error('ImageFormat must be ''png'' or ''jpg''.'); end
% % % 
% % % % Output mode by extension
% % % [outDir, outBase, outExt] = fileparts(outFile);
% % % isKMZ = strcmpi(outExt, '.kmz');
% % % isKML = strcmpi(outExt, '.kml');
% % % if ~(isKMZ || isKML)
% % %     error('Output must end with .kmz or .kml');
% % % end
% % % 
% % % % ---------- Build RGB + alpha ----------
% % % [rgb, alpha] = toUint8RGB(A, opts.CLim, opts.Colormap, opts.TransparentIf);
% % % 
% % % % ---------- Compute LatLonBox (north/south/east/west) ----------
% % % if isprop(R,'LatitudeLimits') % Geographic*Reference (lat/lon)
% % %     latlim = R.LatitudeLimits;  lonlim = R.LongitudeLimits;
% % %     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
% % %         dlat = diff(latlim) / max(R.RasterSize(1)-1, 1);
% % %         dlon = diff(lonlim) / max(R.RasterSize(2)-1, 1);
% % %         latlim = latlim + [-0.5 0.5]*dlat;
% % %         lonlim = lonlim + [-0.5 0.5]*dlon;
% % %     end
% % %     north = latlim(2); south = latlim(1); east = lonlim(2); west = lonlim(1);
% % % else % Map*Reference (projected) -> need CRS to projinv corners
% % %     xw = R.XWorldLimits;  yw = R.YWorldLimits;
% % %     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
% % %         dx = diff(xw) / max(R.RasterSize(2)-1, 1);
% % %         dy = diff(yw) / max(R.RasterSize(1)-1, 1);
% % %         xw = xw + [-0.5 0.5]*dx;
% % %         yw = yw + [-0.5 0.5]*dy;
% % %     end
% % %     x = [xw(1) xw(2) xw(2) xw(1)];
% % %     y = [yw(2) yw(2) yw(1) yw(1)];
% % %     crs = normalizeCRS(opts.CRS);
% % %     if isempty(crs)
% % %         error(['Projected reference detected. Provide ''CRS'' as EPSG code, projcrs, ' ...
% % %                'or a GeoTIFF filename with matching CRS.']);
% % %     end
% % %     [lat,lon] = projinv(crs, x, y);
% % %     north = max(lat); south = min(lat); east = max(lon); west = min(lon);
% % % end
% % % 
% % % % ---------- Write artifacts ----------
% % % tmp = tempname; mkdir(tmp);
% % % 
% % % % Image name used by the KML <href>:
% % % % - KMZ: files/<outBase>.<fmt> inside the archive
% % % % - KML: <outBase>.<fmt> next to the .kml
% % % if isKMZ
% % %     imgDir  = fullfile(tmp,'files'); mkdir(imgDir);
% % %     imgName = [outBase '.' fmt];
% % %     imgRel  = fullfile('files', imgName);   % relative path inside KMZ
% % % else
% % %     imgDir  = tmp;
% % %     imgName = [outBase '.' fmt];
% % %     imgRel  = imgName;
% % % end
% % % imgPath = fullfile(imgDir, imgName);
% % % 
% % % % Write image
% % % if strcmp(fmt,'png')
% % %     imwrite(rgb, imgPath, 'Alpha', alpha);
% % % else
% % %     imwrite(rgb, imgPath, 'Quality', 95);
% % % end
% % % 
% % % % Write KML (always as doc.kml in temp; KMZ expects that as the root file)
% % % kmlTempPath = fullfile(tmp, 'doc.kml');
% % % overlayName = opts.Name; if isempty(overlayName), overlayName = outBase; end
% % % overlayName = char(overlayName);
% % % imgHref = strrep(imgRel, filesep, '/');  % forward slashes for KMZ
% % % 
% % % fid = fopen(kmlTempPath,'w'); assert(fid>0, 'Failed to create KML.');
% % % fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
% % % fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n');
% % % fprintf(fid, '  <name>%s</name>\n', overlayName);
% % % fprintf(fid, '  <GroundOverlay>\n');
% % % fprintf(fid, '    <name>%s</name>\n', overlayName);
% % % fprintf(fid, '    <Icon><href>%s</href></Icon>\n', imgHref);
% % % fprintf(fid, '    <LatLonBox>\n');
% % % fprintf(fid, '      <north>%.10f</north>\n', north);
% % % fprintf(fid, '      <south>%.10f</south>\n', south);
% % % fprintf(fid, '      <east>%.10f</east>\n',  east );
% % % fprintf(fid, '      <west>%.10f</west>\n',  west );
% % % fprintf(fid, '      <rotation>0</rotation>\n');
% % % fprintf(fid, '    </LatLonBox>\n');
% % % fprintf(fid, '  </GroundOverlay>\n</Document>\n</kml>\n');
% % % fclose(fid);
% % % 
% % % % Package to destination
% % % cwd = pwd;
% % % try
% % %     if isKMZ
% % %         % Create a real .zip first, then rename to .kmz
% % %         cd(tmp);
% % %         zipName = fullfile(tmp,'package.zip');
% % %         zip(zipName, {'doc.kml','files'});  % produces package.zip
% % %         cd(cwd);
% % % 
% % %         if ~isempty(outDir) && exist(outDir,'dir')==0
% % %             mkdir(outDir);
% % %         end
% % %         movefile(zipName, outFile);         % rename .zip -> .kmz
% % %     else
% % %         % Plain KML: move KML + image next to each other using outBase
% % %         if ~isempty(outDir) && exist(outDir,'dir')==0
% % %             mkdir(outDir);
% % %         end
% % %         destKML = fullfile(outDir, [outBase '.kml']);
% % %         destIMG = fullfile(outDir, imgName);
% % %         copyfile(imgPath, destIMG);
% % %         movefile(kmlTempPath, destKML);
% % %     end
% % % catch ME
% % %     cd(cwd);
% % %     rethrow(ME);
% % % end
% % % end
% % % 
% % % % ================== helpers ==================
% % % function opts = parseNV(nv, defaults)
% % %     opts = defaults;
% % %     if isempty(nv), return; end
% % %     if mod(numel(nv),2)~=0, error('Options must be name-value pairs.'); end
% % %     for k = 1:2:numel(nv)
% % %         name = nv{k};  val = nv{k+1};
% % %         if ~ischar(name) && ~isstring(name), error('Option names must be strings.'); end
% % %         name = char(name);
% % %         if ~isfield(opts, name), error('Unknown option: %s', name); end
% % %         opts.(name) = val;
% % %     end
% % % end
% % % 
% % % function [rgb, alpha] = toUint8RGB(A, CLim, cmap, TransparentIf)
% % %     if ndims(A)==2
% % %         data = double(A);
% % %         finiteMask = isfinite(data);
% % %         if any(isnan(CLim))
% % %             if any(finiteMask(:)), p = prctile(data(finiteMask), [2 98]); else, p = [0 1]; end
% % %         else
% % %             p = CLim;
% % %         end
% % %         Inorm = mat2gray(data, p);
% % %         rgb   = ind2rgb(gray2ind(Inorm, size(cmap,1)), cmap);
% % %         if isnan(TransparentIf), tmask = isnan(data);
% % %         else, tmask = isnan(data) | (data==TransparentIf); end
% % %         alpha = uint8(255 * (~tmask));
% % %     elseif size(A,3)==3
% % %         if isa(A,'uint8')
% % %             rgb = A;
% % %         else
% % %             mn = min(A(:)); mx = max(A(:));
% % %             if isfloat(A) && isfinite(mn) && mn>=0 && isfinite(mx) && mx<=1
% % %                 rgb = im2uint8(A);
% % %             else
% % %                 rgb = im2uint8(mat2gray(A));
% % %             end
% % %         end
% % %         alpha = uint8(255*ones(size(A,1), size(A,2)));
% % %     else
% % %         error('A must be single-band or 3-band RGB.');
% % %     end
% % % end
% % % 
% % % function crs = normalizeCRS(CRSopt)
% % %     crs = [];
% % %     if isempty(CRSopt), return; end
% % %     if isnumeric(CRSopt), crs = projcrs(CRSopt); return; end
% % %     if isa(CRSopt,'projcrs'), crs = CRSopt; return; end
% % %     if ischar(CRSopt) || isstring(CRSopt)
% % %         f = char(CRSopt);
% % %         if exist(f,'file')==2
% % %             try
% % %                 info = georasterinfo(f);
% % %                 if isfield(info,'CoordinateReferenceSystem') && ~isempty(info.CoordinateReferenceSystem)
% % %                     crs = info.CoordinateReferenceSystem; return;
% % %                 end
% % %             catch, end
% % %             try
% % %                 info = geotiffinfo(f); % older MATLAB fallback
% % %                 if isfield(info,'PCS') && ~isempty(info.PCS)
% % %                     crs = projcrs(info.PCS); return;
% % %                 end
% % %             catch, end
% % %         end
% % %     end
% % % end
% % % % % ---------- Parse options ----------
% % % % defaults.Name          = '';
% % % % defaults.ImageFormat   = 'png';     % 'png' or 'jpg'
% % % % defaults.CLim          = [NaN NaN];
% % % % defaults.Colormap      = parula(256);
% % % % defaults.TransparentIf = NaN;
% % % % defaults.CRS           = [];
% % % % 
% % % % opts = parseNV(varargin, defaults);
% % % % fmt  = lower(char(opts.ImageFormat));
% % % % if ~ismember(fmt, {'png','jpg'}), error('ImageFormat must be ''png'' or ''jpg''.'); end
% % % % 
% % % % % Output mode by extension
% % % % [outDir, outBase, outExt] = fileparts(outFile);
% % % % isKMZ = strcmpi(outExt, '.kmz');
% % % % isKML = strcmpi(outExt, '.kml');
% % % % if ~(isKMZ || isKML)
% % % %     error('Output must end with .kmz or .kml');
% % % % end
% % % % 
% % % % % ---------- Build RGB + alpha ----------
% % % % [rgb, alpha] = toUint8RGB(A, opts.CLim, opts.Colormap, opts.TransparentIf);
% % % % 
% % % % % ---------- Compute LatLonBox (north/south/east/west) ----------
% % % % if isprop(R,'LatitudeLimits') % Geographic*Reference (lat/lon)
% % % %     latlim = R.LatitudeLimits;  lonlim = R.LongitudeLimits;
% % % %     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
% % % %         dlat = diff(latlim) / max(R.RasterSize(1)-1, 1);
% % % %         dlon = diff(lonlim) / max(R.RasterSize(2)-1, 1);
% % % %         latlim = latlim + [-0.5 0.5]*dlat;
% % % %         lonlim = lonlim + [-0.5 0.5]*dlon;
% % % %     end
% % % %     north = latlim(2); south = latlim(1); east = lonlim(2); west = lonlim(1);
% % % % else % Map*Reference (projected) -> need CRS to projinv corners
% % % %     xw = R.XWorldLimits;  yw = R.YWorldLimits;
% % % %     if isprop(R,'RasterInterpretation') && strcmpi(R.RasterInterpretation,'postings')
% % % %         dx = diff(xw) / max(R.RasterSize(2)-1, 1);
% % % %         dy = diff(yw) / max(R.RasterSize(1)-1, 1);
% % % %         xw = xw + [-0.5 0.5]*dx;
% % % %         yw = yw + [-0.5 0.5]*dy;
% % % %     end
% % % %     x = [xw(1) xw(2) xw(2) xw(1)];
% % % %     y = [yw(2) yw(2) yw(1) yw(1)];
% % % %     crs = normalizeCRS(opts.CRS);
% % % %     if isempty(crs)
% % % %         error(['Projected reference detected. Provide ''CRS'' as EPSG code, projcrs, ' ...
% % % %                'or a GeoTIFF filename with matching CRS.']);
% % % %     end
% % % %     [lat,lon] = projinv(crs, x, y);
% % % %     north = max(lat); south = min(lat); east = max(lon); west = min(lon);
% % % % end
% % % % 
% % % % % ---------- Write artifacts ----------
% % % % tmp = tempname; mkdir(tmp);
% % % % 
% % % % % Decide image placement/name used in the KML href
% % % % if isKMZ
% % % %     imgDir  = fullfile(tmp,'files');
% % % %     mkdir(imgDir);
% % % %     imgName = ['overlay.' fmt];
% % % %     imgRel  = fullfile('files', imgName);   % relative path inside KMZ
% % % % else
% % % %     imgDir  = tmp;                          % same temp folder as KML
% % % %     imgName = [outBase '.' fmt];            % same base name as .kml
% % % %     imgRel  = imgName;
% % % % end
% % % % imgPath = fullfile(imgDir, imgName);
% % % % 
% % % % % Write image
% % % % if strcmp(fmt,'png')
% % % %     imwrite(rgb, imgPath, 'Alpha', alpha);
% % % % else
% % % %     imwrite(rgb, imgPath, 'Quality', 95);
% % % % end
% % % % 
% % % % % Write KML
% % % % kmlTempPath = fullfile(tmp, 'doc.kml');
% % % % overlayName = opts.Name; if isempty(overlayName), overlayName = 'overlay'; end
% % % % overlayName = char(overlayName);
% % % % % Ensure forward slashes in href (KMZ internal zips prefer '/')
% % % % imgHref = strrep(imgRel, filesep, '/');
% % % % 
% % % % fid = fopen(kmlTempPath,'w'); assert(fid>0, 'Failed to create KML.');
% % % % fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
% % % % fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n<Document>\n');
% % % % fprintf(fid, '  <name>%s</name>\n', overlayName);
% % % % fprintf(fid, '  <GroundOverlay>\n');
% % % % fprintf(fid, '    <name>%s</name>\n', overlayName);
% % % % fprintf(fid, '    <Icon><href>%s</href></Icon>\n', imgHref);
% % % % fprintf(fid, '    <LatLonBox>\n');
% % % % fprintf(fid, '      <north>%.10f</north>\n', north);
% % % % fprintf(fid, '      <south>%.10f</south>\n', south);
% % % % fprintf(fid, '      <east>%.10f</east>\n',  east );
% % % % fprintf(fid, '      <west>%.10f</west>\n',  west );
% % % % fprintf(fid, '      <rotation>0</rotation>\n');
% % % % fprintf(fid, '    </LatLonBox>\n');
% % % % fprintf(fid, '  </GroundOverlay>\n</Document>\n</kml>\n');
% % % % fclose(fid);
% % % % 
% % % % % Package to destination
% % % % cwd = pwd;
% % % % try
% % % %     if isKMZ
% % % %         % zip as KMZ
% % % %         cd(tmp);
% % % %         zip(outFile, {'doc.kml','files'});
% % % %         cd(cwd);
% % % %     else
% % % %         % move KML + image side-by-side to requested folder
% % % %         if ~isempty(outDir) && exist(outDir,'dir')==0
% % % %             mkdir(outDir);
% % % %         end
% % % %         destKML = fullfile(outDir, [outBase '.kml']);
% % % %         destIMG = fullfile(outDir, imgName);
% % % %         copyfile(imgPath, destIMG);
% % % %         movefile(kmlTempPath, destKML); % rename doc.kml -> <outBase>.kml
% % % %     end
% % % % catch ME
% % % %     cd(cwd);
% % % %     rethrow(ME);
% % % % end
% % % % end
% % % % 
% % % % % ================== helpers ==================
% % % % function opts = parseNV(nv, defaults)
% % % %     opts = defaults;
% % % %     if isempty(nv), return; end
% % % %     if mod(numel(nv),2)~=0
% % % %         error('Options must be name-value pairs.');
% % % %     end
% % % %     for k = 1:2:numel(nv)
% % % %         name = nv{k};  val = nv{k+1};
% % % %         if ~ischar(name) && ~isstring(name)
% % % %             error('Option names must be strings.');
% % % %         end
% % % %         name = char(name);
% % % %         if ~isfield(opts, name)
% % % %             error('Unknown option: %s', name);
% % % %         end
% % % %         opts.(name) = val;
% % % %     end
% % % % end
% % % % 
% % % % function [rgb, alpha] = toUint8RGB(A, CLim, cmap, TransparentIf)
% % % %     if ndims(A)==2
% % % %         data = double(A);
% % % %         finiteMask = isfinite(data);
% % % %         if any(isnan(CLim))
% % % %             if any(finiteMask(:))
% % % %                 p = prctile(data(finiteMask), [2 98]);
% % % %             else
% % % %                 p = [0 1];
% % % %             end
% % % %         else
% % % %             p = CLim;
% % % %         end
% % % %         Inorm = mat2gray(data, p);
% % % %         rgb   = ind2rgb(gray2ind(Inorm, size(cmap,1)), cmap);
% % % % 
% % % %         if isnan(TransparentIf)
% % % %             tmask = isnan(data);
% % % %         else
% % % %             tmask = isnan(data) | (data==TransparentIf);
% % % %         end
% % % %         alpha = uint8(255 * (~tmask));
% % % %     elseif size(A,3)==3
% % % %         if isa(A,'uint8')
% % % %             rgb = A;
% % % %         else
% % % %             mn = min(A(:)); mx = max(A(:));
% % % %             if isfloat(A) && isfinite(mn) && mn>=0 && isfinite(mx) && mx<=1
% % % %                 rgb = im2uint8(A);
% % % %             else
% % % %                 rgb = im2uint8(mat2gray(A));
% % % %             end
% % % %         end
% % % %         alpha = uint8(255*ones(size(A,1), size(A,2)));
% % % %     else
% % % %         error('A must be single-band or 3-band RGB.');
% % % %     end
% % % % end
% % % % 
% % % % function crs = normalizeCRS(CRSopt)
% % % %     crs = [];
% % % %     if isempty(CRSopt), return; end
% % % %     if isnumeric(CRSopt)
% % % %         crs = projcrs(CRSopt);
% % % %         return;
% % % %     end
% % % %     if isa(CRSopt,'projcrs')
% % % %         crs = CRSopt;
% % % %         return;
% % % %     end
% % % %     if ischar(CRSopt) || isstring(CRSopt)
% % % %         f = char(CRSopt);
% % % %         if exist(f,'file')==2
% % % %             try
% % % %                 info = georasterinfo(f);
% % % %                 if isfield(info,'CoordinateReferenceSystem') && ~isempty(info.CoordinateReferenceSystem)
% % % %                     crs = info.CoordinateReferenceSystem;
% % % %                     return;
% % % %                 end
% % % %             catch
% % % %             end
% % % %             try
% % % %                 info = geotiffinfo(f); % older MATLAB fallback
% % % %                 if isfield(info,'PCS') && ~isempty(info.PCS)
% % % %                     crs = projcrs(info.PCS);
% % % %                     return;
% % % %                 end
% % % %             catch
% % % %             end
% % % %         end
% % % %     end
% % % % end

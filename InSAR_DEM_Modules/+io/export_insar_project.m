function exportInfo = export_insar_project(outDir, demData, sarData, geomData, insarData, insarClosureData, opts)
%EXPORT_INSAR_PROJECT Export MATLAB structs and GeoTIFF rasters for an InSAR project.
%
% Standardized exporter for the CAMAS-style InSAR pipeline.
%
% Required inputs:
%   outDir            output root folder
%   demData           DEM/grid struct
%   sarData           SAR acquisition struct array
%   geomData          geometry struct array
%   insarData         interferogram struct array
%   insarClosureData  closure-bias struct array (or [])
%
% Optional opts fields:
%   saveMode          'standard' (default) | 'full' | 'light'
%   exportGeotiff     true (default)
%   exportDem         true (default)
%   exportSar         true (default)
%   exportGeom        true (default)
%   exportInsar       true (default)
%   exportClosure     true (default)
%   overwrite         false (default)
%   saveV73           true (default)
%   saveSingleMats    true (default)
%   includeLonLat     false (default)
%   includePredictors false (default unless full)
%   includeComplex    false (default unless full)
%
% Notes:
% - Heavy raster products are compacted to single/logical/uint16 where appropriate.
% - GeoTIFF export uses demData.R and demData.EPSG as the master reference.
% - Complex rasters are not written directly to GeoTIFF; magnitude/phase products are exported instead.

if nargin < 7 || isempty(opts), opts = struct(); end
opts = populate_export_defaults(opts);
if ~isfield(opts, 'sarDataForNaming') || isempty(opts.sarDataForNaming)
    opts.sarDataForNaming = sarData;
end

if ~isfield(opts, 'geomDataForNaming') || isempty(opts.geomDataForNaming)
    opts.geomDataForNaming = geomData;
end

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

runName = build_export_run_name(opts, insarData);
exportRoot = fullfile(outDir, runName);

if ~exist(exportRoot, 'dir')
    mkdir(exportRoot);
end

matDir = fullfile(exportRoot, 'mats');
geoDir = fullfile(exportRoot, 'geotiff');
metaDir = fullfile(exportRoot, 'metadata');
figDir = fullfile(exportRoot, 'figures');

mkdir_if_needed(matDir);
mkdir_if_needed(geoDir);
mkdir_if_needed(metaDir);
mkdir_if_needed(figDir);


R = demData.R;
epsgCode = demData.EPSG;

% -------- compact copies for saving --------
demSave = compact_dem_data_for_save(demData, opts);
sarSave = compact_sar_data_for_save(sarData, opts);
geomSave = compact_geom_data_for_save(geomData, opts);
insarSave = compact_insar_data_for_save(insarData, opts);
closureSave = compact_closure_data_for_save(insarClosureData, opts);

% -------- save MAT files --------
matFlag = '-v7.3';
if ~opts.saveV73
    matFlag = '-v7';
end

save(fullfile(matDir, sprintf('demData_%s.mat', opts.saveMode)), 'demSave', matFlag);
save(fullfile(matDir, sprintf('sarData_%s.mat', opts.saveMode)), 'sarSave', matFlag);
save(fullfile(matDir, sprintf('geomData_%s.mat', opts.saveMode)), 'geomSave', matFlag);
save(fullfile(matDir, sprintf('insarData_%s.mat', opts.saveMode)), 'insarSave', matFlag);
save(fullfile(matDir, sprintf('insarClosureData_%s.mat', opts.saveMode)), 'closureSave', matFlag);

% -------- metadata / manifest --------
manifest = struct();
manifest.createdOn = datestr(now, 30);
manifest.saveMode = opts.saveMode;
manifest.outDir = exportRoot;
manifest.dem = build_dem_manifest(demSave);
manifest.sar = build_sar_manifest(sarSave);
manifest.geom = build_geom_manifest(geomSave);
manifest.insar = build_insar_manifest(insarSave);
manifest.closure = build_closure_manifest(closureSave);
save(fullfile(metaDir, 'export_manifest.mat'), 'manifest', matFlag);

% -------- GeoTIFF exports --------
if opts.exportGeotiff
    if opts.exportDem
        export_dem_geotiffs(fullfile(geoDir, 'dem'), demSave, R, epsgCode, opts);
    end
    if opts.exportSar
        export_sar_geotiffs(fullfile(geoDir, 'sar'), sarSave, R, epsgCode, opts);
    end
    if opts.exportGeom
        export_geom_geotiffs(fullfile(geoDir, 'geometry'), geomSave, R, epsgCode, opts);
    end
    if opts.exportInsar
        export_insar_geotiffs(fullfile(geoDir, 'insar'), insarSave, R, epsgCode, opts);
    end
    if opts.exportClosure
        export_closure_geotiffs(fullfile(geoDir, 'closure'), closureSave, R, epsgCode, opts);
    end
end

% -------- Figure exports --------
if isfield(opts,'exportFigures') && opts.exportFigures
    export_sar_figures(fullfile(figDir, 'sar'), demSave, sarSave, opts);
    export_insar_figures(fullfile(figDir, 'insar'), demSave, insarSave, opts);
    export_closure_figures(fullfile(figDir, 'closure'), demSave, closureSave, opts);
end

exportInfo = struct();
exportInfo.matDir = matDir;
exportInfo.geoDir = geoDir;
exportInfo.metaDir = metaDir;
exportInfo.manifest = manifest;
exportInfo.figDir = figDir;
end



% ========================================================================
function opts = populate_export_defaults(opts)
set_default = @(s, f, v) setfield_if_missing(s, f, v); %#ok<SFLD>
opts = set_default(opts, 'saveMode', 'standard');
opts = set_default(opts, 'exportGeotiff', true);
opts = set_default(opts, 'exportComplexSLC', true);
opts = set_default(opts, 'exportDem', true);
opts = set_default(opts, 'exportSar', true);
opts = set_default(opts, 'exportGeom', true);
opts = set_default(opts, 'exportInsar', true);
opts = set_default(opts, 'exportClosure', true);
opts = set_default(opts, 'overwrite', false);
opts = set_default(opts, 'saveV73', true);
opts = set_default(opts, 'saveSingleMats', true);
opts = set_default(opts, 'includeLonLat', false);
opts = set_default(opts, 'includePredictors', strcmpi(opts.saveMode, 'full'));
opts = set_default(opts, 'includeComplex', strcmpi(opts.saveMode, 'full'));
opts = set_default(opts, 'exportFigures', false);
opts = set_default(opts, 'figureFormat', 'png');
opts = set_default(opts, 'figureDPI', 300);
opts = set_default(opts, 'figureVisible', 'off');
opts = set_default(opts, 'useKmAxes', true);
end

function s = setfield_if_missing(s, fieldName, value)
if ~isfield(s, fieldName) || isempty(s.(fieldName))
    s.(fieldName) = value;
end
end

function mkdir_if_needed(d)
if ~exist(d, 'dir')
    mkdir(d);
end
end

% ========================================================================
function D = compact_dem_data_for_save(D, opts)
if isempty(D), return; end

if isfield(D,'dem') && ~isempty(D.dem), D.dem = single(D.dem); end
if isfield(D,'X') && ~isempty(D.X), D.X = single(D.X); end
if isfield(D,'Y') && ~isempty(D.Y), D.Y = single(D.Y); end
if isfield(D,'nanMask') && ~isempty(D.nanMask), D.nanMask = logical(D.nanMask); end
if isfield(D,'surfaceNormal') && ~isempty(D.surfaceNormal), D.surfaceNormal = single(D.surfaceNormal); end
if isfield(D,'aspect') && ~isempty(D.aspect), D.aspect = single(D.aspect); end
if isfield(D,'slope') && ~isempty(D.slope), D.slope = single(D.slope); end

if isfield(D,'lon') && ~isempty(D.lon)
    if opts.includeLonLat || strcmpi(opts.saveMode, 'full')
        D.lon = single(D.lon);
    else
        D = rmfield(D, 'lon');
    end
end
if isfield(D,'lat') && ~isempty(D.lat)
    if opts.includeLonLat || strcmpi(opts.saveMode, 'full')
        D.lat = single(D.lat);
    else
        D = rmfield(D, 'lat');
    end
end
end

% ========================================================================
function S = compact_sar_data_for_save(S, opts)
for k = 1:numel(S)

    % Keep struct schema consistent across the struct array.
    % Do not rmfield element-by-element.

    if isfield(S(k),'mli')
        S(k).mli = {};
    end

    if isfield(S(k),'pow')
        if strcmpi(opts.saveMode, 'light')
            S(k).pow = {};
        else
            S(k).pow = cast_cell_numeric(S(k).pow, 'single');
        end
    end

    if isfield(S(k),'slc')
        if strcmpi(opts.saveMode, 'light')
            S(k).slc = {};
        else
            S(k).slc = cast_cell_numeric(S(k).slc, 'single');
        end
    end

    if isfield(S(k),'amp')
        if strcmpi(opts.saveMode, 'light')
            S(k).amp = {};
        else
            S(k).amp = cast_cell_numeric(S(k).amp, 'single');
        end
    end

    if isfield(S(k),'db')
        if strcmpi(opts.saveMode, 'light')
            S(k).db = {};
        else
            S(k).db = cast_cell_numeric(S(k).db, 'single');
        end
    end
end
end

% ========================================================================
function G = compact_geom_data_for_save(G, opts)
for k = 1:numel(G)
    if isfield(G(k),'sarGeometry') && ~isempty(G(k).sarGeometry)
        for j = 1:numel(G(k).sarGeometry)
            if isempty(G(k).sarGeometry{j}), continue; end
            Gi = G(k).sarGeometry{j};
            fieldsSingle = {'Bperp','slant','incidence','slant2','incidence2'};
            for ii = 1:numel(fieldsSingle)
                f = fieldsSingle{ii};
                if isfield(Gi,f) && ~isempty(Gi.(f))
                    Gi.(f) = single(Gi.(f));
                end
            end
            if isfield(Gi,'lookMask') && ~isempty(Gi.lookMask)
                Gi.lookMask = logical(Gi.lookMask);
            end
            if isfield(Gi,'closestIndex_master') && ~isempty(Gi.closestIndex_master)
                if all(isfinite(Gi.closestIndex_master(:))) && all(Gi.closestIndex_master(:) >= 0) && all(Gi.closestIndex_master(:) <= intmax('uint16'))
                    Gi.closestIndex_master = uint16(Gi.closestIndex_master);
                else
                    Gi.closestIndex_master = single(Gi.closestIndex_master);
                end
            end
            if isfield(Gi,'closestIndex_slave') && ~isempty(Gi.closestIndex_slave)
                if all(isfinite(Gi.closestIndex_slave(:))) && all(Gi.closestIndex_slave(:) >= 0) && all(Gi.closestIndex_slave(:) <= intmax('uint16'))
                    Gi.closestIndex_slave = uint16(Gi.closestIndex_slave);
                else
                    Gi.closestIndex_slave = single(Gi.closestIndex_slave);
                end
            end
            G(k).sarGeometry{j} = Gi;
        end
    end

    if strcmpi(opts.saveMode, 'light')
        if isfield(G(k),'sarGeometry')
            G(k).sarGeometry = {};
        end
    end
end
end

% ========================================================================
function I = compact_insar_data_for_save(I, opts)
for k = 1:numel(I)
    fieldsSingle = { ...
        'phzWrapped','phzWrappedReferenced','phzUnwrapped','phzReferenced', ...
        'coherence','refScreen','refScreenWrapped', ...
        'penetration','penetrationSensitivity','penetrationDetectableHeight', ...
        'phaseNoiseStd'};

    for ii = 1:numel(fieldsSingle)
        f = fieldsSingle{ii};
        if isfield(I(k),f) && ~isempty(I(k).(f))
            I(k).(f) = single(I(k).(f));
        end
    end

    if isfield(I(k),'penetrationValidMask') && ~isempty(I(k).penetrationValidMask)
        I(k).penetrationValidMask = logical(I(k).penetrationValidMask);
    end

    if isfield(I(k),'meanPenetrationSensitivity') && ~isempty(I(k).meanPenetrationSensitivity)
        I(k).meanPenetrationSensitivity = single(I(k).meanPenetrationSensitivity);
    end
    if isfield(I(k),'meanDetectableHeight') && ~isempty(I(k).meanDetectableHeight)
        I(k).meanDetectableHeight = single(I(k).meanDetectableHeight);
    end
    if isfield(I(k),'meanBperp') && ~isempty(I(k).meanBperp)
        I(k).meanBperp = single(I(k).meanBperp);
    end
    if isfield(I(k),'minRequiredBperp') && ~isempty(I(k).minRequiredBperp)
        I(k).minRequiredBperp = single(I(k).minRequiredBperp);
    end
    if isfield(I(k),'sensitivityTarget') && ~isempty(I(k).sensitivityTarget)
        I(k).sensitivityTarget = single(I(k).sensitivityTarget);
    end
    if isfield(I(k),'validSensitivityFraction') && ~isempty(I(k).validSensitivityFraction)
        I(k).validSensitivityFraction = single(I(k).validSensitivityFraction);
    end

    % Optional heavy complex fields: clear instead of rmfield
    if isfield(I(k),'complexCoherence')
        if opts.includeComplex || strcmpi(opts.saveMode, 'full')
            if ~isempty(I(k).complexCoherence)
                I(k).complexCoherence = single(I(k).complexCoherence);
            end
        else
            I(k).complexCoherence = [];
        end
    end

    if isfield(I(k),'complexPhase')
        if opts.includeComplex || strcmpi(opts.saveMode, 'full')
            if ~isempty(I(k).complexPhase)
                I(k).complexPhase = single(I(k).complexPhase);
            end
        else
            I(k).complexPhase = [];
        end
    end

    % Optional heavy predictor field: clear instead of rmfield
    if isfield(I(k),'referencePredictors')
        if opts.includePredictors || strcmpi(opts.saveMode, 'full')
            if ~isempty(I(k).referencePredictors)
                I(k).referencePredictors = compact_reference_predictors(I(k).referencePredictors);
            end
        else
            I(k).referencePredictors = [];
        end
    end

    % Light mode: clear large raster fields but preserve schema
    if strcmpi(opts.saveMode, 'light')
        rasterFieldsToClear = { ...
            'phzWrapped','phzWrappedReferenced','phzUnwrapped','phzReferenced', ...
            'coherence','refScreen','refScreenWrapped', ...
            'penetration','penetrationValidMask','penetrationSensitivity', ...
            'penetrationDetectableHeight','phaseNoiseStd', ...
            'complexCoherence','complexPhase','referencePredictors'};

        for ii = 1:numel(rasterFieldsToClear)
            f = rasterFieldsToClear{ii};
            if isfield(I(k), f)
                I(k).(f) = [];
            end
        end
    end
end
end

function pred = compact_reference_predictors(pred)
if isempty(pred), return; end
if isfield(pred,'azCoord') && ~isempty(pred.azCoord), pred.azCoord = single(pred.azCoord); end
if isfield(pred,'rangeCoord') && ~isempty(pred.rangeCoord), pred.rangeCoord = single(pred.rangeCoord); end
if isfield(pred,'azRaw') && ~isempty(pred.azRaw), pred.azRaw = single(pred.azRaw); end
if isfield(pred,'rangeRaw') && ~isempty(pred.rangeRaw), pred.rangeRaw = single(pred.rangeRaw); end
if isfield(pred,'mask') && ~isempty(pred.mask), pred.mask = logical(pred.mask); end
end

% ========================================================================
function C = compact_closure_data_for_save(C, opts)
for k = 1:numel(C)
    if isfield(C(k),'triplet') && ~isempty(C(k).triplet)
        C(k).triplet = uint16(C(k).triplet);
    end
    if isfield(C(k),'burst') && ~isempty(C(k).burst)
        C(k).burst = uint16(C(k).burst);
    end
    if isfield(C(k),'meanBias') && ~isempty(C(k).meanBias)
        C(k).meanBias = single(C(k).meanBias);
    end
    if isfield(C(k),'stdBias') && ~isempty(C(k).stdBias)
        C(k).stdBias = single(C(k).stdBias);
    end
    if isfield(C(k),'validFraction') && ~isempty(C(k).validFraction)
        C(k).validFraction = single(C(k).validFraction);
    end
    if isfield(C(k),'closureMap')
        if strcmpi(opts.saveMode, 'light')
            C(k).closureMap = [];
        elseif ~isempty(C(k).closureMap)
            C(k).closureMap = single(C(k).closureMap);
        end
    end
end
end

% ========================================================================
function export_dem_geotiffs(outDir, D, R, epsgCode, opts)
mkdir_if_needed(outDir);
write_geotiff_if_present(fullfile(outDir, 'dem.tif'), D, 'dem', R, epsgCode, opts);
write_geotiff_if_present(fullfile(outDir, 'slope.tif'), D, 'slope', R, epsgCode, opts);
write_geotiff_if_present(fullfile(outDir, 'aspect.tif'), D, 'aspect', R, epsgCode, opts);
if isfield(D,'nanMask') && ~isempty(D.nanMask)
    write_geotiff_raster(fullfile(outDir, 'nanmask.tif'), uint8(D.nanMask), R, epsgCode, opts);
end
end

function export_sar_geotiffs(outDir, S, R, epsgCode, opts)
mkdir_if_needed(outDir);

for ii = 1:numel(S)
    nCells = 0;

    if isfield(S(ii), 'slcFiles') && ~isempty(S(ii).slcFiles)
        nCells = numel(S(ii).slcFiles);
    end

    for jj = 1:nCells
        sd = fullfile(outDir, build_sar_dir_name(S(ii), ii, jj));
        mkdir_if_needed(sd);

        % Export refreshed power as mature MLI product
        write_cell_geotiff_if_present(fullfile(sd, 'mli.tif'), S(ii), 'pow', jj, R, epsgCode, opts);
        write_cell_geotiff_if_present(fullfile(sd, 'amp.tif'), S(ii), 'amp', jj, R, epsgCode, opts);
        write_cell_geotiff_if_present(fullfile(sd, 'db.tif'),  S(ii), 'db',  jj, R, epsgCode, opts);

        % Export complex SLC as 2-band real/imag GeoTIFF
        if isfield(opts,'exportComplexSLC') && opts.exportComplexSLC && ...
           isfield(S(ii),'slc') && numel(S(ii).slc) >= jj && ~isempty(S(ii).slc{jj})

            write_complex_geotiff( ...
                fullfile(sd, 'slc_complex.tif'), ...
                S(ii).slc{jj}, R, epsgCode, opts);
        end
    end
end
end

function export_geom_geotiffs(outDir, G, R, epsgCode, opts)
mkdir_if_needed(outDir);

for k = 1:numel(G)
    if ~isfield(G(k),'sarGeometry') || isempty(G(k).sarGeometry), continue; end

    for j = 1:numel(G(k).sarGeometry)
        Gi = G(k).sarGeometry{j};
        if isempty(Gi), continue; end

        gd = fullfile(outDir, build_geom_dir_name(G(k), j, opts));
        mkdir_if_needed(gd);

        write_geotiff_if_present(fullfile(gd, 'bperp.tif'), Gi, 'Bperp', R, epsgCode, opts);
        write_geotiff_if_present(fullfile(gd, 'slant_range.tif'), Gi, 'slant', R, epsgCode, opts);
        write_geotiff_if_present(fullfile(gd, 'incidence_angle.tif'), Gi, 'incidence', R, epsgCode, opts);
        write_geotiff_if_present(fullfile(gd, 'lookmask.tif'), Gi, 'lookMask', R, epsgCode, opts, 'logicalMask');
        write_geotiff_if_present(fullfile(gd, 'slant_range_secondary.tif'), Gi, 'slant2', R, epsgCode, opts);
        write_geotiff_if_present(fullfile(gd, 'incidence_angle_secondary.tif'), Gi, 'incidence2', R, epsgCode, opts);
    end
end
end

function export_insar_geotiffs(outDir, I, R, epsgCode, opts)
mkdir_if_needed(outDir);

for k = 1:numel(I)
    id = fullfile(outDir, build_insar_pair_dir_name(I(k), k, opts));
    mkdir_if_needed(id);

    write_geotiff_if_present(fullfile(id, 'phz_wrapped.tif'), I(k), 'phzWrapped', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'phz_wrapped_referenced.tif'), I(k), 'phzWrappedReferenced', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'phz_unwrapped.tif'), I(k), 'phzUnwrapped', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'phz_referenced.tif'), I(k), 'phzReferenced', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'coherence.tif'), I(k), 'coherence', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'ref_screen.tif'), I(k), 'refScreen', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'ref_screen_wrapped.tif'), I(k), 'refScreenWrapped', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'penetration.tif'), I(k), 'penetration', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'penetration_sensitivity.tif'), I(k), 'penetrationSensitivity', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'penetration_detectable_height.tif'), I(k), 'penetrationDetectableHeight', R, epsgCode, opts);
    write_geotiff_if_present(fullfile(id, 'phase_noise_std.tif'), I(k), 'phaseNoiseStd', R, epsgCode, opts);

    if isfield(I(k),'penetrationValidMask') && ~isempty(I(k).penetrationValidMask)
        write_geotiff_raster(fullfile(id, 'penetration_valid_mask.tif'), uint8(I(k).penetrationValidMask), R, epsgCode, opts);
    end

    if isfield(I(k),'complexCoherence') && ~isempty(I(k).complexCoherence)
        cc = I(k).complexCoherence;
        write_geotiff_raster(fullfile(id, 'complex_coherence_mag.tif'), single(abs(cc)), R, epsgCode, opts);
        write_geotiff_raster(fullfile(id, 'complex_coherence_phase.tif'), single(angle(cc)), R, epsgCode, opts);
    end

    if isfield(I(k),'complexPhase') && ~isempty(I(k).complexPhase)
        cp = I(k).complexPhase;
        write_geotiff_raster(fullfile(id, 'complex_phase_angle.tif'), single(angle(cp)), R, epsgCode, opts);
    end
end
end

function export_closure_geotiffs(outDir, C, R, epsgCode, opts)
mkdir_if_needed(outDir);
for k = 1:numel(C)
    if ~isfield(C(k),'closureMap') || isempty(C(k).closureMap), continue; end
    td = fullfile(outDir, build_closure_dir_name(C(k), k, opts));
    mkdir_if_needed(td);
    write_geotiff_raster(fullfile(td, 'closure_map.tif'), single(C(k).closureMap), R, epsgCode, opts);
end
end
% ========================================================================

function write_cell_geotiff_if_present(filename, S, fieldName, jj, R, epsgCode, opts)
if ~isfield(S, fieldName), return; end

v = S.(fieldName);
if isempty(v) || ~iscell(v) || numel(v) < jj || isempty(v{jj})
    return;
end

write_geotiff_raster(filename, v{jj}, R, epsgCode, opts);
end

function write_geotiff_if_present(filename, S, fieldName, R, epsgCode, opts, mode)
if nargin < 7, mode = ''; end
if ~isfield(S, fieldName), return; end
v = S.(fieldName);
if isempty(v), return; end
switch lower(mode)
    case 'logicalmask'
        write_geotiff_raster(filename, uint8(v), R, epsgCode, opts);
    otherwise
        write_geotiff_raster(filename, v, R, epsgCode, opts);
end
end

function write_geotiff_raster(filename, A, R, epsgCode, opts)
if exist(filename, 'file') && ~opts.overwrite
    return;
end

A = prepare_raster_for_geotiff(A);

try
    if isfield(opts, 'useCompression') && opts.useCompression
        tiffTags = struct();
        tiffTags.Compression = Tiff.Compression.Deflate;
        geotiffwrite(filename, A, R, 'CoordRefSysCode', epsgCode, 'TiffTags', tiffTags);
    else
        geotiffwrite(filename, A, R, 'CoordRefSysCode', epsgCode);
    end
catch ME
    warning('Failed to write GeoTIFF: %s\n%s', filename, ME.message);
end
end

function write_complex_geotiff(filename, slc, R, epsgCode, opts)
if exist(filename, 'file') && ~opts.overwrite
    return;
end

slc = single(slc);

A = zeros([size(slc), 2], 'single');
A(:,:,1) = real(slc);
A(:,:,2) = imag(slc);

try
    if isfield(opts, 'useCompression') && opts.useCompression
        tiffTags = struct();
        tiffTags.Compression = Tiff.Compression.Deflate;
        geotiffwrite(filename, A, R, 'CoordRefSysCode', epsgCode, 'TiffTags', tiffTags);
    else
        geotiffwrite(filename, A, R, 'CoordRefSysCode', epsgCode);
    end
catch ME
    warning('Failed to write complex GeoTIFF: %s\n%s', filename, ME.message);
end
end

function A = prepare_raster_for_geotiff(A)
if islogical(A)
    A = uint8(A);
elseif isnumeric(A)
    if isa(A, 'double')
        A = single(A);
    end
end
end

% ========================================================================
function manifest = build_dem_manifest(D)
manifest = struct('fields', fieldnames(D));
end
function manifest = build_sar_manifest(S)
manifest = struct('numAcquisitions', numel(S));
end
function manifest = build_geom_manifest(G)
manifest = struct('numGeomGroups', numel(G));
end
function manifest = build_insar_manifest(I)
manifest = struct('numPairs', numel(I));
end
function manifest = build_closure_manifest(C)
manifest = struct('numTriplets', numel(C));
end

% ========================================================================
function out = cast_cell_numeric(in, targetClass)
out = in;
if ~iscell(out), return; end
for k = 1:numel(out)
    if isnumeric(out{k}) || islogical(out{k})
        switch lower(targetClass)
            case 'single'
                out{k} = single(out{k});
            case 'double'
                out{k} = double(out{k});
            case 'logical'
                out{k} = logical(out{k});
        end
    end
end
end

function s = remove_if_exists(s, fieldName)
if isfield(s, fieldName)
    s = rmfield(s, fieldName);
end
end

function s = keep_selected_fields(s, keepFields)
allFields = fieldnames(s);
rm = setdiff(allFields, keepFields);
if ~isempty(rm)
    s = rmfield(s, rm);
end
end

function v = get_field_or_default(s, fieldName, defaultValue)
if isfield(s, fieldName) && ~isempty(s.(fieldName))
    v = s.(fieldName);
else
    v = defaultValue;
end
end

function name = build_sar_dir_name(Sk, ii, jj)
if nargin < 3 || isempty(jj)
    jj = 1;
end

slcFile = '';

if isfield(Sk, 'slcFiles') && ~isempty(Sk.slcFiles) && numel(Sk.slcFiles) >= jj
    slcFile = Sk.slcFiles{jj};
end

if isempty(slcFile)
    name = sprintf('acq_%03d_%03d', ii, jj);
else
    name = build_scene_short_name(slcFile);
end
end

function name = build_geom_dir_name(geomData, j, opts)

% ===== PRIMARY SOURCE: pairList =====
if isfield(geomData,'pairList') && size(geomData.pairList,1) >= j

    pl = geomData.pairList;

    ii1 = pl(j,1);
    jj1 = pl(j,2);
    ii2 = pl(j,3);
    jj2 = pl(j,4);

    name = build_pair_name_from_indices(ii1, jj1, ii2, jj2, opts);
else
    name = sprintf('geom_%03d', j);
end
end

function name = build_insar_pair_dir_name(Ik, k, opts)

ii1 = [];
jj1 = [];
ii2 = [];
jj2 = [];

% ===== PRIMARY SOURCE: pairList =====
if isfield(opts,'geomDataForNaming') && ...
   isfield(opts.geomDataForNaming,'pairList') && ...
   size(opts.geomDataForNaming.pairList,1) >= k

    pl = opts.geomDataForNaming.pairList;

    ii1 = pl(k,1);
    jj1 = pl(k,2);
    ii2 = pl(k,3);
    jj2 = pl(k,4);
end

% ===== FALLBACK: stored lineage (future-proof) =====
if isempty(ii1) && isfield(Ik,'ii1')
    ii1 = Ik.ii1;
    jj1 = Ik.jj1;
    ii2 = Ik.ii2;
    jj2 = Ik.jj2;
end

% ===== FINAL fallback =====
if isempty(ii1) && isfield(Ik,'pair') && numel(Ik.pair) >= 2
    ii1 = Ik.pair(1); jj1 = 1;
    ii2 = Ik.pair(2); jj2 = 1;
end

if ~isempty(ii1) && ~isempty(ii2)
    name = build_pair_name_from_indices(ii1, jj1, ii2, jj2, opts);
else
    name = sprintf('pair_%03d', k);
end
end

function name = build_closure_dir_name(Ck, k, opts)

triplet = get_field_or_default(Ck, 'triplet', []);

if isempty(triplet)
    name = sprintf('triplet_%03d', k);
    return;
end

% Try to resolve each triplet entry as a scene index into sarData
parts = strings(1, numel(triplet));

for i = 1:numel(triplet)
    idx = triplet(i);

    if isfield(opts,'sarDataForNaming') && ...
       numel(opts.sarDataForNaming) >= idx && ...
       isfield(opts.sarDataForNaming(idx),'slcFiles') && ...
       ~isempty(opts.sarDataForNaming(idx).slcFiles)

        % Use first SLC file as representative scene name
        parts(i) = string(build_scene_short_name(opts.sarDataForNaming(idx).slcFiles{1}));
    else
        parts(i) = sprintf('node_%03d', idx);
    end
end

name = sanitize_name(char(strjoin(parts, '__')));
end

function name = build_pair_name_from_indices(ii1, jj1, ii2, jj2, opts)
name1 = build_scene_short_name_from_indices(ii1, jj1, opts);
name2 = build_scene_short_name_from_indices(ii2, jj2, opts);

if isempty(name1), name1 = sprintf('acq_%03d_%03d', ii1, jj1); end
if isempty(name2), name2 = sprintf('acq_%03d_%03d', ii2, jj2); end

name = sanitize_name(sprintf('%s__%s', name1, name2));
end

function shortName = build_scene_short_name_from_indices(ii, jj, opts)
shortName = '';

if ~isfield(opts,'sarDataForNaming') || isempty(opts.sarDataForNaming), return; end
if numel(opts.sarDataForNaming) < ii, return; end

S = opts.sarDataForNaming(ii);
if ~isfield(S,'slcFiles') || isempty(S.slcFiles) || numel(S.slcFiles) < jj || isempty(S.slcFiles{jj})
    return;
end

shortName = build_scene_short_name(S.slcFiles{jj});
end

function shortName = build_scene_short_name(slcFile)
[~, base, ~] = fileparts(char(slcFile));
parts = split(string(base), '_');

idxGLSAR = find(parts == "GLSAR", 1, 'first');
idxLift  = find(startsWith(parts, "lift"), 1, 'first');
idxTime  = find(~cellfun('isempty', regexp(cellstr(parts), '^\d{8}T\d{6}$', 'once')), 1, 'first');

pol = "XX";
if numel(parts) >= 3
    p1 = parts(end-2);
    p2 = parts(end-1);
    if strlength(p1) == 1 && strlength(p2) == 1
        pol = append(p1, p2);
    end
end

if ~isempty(idxGLSAR) && ~isempty(idxLift) && ~isempty(idxTime) && idxGLSAR < numel(parts)
    shortName = sprintf('%s_%s_%s_%s_%s', ...
        char(parts(idxGLSAR)), ...
        char(parts(idxGLSAR+1)), ...
        char(parts(idxLift)), ...
        char(parts(idxTime)), ...
        char(pol));
else
    shortName = char(base);
end

shortName = sanitize_name(shortName);
end

function runName = build_export_run_name(opts, insarData)
%BUILD_EXPORT_RUN_NAME Create a unique export folder name tied to the data.

% Manual override always wins
if isfield(opts, 'exportSubdir') && ~isempty(opts.exportSubdir)
    runName = sanitize_name(opts.exportSubdir);
    return;
end

timestamp = datestr(now, 'yyyymmddTHHMMSS');

if isfield(opts, 'saveMode') && ~isempty(opts.saveMode)
    modeStr = sanitize_name(opts.saveMode);
else
    modeStr = 'standard';
end

% Optional project label
if isfield(opts, 'runLabel') && ~isempty(opts.runLabel)
    labelStr = sanitize_name(opts.runLabel);
else
    labelStr = '';
end

% Build a data-tied name from the first InSAR pair if available
if ~isempty(insarData)
    try
        pairBase = build_insar_pair_dir_name(insarData(1), 1, opts);
    catch
        pairBase = 'export';
    end
else
    pairBase = 'export';
end

if isempty(labelStr)
    runName = sprintf('%s__%s__%s', pairBase, modeStr, timestamp);
else
    runName = sprintf('%s__%s__%s__%s', labelStr, pairBase, modeStr, timestamp);
end

runName = sanitize_name(runName);
end

function s = sanitize_name(s)
if isstring(s), s = char(s); end
s = regexprep(s, '[^A-Za-z0-9_\-]', '_');
end

%=======================FIGURES======================================
function hFig = plot_raster_map(demData, Z, plotOpts)

if nargin < 3, plotOpts = struct(); end
if ~isfield(plotOpts,'titleText'),      plotOpts.titleText = ''; end
if ~isfield(plotOpts,'cbarLabel'),      plotOpts.cbarLabel = ''; end
Zvalid = double(Z(isfinite(Z)));
if ~isfield(plotOpts,'climVals') || isempty(plotOpts.climVals)
    if ~isempty(Zvalid)
        if isfield(plotOpts,'phasePiLocked') && plotOpts.phasePiLocked
            q = quantile(Zvalid, [0.05 0.95]);
            maxAbsQ = max(abs(q));
            nPi = max(1, ceil(maxAbsQ / pi));
            plotOpts.climVals = [-nPi*pi, nPi*pi];
        else
            plotOpts.climVals = quantile(Zvalid, [0.05 0.95]);
        end
    else
        plotOpts.climVals = [];
    end
end
if ~isfield(plotOpts,'cmap'),           plotOpts.cmap = parula(256); end
if ~isfield(plotOpts,'visible'),        plotOpts.visible = 'on'; end
if ~isfield(plotOpts,'useKm'),          plotOpts.useKm = true; end
if ~isfield(plotOpts,'showHillshade'),  plotOpts.showHillshade = false; end
if ~isfield(plotOpts,'alphaVal'),       plotOpts.alphaVal = 0.625; end

if plotOpts.useKm
    x = double(demData.X(1,:)) ./ 1000;
    y = double(demData.Y(:,1)) ./ 1000;
    xlab = 'Easting (km)';
    ylab = 'Northing (km)';
else
    x = double(demData.X(1,:));
    y = double(demData.Y(:,1));
    xlab = 'Easting (m)';
    ylab = 'Northing (m)';
end

hFig = figure('Visible', plotOpts.visible, 'Color', 'w');
ax = axes('Parent', hFig);

if plotOpts.showHillshade && ...
   isfield(demData,'aspect') && isfield(demData,'slope') && ...
   ~isempty(demData.aspect) && ~isempty(demData.slope)

    % --- Hillshade background ---
    bg = ((cosd(double(demData.aspect) + 45) + sind(double(demData.aspect) + 45))) .* ...
         sind(2.5 .* double(demData.slope));

    imagesc(ax, x, y, bg);colormap("bone");utils.freezeColors;
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    % colormap(ax, bone(256));

    % Freeze background colors
    axes(ax);
    % utils.freezeColors;
    hold(ax, 'on');

    % --- Overlay raster ---
    Zplot = double(Z);
    alphaMask = isfinite(Zplot);

    imagesc(ax, x, y, Zplot, ...
        'AlphaData', plotOpts.alphaVal);

    colormap(ax, plotOpts.cmap);

else
    imagesc(ax, x, y, double(Z));
    axis(ax, 'image');
    set(ax, 'YDir', 'normal');
    colormap(ax, plotOpts.cmap);
end

xlabel(ax, xlab);
ylabel(ax, ylab);
title(ax, plotOpts.titleText);

hc = colorbar(ax);
ylabel(hc, plotOpts.cbarLabel, 'FontWeight','bold','FontSize',14);

% --- Safe clim handling ---
climVals = plotOpts.climVals;
if ~isempty(climVals)
    if isnumeric(climVals) && numel(climVals) == 2 && ...
            all(isfinite(climVals)) && climVals(1) < climVals(2)
        clim(ax, double(climVals(:)'));
    else
        warning('Ignoring invalid climVals. Expected [min max] with min < max.');
    end
end

set(ax, 'FontName', 'serif', 'FontWeight', 'bold', 'FontSize', 12);
daspect([1,1,1])
end

% Figure Save Helper
function save_figure(hFig, filename, opts)
if ~isfield(opts,'figureDPI') || isempty(opts.figureDPI)
    dpi = 300;
else
    dpi = opts.figureDPI;
end

exportgraphics(hFig, filename, 'Resolution', dpi);
close(hFig);
end

% Export SAR Figures
function export_sar_figures(outDir, demData, sarData, opts)
mkdir_if_needed(outDir);

for ii = 1:numel(sarData)

    if ~isfield(sarData(ii),'slcFiles') || isempty(sarData(ii).slcFiles)
        continue;
    end

    for jj = 1:numel(sarData(ii).slcFiles)

        sceneDir = fullfile(outDir, build_sar_dir_name(sarData(ii), ii, jj));
        mkdir_if_needed(sceneDir);

        % --- MLI (pow) ---
        if isfield(sarData(ii),'pow') && numel(sarData(ii).pow) >= jj && ~isempty(sarData(ii).pow{jj})
            plotOpts = get_plot_opts(opts, 'sar', 'mli');
            hFig = plot_raster_map(demData, sarData(ii).pow{jj}, plotOpts);
            save_figure(hFig, fullfile(sceneDir, ['mli.' opts.figureFormat]), opts);
        end

        % --- dB ---
        if isfield(sarData(ii),'db') && numel(sarData(ii).db) >= jj && ~isempty(sarData(ii).db{jj})
            plotOpts = get_plot_opts(opts, 'sar', 'db');
            hFig = plot_raster_map(demData, sarData(ii).db{jj}, plotOpts);
            save_figure(hFig, fullfile(sceneDir, ['db.' opts.figureFormat]), opts);
        end

    end
end
end

% Export InSAR Figures
function export_insar_figures(outDir, demData, insarData, opts)
mkdir_if_needed(outDir);

for k = 1:numel(insarData)

    pairDir = fullfile(outDir, build_insar_pair_dir_name(insarData(k), k, opts));
    mkdir_if_needed(pairDir);

    % --- Wrapped Phase ---
    if isfield(insarData(k),'phzWrapped') && ~isempty(insarData(k).phzWrapped)
        plotOpts = get_plot_opts(opts, 'insar', 'phzWrapped');
        hFig = plot_raster_map(demData, insarData(k).phzWrapped, plotOpts);
        save_figure(hFig, fullfile(pairDir, ['phz_wrapped.' opts.figureFormat]), opts);
    end

    % --- Referenced Phase ---
    if isfield(insarData(k),'phzReferenced') && ~isempty(insarData(k).phzReferenced)
        plotOpts = get_plot_opts(opts, 'insar', 'phzReferenced');
        hFig = plot_raster_map(demData, insarData(k).phzReferenced, plotOpts);
        save_figure(hFig, fullfile(pairDir, ['phz_referenced.' opts.figureFormat]), opts);
    end

    % --- Coherence ---
    if isfield(insarData(k),'coherence') && ~isempty(insarData(k).coherence)
        plotOpts = get_plot_opts(opts, 'insar', 'coherence');
        hFig = plot_raster_map(demData, insarData(k).coherence, plotOpts);
        save_figure(hFig, fullfile(pairDir, ['coherence.' opts.figureFormat]), opts);
    end

    % --- Penetration  ---
    if isfield(insarData(k),'penetrationValid') && insarData(k).penetrationValid && ...
       isfield(insarData(k),'penetration') && ~isempty(insarData(k).penetration)

        plotOpts = get_plot_opts(opts, 'insar', 'penetration');
        hFig = plot_raster_map(demData, insarData(k).penetration, plotOpts);
        save_figure(hFig, fullfile(pairDir, ['penetration.' opts.figureFormat]), opts);
    end

    % --- Phase Noise ---
    if isfield(insarData(k),'phaseNoiseStd') && ~isempty(insarData(k).phaseNoiseStd)
        plotOpts = get_plot_opts(opts, 'insar', 'phaseNoiseStd');
        hFig = plot_raster_map(demData, insarData(k).phaseNoiseStd, plotOpts);
        save_figure(hFig, fullfile(pairDir, ['phase_noise_std.' opts.figureFormat]), opts);
    end

end
end

% Export Closure Bias Figures
function export_closure_figures(outDir, demData, closureData, opts)
mkdir_if_needed(outDir);

for k = 1:numel(closureData)

    if ~isfield(closureData(k),'closureMap') || isempty(closureData(k).closureMap)
        continue;
    end

    closureDir = fullfile(outDir, build_closure_dir_name(closureData(k), k, opts));
    mkdir_if_needed(closureDir);

    plotOpts = get_plot_opts(opts, 'closure', 'closureMap');

    hFig = plot_raster_map(demData, closureData(k).closureMap, plotOpts);
    save_figure(hFig, fullfile(closureDir, ['closure_map.' opts.figureFormat]), opts);
end
end

function plotOpts = get_plot_opts(opts, familyName, productName)
plotOpts = struct();

% Common defaults
if isfield(opts, 'plotOpts') && isfield(opts.plotOpts, 'common')
    plotOpts = merge_structs(plotOpts, opts.plotOpts.common);
end

% Family/product-specific overrides
if isfield(opts, 'plotOpts') && isfield(opts.plotOpts, familyName)
    family = opts.plotOpts.(familyName);
    if isfield(family, productName)
        plotOpts = merge_structs(plotOpts, family.(productName));
    end
end

% Fallbacks
if ~isfield(plotOpts, 'visible') && isfield(opts, 'figureVisible')
    plotOpts.visible = opts.figureVisible;
end
if ~isfield(plotOpts, 'useKm') && isfield(opts, 'useKmAxes')
    plotOpts.useKm = opts.useKmAxes;
end
if ~isfield(plotOpts, 'showHillshade')
    plotOpts.showHillshade = false;
end
if ~isfield(plotOpts, 'climVals')
    plotOpts.climVals = [];
end
if ~isfield(plotOpts, 'titleText')
    plotOpts.titleText = '';
end
if ~isfield(plotOpts, 'cbarLabel')
    plotOpts.cbarLabel = '';
end
if ~isfield(plotOpts, 'cmap')
    plotOpts.cmap = parula(256);
end
if ~isfield(plotOpts, 'phasePiLocked')
    plotOpts.phasePiLocked = false;
end
end

function out = merge_structs(a, b)
out = a;
if isempty(b), return; end
f = fieldnames(b);
for k = 1:numel(f)
    out.(f{k}) = b.(f{k});
end
end
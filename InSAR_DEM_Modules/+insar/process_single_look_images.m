function sarData = process_single_look_images(sarData, demData, geomData, CR, lambda, correction, filterSize)
%PROCESS_SINGLE_LOOK_IMAGES Compute SLI products for each SLC in sarData
% with optional CR-based RCS calibration, terrain, and range corrections.
%
% Inputs:
%   sarData     : Struct array with fields slc{j}
%   geomData    : Struct with slcGeomIndex and sarGeometry
%   demData     : DEM struct with pixelArea, slope, etc.
%   lambda      : Radar wavelength [m]
%   crResult    : (optional) struct array with .Power and .RCS per CR
%   correction  : Cell array of corrections (e.g., {'gamma0','pixelarea','rangecorrect'})
%   filterSize  : Size of median filter (default = 3)
%
% Output:
%   sarData     : Updated struct with .amp, .pow, .db fields per SLC

if nargin < 7 || isempty(filterSize), filterSize = 3; end

hasCR = exist('CR','var') && ~isempty(CR);
numDirs = numel(sarData);

for ii = 1:numDirs
    for jj = 1:numel(sarData(ii).slc)
        slc = sarData(ii).slc{jj};
        slcGeom = geomData.slcGeomIndex(ii,jj);

        if isnan(slcGeom.idx)
            error('Missing geometry index for sarData(%d).slc{%d}', ii, jj);
        end

        geom = geomData.sarGeometry{slcGeom.idx};
        if strcmpi(slcGeom.pass, 'master')
            slantRange = geom.slant;
            incidence = geom.incidence;
        elseif strcmpi(slcGeom.pass, 'slave')
            slantRange = geom.slant2;
            incidence = geom.incidence2;
        else
            error('Unknown pass type: %s', slcGeom.pass);
        end

        slope = demData.slope;
        pixelArea = demData.pixelArea;

        % Try to get CR power for this burst robustly
        if hasCR
            CR_power_all = nan(numel(CR),1);
            for c = 1:numel(CR)
                if jj <= numel(CR(c).Power{ii})
                    CR_power_all(c) = CR(c).Power{ii}(jj);
                end
            end
            P = median(CR_power_all(~isnan(CR_power_all)), 'omitnan');

            a_vals = arrayfun(@(cr) ...
                conditional_a(cr), CR);

            a = median(a_vals(~isnan(a_vals)), 'omitnan');
        else
            P = [];
            a = 1;
        end

        % Compute SLI
        [amp, pow, db] = insar.compute_sli(slc, P, a, lambda, pixelArea, incidence, slope, slantRange, correction);

        % Multi-look Median filtering
        pow = insar.nanmedfilt2(pow, [filterSize, filterSize]);
        db  = insar.nanmedfilt2(db, [filterSize, filterSize]);

        sarData(ii).amp{jj} = amp;
        sarData(ii).pow{jj} = pow;
        sarData(ii).db{jj}  = db;
    end
end
end

function a = conditional_a(cr)
    if isfield(cr, 'a') && ~isempty(cr.a)
        a = cr.a;
    else
        a = NaN;
    end
end


% function sarData = process_single_look_images(sarData, demData, geomData, CR, lambda, correction, filterSize)
% %PROCESS_SINGLE_LOOK_IMAGES Compute SLI products for each SLC in sarData
% % with optional CR-based RCS calibration, terrain, and range corrections.
% %
% % Inputs:
% %   sarData     : Struct array with fields slc{j}
% %   geomData    : Struct with slcGeomIndex and sarGeometry
% %   demData     : DEM struct with pixelArea, slope, etc.
% %   lambda      : Radar wavelength [m]
% %   crResult    : (optional) struct array with .Power and .RCS per CR
% %   correction  : Cell array of corrections (e.g., {'gamma0','pixelarea','rangecorrect'})
% %   filterSize  : Size of median filter (default = 3)
% %
% % Output:
% %   sarData     : Updated struct with .amp, .pow, .db fields per SLC
% 
% if nargin < 7 || isempty(filterSize), filterSize = 3; end
% 
% hasCR = exist('CR','var') && ~isempty(CR);
% numDirs = numel(sarData);
% 
% for ii = 1:numDirs
%     for jj = 1:numel(sarData(ii).slc)
%         slc = sarData(ii).slc{jj};
%         geomIdx = geomData.slcGeomIndex(ii,jj);
% 
%         if isnan(geomIdx)
%             error('Missing geometry index for sarData(%d).slc{%d}', ii, jj);
%         end
% 
%         geom = geomData.sarGeometry{geomIdx};
%         slantRange = geom.slant;
%         incidence = geom.incidence;
%         slope = demData.slope;
%         pixelArea = demData.pixelArea;
% 
%         % Try to get CR power for this burst
%         if hasCR
%             CR_power_all = cellfun(@(p) p(jj), {CR.Power}, 'UniformOutput', true);
%             P = median(cell2mat(CR_power_all(:)));
%             RCS_vals = [CR.RCS];
%             a_vals = [CR.a];
%             a = median(a_vals(~isnan(a_vals)));
%         else
%             P = [];
%             a = 1;
%         end
% 
%         % Compute SLI
%         [amp, pow, db] = insar.compute_sli(slc, P, a, lambda, pixelArea, incidence, slope, slantRange, correction);
% 
%         % Multi-look Median filtering
%         pow = insar.nanmedfilt2(pow, [filterSize, filterSize]);
%         db  = insar.nanmedfilt2(db, [filterSize, filterSize]);
% 
%         sarData(ii).amp{jj} = amp;
%         sarData(ii).pow{jj} = pow;
%         sarData(ii).db{jj}  = db;
%     end
% end
% end

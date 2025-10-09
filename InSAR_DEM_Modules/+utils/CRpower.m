function crResult = CRpower(crFile, crList, demData, sarData, geomData, lambda, k, weightType, sigma)
%CRPOWER Computes slant-range-corrected CR power across SLCs using amplitude+distance weighting.
%
% crResult(c) fields:
%   .Name
%   .Easting, .Northing, .Elevation      % refined (combined across all SLCs)
%   .Index                               % DEM linear index nearest to refined XY
%   .Power{dir}(slc)                     % R^2-corrected power per SLC
%   .Weights{dir}{slc}.indices/.value    % per-SLC weights over k DEM pixels
%   .BestSLCPerDir(dir)                  % argmax over SLCs per dir (optional)
%   .a, .RCS                             % if available from CR file

if nargin < 6 || isempty(k), k = 5; end
if nargin < 7 || isempty(weightType), weightType = 'gaussian'; end
if nargin < 8 || isempty(sigma), sigma = 1; end

CR = io.read_corner_reflectors(crFile);
numCRs  = numel(crList);
numDirs = numel(sarData);

CR_power   = cell(numCRs, numDirs);    % each: row vector over SLCs for that dir
CR_weights = cell(numCRs, numDirs);    % each: 1xNSLC cell; each {j} has struct indices/value
CR_idxList = cell(numCRs, 1);          % k DEM candidate indices from DEM-neighborhood search

% keep per-SLC refined coordinates to combine robustly at end
CR_refined_samples = repmat(struct('x',[],'y',[],'z',[]), numCRs, 1);

for c = 1:numCRs
    [crIdxs, ~, xCR0, yCR0, zCR0] = utils.get_cr_dem_indices(CR, crList{c}, demData, k, weightType, sigma);
    CR_idxList{c} = crIdxs;

    % initialize with initial guess so we always have at least one sample
    CR_refined_samples(c).x(end+1) = xCR0;
    CR_refined_samples(c).y(end+1) = yCR0;
    CR_refined_samples(c).z(end+1) = zCR0;
end

for i = 1:numDirs
    NSLC = numel(sarData(i).slc);
    for c = 1:numCRs
        % preallocate holders for this dir
        CR_power{c,i}   = nan(1, NSLC);
        CR_weights{c,i} = cell(1, NSLC);
    end

    for j = 1:NSLC
        slc = sarData(i).slc{j};

        % geometry selection for this SLC
        slcGeom = geomData.slcGeomIndex(i,j);
        if isnan(slcGeom.idx)
            error('Missing geometry index for sarData(%d).slc{%d}', i, j);
        end
        geom = geomData.sarGeometry{slcGeom.idx};
        switch lower(slcGeom.pass)
            case 'master'
                slantRange = geom.slant;
            case 'slave'
                slantRange = geom.slant2;
            otherwise
                error('Unknown pass type: %s', slcGeom.pass);
        end

        for c = 1:numCRs
            % candidate DEM pixels near nominal CR
            cand = CR_idxList{c};
            amps = abs(slc(cand));
            if all(~isfinite(amps))
                % no usable data
                CR_power{c,i}(j) = NaN;
                CR_weights{c,i}{j} = struct('indices', cand(:), 'value', zeros(size(cand(:))));
                continue
            end

            % pick up to k strongest pixels
            Kuse = min(k, numel(cand));
            [~, maxIdxs] = maxk(amps, Kuse);
            idxs = cand(maxIdxs);
            ampWeights = amps(maxIdxs);
            ampWeights = ampWeights / max(sum(ampWeights), eps);

            % distance weighting to nominal CR (or last refined sample)
            % use most recent refined XY (robust), else nominal CR
            x0 = CR_refined_samples(c).x(end);
            y0 = CR_refined_samples(c).y(end);

            dx = demData.X(idxs) - x0;
            dy = demData.Y(idxs) - y0;
            dists = hypot(dx, dy);

            switch lower(weightType)
                case 'inverse'
                    distWeights = 1 ./ max(dists, 1e-6);
                case 'gaussian'
                    distWeights = exp(-(dists.^2) / (2*sigma^2));
                otherwise
                    error('Unknown weightType: %s', weightType);
            end

            w = ampWeights .* distWeights;
            w = w / max(sum(w), eps);

            % refined CR coordinate for this SLC
            xRef = sum(w .* demData.X(idxs));
            yRef = sum(w .* demData.Y(idxs));
            zRef = sum(w .* demData.dem(idxs));

            % store sample to be combined later
            CR_refined_samples(c).x(end+1) = xRef;
            CR_refined_samples(c).y(end+1) = yRef;
            CR_refined_samples(c).z(end+1) = zRef;

            % range-corrected power
            linPower = abs(slc(idxs)).^2 .* slantRange(idxs).^2;
            CR_power{c,i}(j) = sum(w .* linPower);

            % store weights for referencing
            CR_weights{c,i}{j} = struct('indices', idxs, 'value', w);
        end
    end
end

% Finalize CR results per corner reflector
for c = 1:numCRs
    crResult(c).Name = crList{c};

    % combine all refined coordinate samples robustly (median)
    xRefAll = CR_refined_samples(c).x(:);
    yRefAll = CR_refined_samples(c).y(:);
    zRefAll = CR_refined_samples(c).z(:);

    crResult(c).Easting  = median(xRefAll, 'omitnan');
    crResult(c).Northing = median(yRefAll, 'omitnan');
    crResult(c).Elevation= median(zRefAll, 'omitnan');

    % choose DEM index nearest to the combined refined coordinate
    cand = CR_idxList{c};
    dx = demData.X(cand) - crResult(c).Easting;
    dy = demData.Y(cand) - crResult(c).Northing;
    [~, rel] = min(dx.^2 + dy.^2);
    crResult(c).Index = cand(rel);

    % store power/weights as collected
    crResult(c).Power   = CR_power(c, :);    % cell per dir, vector per slc
    crResult(c).Weights = CR_weights(c, :);  % cell per dir, cell per slc of structs

    % best SLC per dir (optional, keeps argmax index j)
    bestPerDir = nan(1, numDirs);
    for i = 1:numDirs
        p = crResult(c).Power{i};
        if ~isempty(p) && any(isfinite(p))
            [~, bestPerDir(i)] = max(p);
        else
            bestPerDir(i) = NaN;
        end
    end
    crResult(c).BestSLCPerDir = bestPerDir;

    % add RCS if 'a' present
    crIdx = find(strcmp(CR.Name, crList{c}), 1);
    if ~isempty(crIdx) && ismember('a', CR.Properties.VariableNames)
        a = CR.a(crIdx);
        crResult(c).a   = a;
        crResult(c).RCS = (4 * pi * a^4) / (3 * lambda^2);
    else
        crResult(c).a   = NaN;
        crResult(c).RCS = NaN;
    end
end
end

% function crResult = CRpower(crFile, crList, demData, sarData, geomData, lambda, k, weightType, sigma)
% %CRPOWER Computes slant-range-corrected CR power across SLCs using amplitude-based weighting.
% 
% if nargin < 6 || isempty(k), k = 5; end
% if nargin < 7 || isempty(weightType), weightType = 'gaussian'; end
% if nargin < 8 || isempty(sigma), sigma = 1; end
% 
% CR = io.read_corner_reflectors(crFile);
% numCRs = numel(crList);
% numDirs = numel(sarData);
% 
% CR_power = cell(numCRs, numDirs);
% CR_idxList = cell(numCRs, 1);
% CR_weights = cell(numCRs, numDirs);  % Updated shape to match slc indexing
% CR_x = cell(numCRs, 1);
% CR_y = cell(numCRs, 1);
% CR_z = cell(numCRs, 1);
% 
% for c = 1:numCRs
%     [crIdxs, weights, xCR, yCR, zCR] = utils.get_cr_dem_indices(CR, crList{c}, demData, k, weightType, sigma);
%     CR_idxList{c} = crIdxs;
%     CR_x{c} = xCR;
%     CR_y{c} = yCR;
%     CR_z{c} = zCR;
% end
% 
% for i = 1:numDirs
%     for j = 1:numel(sarData(i).slc)
%         slc = sarData(i).slc{j};
%         slcGeom = geomData.slcGeomIndex(i,j);
% 
%         if isnan(slcGeom.idx)
%             error('Missing geometry index for sarData(%d).slc{%d}', i, j);
%         end
% 
%         geom = geomData.sarGeometry{slcGeom.idx};
%         if strcmpi(slcGeom.pass, 'master')
%             slantRange = geom.slant;
%         elseif strcmpi(slcGeom.pass, 'slave')
%             slantRange = geom.slant2;
%         else
%             error('Unknown pass type: %s', slcGeom.pass);
%         end
% 
%         for c = 1:numCRs
%             amps = abs(slc(CR_idxList{c}));
%             [~, maxIdxs] = maxk(amps, min(k, numel(amps)));
%             idxs = CR_idxList{c}(maxIdxs);
%             ampWeights = amps(maxIdxs);
%             ampWeights = ampWeights / sum(ampWeights);
% 
%             dists = sqrt((demData.X(idxs) - CR_x{c}).^2 + (demData.Y(idxs) - CR_y{c}).^2);
%             if strcmpi(weightType, 'inverse')
%                 distWeights = 1 ./ (dists + 1e-6);
%             else
%                 distWeights = exp(-dists.^2 / (2 * sigma^2));
%             end
% 
%             weights = ampWeights .* distWeights;
%             weights = weights / sum(weights);
% 
%             % Refine CR coordinates
%             xRefined = sum(weights .* demData.X(idxs));
%             yRefined = sum(weights .* demData.Y(idxs));
%             zRefined = sum(weights .* demData.dem(idxs));
%             CR_x{c} = xRefined;
%             CR_y{c} = yRefined;
%             CR_z{c} = zRefined;
% 
%             % Range-corrected power (multiply by R^2)
%             linPower = abs(slc(idxs)).^2 .* slantRange(idxs).^2;
%             CR_power{c, i}(j) = sum(weights .* linPower);
% 
%             CR_weights{c, i}{j}.value = weights; % Save per-SLC weights for referencing
%             CR_weights{c, i}{j}.indices = idxs;
%         end
%     end
% end
% 
% for c = 1:numCRs
%     crResult(c).Name = crList{c};
%     crResult(c).Easting = CR_x{c};
%     crResult(c).Northing = CR_y{c};
%     crResult(c).Elevation = CR_z{c};
%     [~, ix] = max(cellfun(@(x) max([x{:}.value]), CR_weights(c,:)));
%     crResult(c).Index = CR_idxList{c}(ix);
%     crResult(c).Power = CR_power(c, :);
%     crResult(c).Weights = CR_weights(c, :);
% 
%     crIdx = find(strcmp(CR.Name, crList{c}), 1);
%     if ~isempty(crIdx) && ismember('a', CR.Properties.VariableNames)
%         a = CR.a(crIdx);
%         crResult(c).a = a;
%         crResult(c).RCS = (4 * pi * a^4) / (3 * lambda^2);
%     else
%         crResult(c).a = NaN;
%         crResult(c).RCS = NaN;
%     end
% end
% 
% end
% 
% % function crResult = CRpower(crFile, crList, demData, sarData, geomData, lambda, k, weightType, sigma)
% % %CRPOWER Computes slant-range-corrected CR power across SLCs using amplitude-based weighting.
% % 
% % if nargin < 6 || isempty(k), k = 5; end
% % if nargin < 7 || isempty(weightType), weightType = 'gaussian'; end
% % if nargin < 8 || isempty(sigma), sigma = 1; end
% % 
% % CR = io.read_corner_reflectors(crFile);
% % numCRs = numel(crList);
% % numDirs = numel(sarData);
% % 
% % CR_power = cell(numCRs, numDirs);
% % CR_idxList = cell(numCRs, 1);
% % CR_weights = cell(numCRs, numDirs);  % <-- fix: use full {c,i} addressing
% % CR_x = cell(numCRs, 1);
% % CR_y = cell(numCRs, 1);
% % CR_z = cell(numCRs, 1);
% % 
% % for c = 1:numCRs
% %     [crIdxs, weights, xCR, yCR, zCR] = utils.get_cr_dem_indices(CR, crList{c}, demData, k, weightType, sigma);
% %     CR_idxList{c} = crIdxs;
% %     CR_x{c} = xCR;
% %     CR_y{c} = yCR;
% %     CR_z{c} = zCR;
% % end
% % 
% % for i = 1:numDirs
% %     for j = 1:numel(sarData(i).slc)
% %         slc = sarData(i).slc{j};
% %         slcGeom = geomData.slcGeomIndex(i,j);
% % 
% %         if isnan(slcGeom.idx)
% %             error('Missing geometry index for sarData(%d).slc{%d}', i, j);
% %         end
% % 
% %         geom = geomData.sarGeometry{slcGeom.idx};
% %         if strcmpi(slcGeom.pass, 'master')
% %             slantRange = geom.slant;
% %         elseif strcmpi(slcGeom.pass, 'slave')
% %             slantRange = geom.slant2;
% %         else
% %             error('Unknown pass type: %s', slcGeom.pass);
% %         end
% % 
% %         for c = 1:numCRs
% %             amps = abs(slc(CR_idxList{c}));
% %             [~, maxIdxs] = maxk(amps, min(k, numel(amps)));
% %             idxs = CR_idxList{c}(maxIdxs);
% %             ampWeights = amps(maxIdxs);
% %             ampWeights = ampWeights / sum(ampWeights);
% % 
% %             dists = sqrt((demData.X(idxs) - CR_x{c}).^2 + (demData.Y(idxs) - CR_y{c}).^2);
% %             if strcmpi(weightType, 'inverse')
% %                 distWeights = 1 ./ (dists + 1e-6);
% %             else
% %                 distWeights = exp(-dists.^2 / (2 * sigma^2));
% %             end
% % 
% %             weights = ampWeights .* distWeights;
% %             weights = weights / sum(weights);
% % 
% %             % Refine CR coordinates
% %             xRefined = sum(weights .* demData.X(idxs));
% %             yRefined = sum(weights .* demData.Y(idxs));
% %             zRefined = sum(weights .* demData.dem(idxs));
% %             CR_x{c} = xRefined;
% %             CR_y{c} = yRefined;
% %             CR_z{c} = zRefined;
% % 
% %             % Range-corrected power (multiply by R^2)
% %             linPower = abs(slc(idxs)).^2 .* slantRange(idxs).^2;
% %             CR_power{c, i}(j) = sum(weights .* linPower);
% % 
% %             % Assign weights with struct
% %             entry.value = weights;
% %             entry.indices = idxs;
% %             if numel(CR_weights{c, i}) < j
% %                 CR_weights{c, i}{j} = entry;
% %             else
% %                 CR_weights{c, i}{j} = entry;
% %             end
% %         end
% %     end
% % end
% % 
% % for c = 1:numCRs
% %     crResult(c).Name = crList{c};
% %     crResult(c).Easting = CR_x{c};
% %     crResult(c).Northing = CR_y{c};
% %     crResult(c).Elevation = CR_z{c};
% %     [~, ix] = max(CR_weights{c});
% %     crResult(c).Index = CR_idxList{c}(ix);
% %     crResult(c).Power = CR_power(c, :);
% %     crResult(c).Weights = CR_weights(c, :);
% % 
% %     crIdx = find(strcmp(CR.Name, crList{c}), 1);
% %     if ~isempty(crIdx) && ismember('a', CR.Properties.VariableNames)
% %         a = CR.a(crIdx);
% %         crResult(c).a = a;
% %         crResult(c).RCS = (4 * pi * a^4) / (3 * lambda^2);
% %     else
% %         crResult(c).a = NaN;
% %         crResult(c).RCS = NaN;
% %     end
% % end
% % 
% % end
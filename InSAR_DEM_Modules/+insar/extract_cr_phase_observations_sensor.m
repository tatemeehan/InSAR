function obs = extract_cr_phase_observations_sensor( ...
    phzWrapped, phzUnwrapped, cor, cphz, crResult, pred, i, b1, j, b2, opts)

imgForCR = opts.referenceImgForCR;

obs = repmat(struct( ...
    'name', '', ...
    'phiWrapped', NaN, ...
    'phiUnwrapped', NaN, ...
    'weight', NaN, ...
    'coh', NaN, ...
    'nPx', 0, ...
    'az', NaN, ...
    'rg', NaN, ...
    'supportPx', [], ...
    'supportWt', [], ...
    'isValid', false), 1, numel(crResult));

for c = 1:numel(crResult)

    switch lower(imgForCR)
        case 'first'
            wEntry = utils.pickCRWeights(crResult(c), i, b1, opts.crSlcPolicy);
        case 'second'
            wEntry = utils.pickCRWeights(crResult(c), j, b2, opts.crSlcPolicy);
        case 'both'
            e1 = utils.pickCRWeights(crResult(c), i, b1, opts.crSlcPolicy);
            e2 = utils.pickCRWeights(crResult(c), j, b2, opts.crSlcPolicy);
            wEntry = utils.combineWeightEntries(e1, e2);
        otherwise
            wEntry = utils.pickCRWeights(crResult(c), i, b1, opts.crSlcPolicy);
    end

    obs(c).name = crResult(c).Name;

    if isempty(wEntry) || ~isfield(wEntry,'indices') || ~isfield(wEntry,'value')
        continue;
    end

    px = wEntry.indices(:);
    wt = wEntry.value(:);

    good = isfinite(px) & px >= 1 & px <= numel(phzWrapped);
    px = px(good);
    wt = wt(good);

    if isempty(px)
        continue;
    end

    wt = wt / max(sum(wt), eps);

    validU = isfinite(phzUnwrapped(px));
    validW = isfinite(phzWrapped(px));
    validC = isfinite(cor(px));
    validA = isfinite(pred.azCoord(px));
    validR = isfinite(pred.rangeCoord(px));

    if ~any(validU) || ~any(validW)
        continue;
    end

    wU = wt(validU);
    wU = wU / max(sum(wU), eps);
    obs(c).phiUnwrapped = sum(wU .* phzUnwrapped(px(validU)));

    wW = wt(validW);
    wW = wW / max(sum(wW), eps);
    obs(c).phiWrapped = angle(sum(wW .* exp(1i * phzWrapped(px(validW)))));

    if any(validC)
        obs(c).coh = mean(cor(px(validC)), 'omitnan');
    end

    if any(validA)
        wA = wt(validA); wA = wA / max(sum(wA), eps);
        obs(c).az = sum(wA .* pred.azCoord(px(validA)));
    end
    if any(validR)
        wR = wt(validR); wR = wR / max(sum(wR), eps);
        obs(c).rg = sum(wR .* pred.rangeCoord(px(validR)));
    end

    obs(c).nPx = nnz(validU);
    obs(c).weight = max(obs(c).coh, 0);

    % Save support for later model scoring
    obs(c).supportPx = px;
    obs(c).supportWt = wt;

    obs(c).isValid = isfinite(obs(c).phiUnwrapped) && ...
                     isfinite(obs(c).phiWrapped)   && ...
                     isfinite(obs(c).coh)          && ...
                     obs(c).coh >= opts.referenceMinCoherence && ...
                     isfinite(obs(c).az) && isfinite(obs(c).rg);
end
end
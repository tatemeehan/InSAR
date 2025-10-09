function out = fit_sigma_gridsearch_single(gamma_meas, kz, mask, opts)
% Grid/line search for sigma_z (and optional rho) minimizing loss between
% measured coherence and MODEL * gamma_snr.
%
% Model: gamma_model = rho * exp(-0.5 * (kz*sigma)^2)
% Compared as:  yData  vs  yModelEff, where yModelEff = gamma_model .* gamma_snr
%
% Inputs
%   gamma_meas : measured coherence (0..1)
%   kz         : vertical wavenumber (1/m)
%   mask       : logical valid mask
%   opts (all optional)
%     .gamma_snr    : SNR ceiling map (same size). If empty, uses ones().
%     .floorC       : ignore points with gamma_meas < floorC  (default 0.05)
%     .useLog       : compare in log-domain (default true)
%     .loss         : 'huber' | 'rmse' (default 'huber')
%     .huberDelta   : Huber delta (default 0.05 if useLog else 0.02)
%     .weights      : extra weights map (same size) (default [])
%     .snrWeightExp : exponent p for SNR-based weighting, w *= gamma_snr^p (default 0)
%     .decimateN    : max samples used (default 2e5)
%     .binNbins     : if >0, bin in kz^2 and use bin medians (fast/robust)
%     .sigmaGrid    : vector of sigma_z candidates [m] (default linspace(0.01,0.50,60))
%     .rhoGrid      : [] (fix rho=1) OR vector in (0,1] to fit temporal/other loss
%
% Output
%   out.sigma_z, out.rho, out.bestLoss, out.lossCurve, out.sigmaGrid, out.rhoGrid, out.nUsed

    if nargin<4 || isempty(opts), opts = struct; end
    gamma_snr   = getOpt(opts,'gamma_snr',[]);
    floorC      = getOpt(opts,'floorC',0.05);
    useLog      = getOpt(opts,'useLog',false);
    lossName    = getOpt(opts,'loss','huber');
    huberDelta  = getOpt(opts,'huberDelta',[]);
    Wext        = getOpt(opts,'weights',[]);
    snrWexp     = getOpt(opts,'snrWeightExp',0);
    Nmax        = getOpt(opts,'decimateN',2e5);
    nbins       = getOpt(opts,'binNbins',0);
    sigmaGrid   = getOpt(opts,'sigmaGrid', linspace(0.01,1,100));
    rhoGrid     = getOpt(opts,'rhoGrid', []);

    if isempty(huberDelta)
        if useLog, huberDelta = 0.05; else, huberDelta = 0.02; end
    end
    if isempty(gamma_snr), gamma_snr = ones(size(gamma_meas)); end

    % -------- collect samples (no SNR division here!) --------
    valid = isfinite(gamma_meas) & isfinite(kz) & isfinite(gamma_snr) & mask & gamma_meas>0;
    y  = gamma_meas(valid);
    ks = kz(valid);
    gs = max(gamma_snr(valid), 1e-6);         % safe floor
    w  = ones(size(y));
    if ~isempty(Wext), w = w .* Wext(valid); end
    if snrWexp ~= 0,    w = w .* (gs.^snrWexp); end

    keep = (y >= floorC) & (y <= 0.999);
    y = y(keep); ks = ks(keep); gs = gs(keep); w = w(keep);

    % decimate if needed
    n = numel(y);
    if n > Nmax
        ix = randperm(n, Nmax);
        y = y(ix); ks = ks(ix); gs = gs(ix); w = w(ix);
    end

    % optional binning in kz^2
    x = ks.^2;
    if nbins>0
        edges = linspace(min(x), max(x), nbins+1);
        xc = 0.5*(edges(1:end-1)+edges(2:end));
        ym = nan(nbins,1); gm = nan(nbins,1); wm = zeros(nbins,1);
        for i=1:nbins
            sel = x>=edges(i) & x<edges(i+1);
            if nnz(sel) > 50
                ym(i) = median(y(sel),'omitnan');
                gm(i) = median(gs(sel),'omitnan');
                wm(i) = nnz(sel);
            end
        end
        ok = isfinite(ym) & isfinite(gm) & wm>0;
        x = xc(ok); y = ym(ok); gs = gm(ok); w = wm(ok)/max(wm(ok));
    end

    % set comparison domain
    if useLog
        yData = log(y);
    else
        yData = y;
    end
    x = x(:).'; w  = w(:); gs = gs(:);
    if isempty(rhoGrid), rhoGrid = 1; end
    lossCurve = nan(numel(sigmaGrid), numel(rhoGrid));

    % -------- evaluate loss over the grid (vectorized where helpful) --------
    for jr = 1:numel(rhoGrid)
        rho = rhoGrid(jr);
        for is = 1:numel(sigmaGrid)
            m = rho * exp(-0.5 * (sigmaGrid(is)^2) * x);   % 1 × Nx
            mEff = max(1e-9, m.' .* gs);                   % Nx × 1
            if useLog, yM = log(mEff); else, yM = mEff; end
            lossCurve(is,jr) = weighted_loss(yData - yM, w, lossName, huberDelta);
        end
    end
    % for jr = 1:numel(rhoGrid)
    %     rho = rhoGrid(jr);
    %     for is = 1:numel(sigmaGrid)
    %         m = rho * exp(-0.5 * (sigmaGrid(is)^2) .* x);   % model (no SNR yet)
    %         mEff = max(1e-9, m .* gs);                      % apply SNR to model
    %         if useLog
    %             yM = log(mEff);
    %         else
    %             yM = mEff;
    %         end
    %         lossCurve(is,jr) = weighted_loss(yData - yM, w, lossName, huberDelta);
    %     end
    % end

    [bestLoss, id] = min(lossCurve(:));
    [is, jr] = ind2sub(size(lossCurve), id);

    out = struct('sigma_z', sigmaGrid(is), ...
                 'rho',     rhoGrid(jr), ...
                 'bestLoss',bestLoss, ...
                 'lossCurve',lossCurve, ...
                 'sigmaGrid',sigmaGrid, ...
                 'rhoGrid',  rhoGrid, ...
                 'nUsed',    numel(y));
end

function L = weighted_loss(res, w, name, delta)
    res = res(:); w = w(:);
    switch lower(name)
        case 'rmse'
            L = sqrt( sum(w .* (res.^2)) / max(sum(w),eps) );
        case 'huber'
            a  = abs(res);
            Lh = (a<=delta) .* (0.5*res.^2) + (a>delta) .* (delta*(a - 0.5*delta));
            L  = sum(w .* Lh) / max(sum(w),eps);
        otherwise
            error('Unknown loss: %s', name);
    end
end

function v = getOpt(S,f,def)
    if isstruct(S) && isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = def; end
end

% function out = fit_sigma_gridsearch_single(gamma_meas, kz, mask, opts)
% % Grid/line search for sigma_z (and optional rho) minimizing a robust loss
% % between measured and modeled coherence.
% %
% % Inputs
% %   gamma_meas : measured coherence map (0..1)
% %   kz         : vertical wavenumber map (1/m)
% %   mask       : logical mask of valid pixels
% %   opts (all optional)
% %     .gamma_snr   : SNR ceiling map (same size); if present, gamma_meas is divided by it
% %     .floorC      : ignore points with gamma_meas < floorC (default 0.05)
% %     .useLog      : compare in log-domain (default true)
% %     .loss        : 'huber' | 'rmse' (default 'huber')
% %     .huberDelta  : Huber delta (default 0.05 if useLog else 0.02)
% %     .weights     : extra weights map (same size) (default [])
% %     .decimateN   : max samples to use (default 200e3)
% %     .binNbins    : if set (>0), bin in kz^2 first and use bin medians
% %     .sigmaGrid   : vector of sigma_z candidates [m] (default linspace(0.01,0.50,60))
% %     .rhoGrid     : [] (fix rho=1) OR vector in (0,1] if you want to fit temporal/other loss
% %
% % Output (struct)
% %   .sigma_z, .rho, .bestLoss
% %   .lossCurve       (per sigma if rho fixed, else matrix [sigma x rho])
% %   .predCurve       (model vs unique kz^2 grid if binned)
% %   .usedIdx, .nUsed
% 
% if nargin<4, opts = struct; end
% g = @(f,d) (isfield(opts,f) && ~isempty(opts.(f))) * opts.(f) + (~(isfield(opts,f) && ~isempty(opts.(f)))) * d;
% 
% gamma_snr = g('gamma_snr',[]);
% floorC    = g('floorC',0.05);
% useLog    = g('useLog',true);
% lossName  = g('loss','huber');
% hDelta    = g('huberDelta',[]);
% Wext      = g('weights',[]);
% Nmax      = g('decimateN',200e3);
% nbins     = g('binNbins',0);
% sigGrid   = g('sigmaGrid', linspace(0.01,0.50,60));       % meters
% rhoGrid   = g('rhoGrid', []);                             % [] => fix rho = 1
% 
% if isempty(hDelta)
%     if useLog
%         hDelta = 0.05;
%     else
%         hDelta = 0.02;
%     end
% end
% 
% % ----- gather samples
% valid = isfinite(gamma_meas) & isfinite(kz) & mask;
% if ~isempty(gamma_snr)
%     valid = valid & isfinite(gamma_snr) & (gamma_snr > 1e-3);
% end
% 
% gY = gamma_meas(valid);
% k  = kz(valid);
% w  = ones(size(gY));
% if ~isempty(Wext), w = w .* Wext(valid); end
% if ~isempty(gamma_snr)
%     % remove SNR ceiling (cap for safety)
%     gY = gY ./ max(gamma_snr(valid), 1e-3);
%     gY = min(gY, 0.999);
%     % make SNR a weight if you like (often helps)
%     w  = w .* (gamma_snr(valid).^2);
% end
% keep = (gY >= floorC) & (gY <= 0.999);
% gY = gY(keep); k = k(keep); w = w(keep);
% 
% % Decimate if needed
% n = numel(gY);
% if n > Nmax
%     ix = randperm(n, Nmax);
%     gY = gY(ix); k = k(ix); w = w(ix);
% end
% 
% % Optional binning in kz^2 to avoid heavy samples
% if nbins>0
%     x = (k.^2);
%     edges = linspace(min(x), max(x), nbins+1);
%     xc = 0.5*(edges(1:end-1)+edges(2:end));
%     ym = nan(nbins,1); wm = zeros(nbins,1);
%     for i=1:nbins
%         sel = x>=edges(i) & x<edges(i+1);
%         if nnz(sel) > 50
%             ym(i) = median(gY(sel),'omitnan');
%             wm(i) = nnz(sel);
%         end
%     end
%     ok  = isfinite(ym) & wm>0;
%     x   = xc(ok); gY = ym(ok); w = wm(ok)/max(wm(ok));
% else
%     x = (k.^2);
% end
% 
% % Precompute transforms
% if useLog
%     yData = log(gY);
% else
%     yData = gY;
% end
% 
% rhoList = isempty(rhoGrid) * 1 + (~isempty(rhoGrid)) * rhoGrid;  %#ok<NASGU>
% if isempty(rhoGrid), rhoGrid = 1; end
% 
% % Vectorized loss eval: loop rho, loop sigma
% lossCurve = nan(numel(sigGrid), numel(rhoGrid));
% 
% for jr = 1:numel(rhoGrid)
%     rho = rhoGrid(jr);
%     % model = rho * exp(-0.5 * sigma^2 * x)
%     for is = 1:numel(sigGrid)
%         m = rho * exp(-0.5 * (sigGrid(is)^2) .* x);
%         if useLog
%             yM = log(m);
%         else
%             yM = m;
%         end
%         lossCurve(is,jr) = weighted_loss(yData - yM, w, lossName, hDelta);
%     end
% end
% 
% % Pick argmin
% [bestLoss, idx] = min(lossCurve(:));
% [is, jr] = ind2sub(size(lossCurve), idx);
% sigma_z = sigGrid(is);
% rho     = rhoGrid(jr);
% 
% out = struct('sigma_z',sigma_z,'rho',rho,'bestLoss',bestLoss, ...
%              'lossCurve',lossCurve,'sigmaGrid',sigGrid,'rhoGrid',rhoGrid, ...
%              'nUsed',numel(yData));
% end
% 
% function L = weighted_loss(res, w, name, delta)
%     res = res(:); w = w(:);
%     switch lower(name)
%         case 'rmse'
%             L = sqrt( sum(w .* (res.^2)) / max(sum(w),eps) );
%         case 'huber'
%             a = abs(res);
%             Lh = (a<=delta) .* (0.5 * res.^2) + (a>delta) .* (delta*(a - 0.5*delta));
%             L = sum(w .* Lh) / max(sum(w),eps);
%         otherwise
%             error('Unknown loss: %s', name);
%     end
% end

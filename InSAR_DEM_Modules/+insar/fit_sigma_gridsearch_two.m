function out = fit_sigma_gridsearch_two(gamma_meas, kz, mask, opts)
% Grid search for two-component model:
%   gamma_model = rho * [ f*exp(-0.5*(kz*s1)^2) + (1-f)*exp(-0.5*(kz*s2)^2) ]
% Compare against measured as MODEL*gamma_snr (no division of measurements).

    if nargin<4 || isempty(opts), opts = struct; end
    % shared options
    gamma_snr   = getOpt(opts,'gamma_snr',[]);
    floorC      = getOpt(opts,'floorC',0.05);
    useLog      = getOpt(opts,'useLog',false);
    lossName    = getOpt(opts,'loss','huber');
    huberDelta  = getOpt(opts,'huberDelta',[]);
    Wext        = getOpt(opts,'weights',[]);
    snrWexp     = getOpt(opts,'snrWeightExp',0);
    Nmax        = getOpt(opts,'decimateN',2e5);
    nbins       = getOpt(opts,'binNbins',0);
    % two-comp grids
    s1Grid      = getOpt(opts,'sigma1Grid', linspace(0.005,0.5,51));
    s2Grid      = getOpt(opts,'sigma2Grid', linspace(0.005,0.05,11));
    fGrid       = getOpt(opts,'fracGrid',   linspace(0.05,0.5,11));
    rhoGrid     = getOpt(opts,'rhoGrid',    1);

    if isempty(huberDelta)
        if useLog, huberDelta = 0.05; else, huberDelta = 0.02; end
    end
    if isempty(gamma_snr), gamma_snr = ones(size(gamma_meas)); end

    % -------- collect samples (no SNR division!) --------
    valid = isfinite(gamma_meas) & isfinite(kz) & isfinite(gamma_snr) & mask & gamma_meas>0;
    y  = gamma_meas(valid);
    ks = kz(valid);
    gs = max(gamma_snr(valid), 1e-6);
    w  = ones(size(y));
    if ~isempty(Wext), w = w .* Wext(valid); end
    if snrWexp ~= 0,    w = w .* (gs.^snrWexp); end

    keep = (y >= floorC) & (y <= 0.999);
    y = y(keep); ks = ks(keep); gs = gs(keep); w = w(keep);

    n = numel(y);
    if n > Nmax
        ix = randperm(n, Nmax);
        y = y(ix); ks = ks(ix); gs = gs(ix); w = w(ix);
    end

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

    x = x(:).';    % 1 × Nx
    w  = w(:);               % Nx×1
    gs = gs(:);              % Nx×1

    if useLog
        yData = log(y);
    else
        yData = y;
    end

    % % precompute exponentials (speed)
    % E1 = exp(-0.5 * (x .* reshape(s1Grid.^2,[],1)));   % [Ns1 x Nx]
    % E2 = exp(-0.5 * (x .* reshape(s2Grid.^2,[],1)));   % [Ns2 x Nx]

    % --- precompute exponentials (outer products) ---
    % E1(i,n) = exp(-0.5 * s1Grid(i)^2 * x(n))
    % E2(j,n) = exp(-0.5 * s2Grid(j)^2 * x(n))
    s1sq = (s1Grid(:).^2);        % Ns1 × 1
    s2sq = (s2Grid(:).^2);        % Ns2 × 1

    E1 = exp(-0.5 * (s1sq * x));  % Ns1 × Nx  (outer product)
    E2 = exp(-0.5 * (s2sq * x));  % Ns2 × Nx

    best = struct('L',inf,'s1',NaN,'s2',NaN,'f',NaN,'rho',NaN);

    for ir = 1:numel(rhoGrid)
        rho = rhoGrid(ir);
        for i1 = 1:numel(s1Grid)
            e1 = E1(i1,:);                % 1 × Nx
            for i2 = 1:numel(s2Grid)
                if s2Grid(i2) < s1Grid(i1), continue; end
                e2 = E2(i2,:);            % 1 × Nx
                for jf = 1:numel(fGrid)
                    f = fGrid(jf);
                    m = rho * ( f*e1 + (1-f)*e2 );   % 1 × Nx
                    mEff = max(1e-9, m.' .* gs);     % Nx × 1  (apply SNR to model)
                    if useLog
                        yM = log(mEff);
                    else
                        yM = mEff;
                    end
                    lambda_sigma = 1e-3;  sigma_ref = 0.05;   % meters
                    lambda_frac  = 1e-3;

                    reg = lambda_sigma.*(s2Grid(i2)./sigma_ref).^2 + lambda_frac.*(1 - f).^2;
                    L = weighted_loss(yData - yM, w, lossName, huberDelta) + reg;
                    % L = weighted_loss(yData - yM, w, lossName, huberDelta);
                    if L < best.L
                        best.L = L; best.s1 = s1Grid(i1); best.s2 = s2Grid(i2);
                        best.f = f; best.rho = rho;
                    end
                end
            end
        end
    end
    % for ir = 1:numel(rhoGrid)
    %     rho = rhoGrid(ir);
    %     for i1 = 1:numel(s1Grid)
    %         e1 = E1(i1,:);   % 1 x Nx
    %         for i2 = 1:numel(s2Grid)
    %             if s2Grid(i2) < s1Grid(i1), continue; end
    %             e2 = E2(i2,:);
    %             for jf = 1:numel(fGrid)
    %                 f = fGrid(jf);
    %                 m = rho * ( f*e1 + (1-f)*e2 );    % model (no SNR yet)
    %                 mEff = max(1e-9, m .* gs');       % apply SNR to model; gs is Nx x 1
    %                 if useLog
    %                     yM = log(mEff);
    %                 else
    %                     yM = mEff;
    %                 end
    %                 L = weighted_loss(yData - yM', w, lossName, huberDelta);
    %                 if L < best.L
    %                     best.L = L; best.s1 = s1Grid(i1); best.s2 = s2Grid(i2);
    %                     best.f = f; best.rho = rho;
    %                 end
    %             end
    %         end
    %     end
    % end

    out = struct('sigma1',best.s1,'sigma2',best.s2,'frac',best.f,'rho',best.rho, ...
                 'bestLoss',best.L, 'nUsed', numel(y));
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

% function out = fit_sigma_gridsearch_two(gamma_meas, kz, mask, opts)
% % Grid search for two-component model:
% % gamma_model = rho * [ f*exp(-0.5*(kz*s1)^2) + (1-f)*exp(-0.5*(kz*s2)^2) ]
% % with constraints 0 < s1 <= s2, 0 <= f <= 1.
% %
% % Key opts (others same as single version):
% %   .sigma1Grid = linspace(0.02,0.25,30)
% %   .sigma2Grid = linspace(0.10,0.80,30)
% %   .fracGrid   = linspace(0.1,0.9,9)
% %   .rhoGrid    = [] or linspace(0.7,1,7)
% 
% if nargin<4, opts = struct; end
% % reuse most options from single
% base = fit_sigma_gridsearch_single(ones(2), ones(2), true(2), struct()); %#ok<NASGU>
% get = @(f,d) (isfield(opts,f) && ~isempty(opts.(f))) * opts.(f) + (~(isfield(opts,f) && ~isempty(opts.(f)))) * d;
% 
% optsSingle = opts;  % will share parser
% optsSingle.sigmaGrid = []; % not used
% 
% % collect samples (shared routine)
% tmp = fit_sigma_gridsearch_single(gamma_meas, kz, mask, rmfield_if(optsSingle,{'sigmaGrid','rhoGrid'}));
% % Rebuild x,y,w with the same rules used inside:
% gamma_snr = get('gamma_snr',[]);
% floorC    = get('floorC',0.05);
% useLog    = get('useLog',true);
% lossName  = get('loss','huber');
% hDelta    = get('huberDelta',useLog*0.05 + (~useLog)*0.02);
% Wext      = get('weights',[]);
% Nmax      = get('decimateN',200e3);
% nbins     = get('binNbins',0);
% 
% valid = isfinite(gamma_meas) & isfinite(kz) & mask;
% if ~isempty(gamma_snr), valid = valid & isfinite(gamma_snr) & (gamma_snr>1e-3); end
% gY = gamma_meas(valid);
% k  = kz(valid);
% w  = ones(size(gY));
% if ~isempty(Wext), w = w .* Wext(valid); end
% if ~isempty(gamma_snr)
%     gY = gY ./ max(gamma_snr(valid),1e-3);
%     gY = min(gY,0.999);
%     w  = w .* (gamma_snr(valid).^2);
% end
% keep = (gY >= floorC) & (gY <= 0.999);
% gY = gY(keep); k = k(keep); w = w(keep);
% 
% if numel(gY) > Nmax
%     ix = randperm(numel(gY), Nmax);
%     gY = gY(ix); k = k(ix); w = w(ix);
% end
% 
% x = (k.^2);
% if nbins>0
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
%     ok = isfinite(ym) & wm>0;
%     x  = xc(ok); gY = ym(ok); w = wm(ok)/max(wm(ok));
% end
% 
% if useLog
%     yData = log(gY);
% else
%     yData = gY;
% end
% 
% % parameter grids
% s1Grid = get('sigma1Grid', linspace(0.02,0.25,30));
% s2Grid = get('sigma2Grid', linspace(0.10,0.80,30));
% fGrid  = get('fracGrid',   linspace(0.20,0.80,9));
% rhoGrid = get('rhoGrid', 1);
% 
% % precompute exp terms for speed
% E1 = exp(-0.5 * (x .* reshape(s1Grid.^2,[],1)));  % [Ns1 x Nx]
% E2 = exp(-0.5 * (x .* reshape(s2Grid.^2,[],1)));  % [Ns2 x Nx]
% 
% best = struct('L',inf,'s1',NaN,'s2',NaN,'f',NaN,'rho',NaN);
% 
% for ir = 1:numel(rhoGrid)
%     rho = rhoGrid(ir);
%     for i1 = 1:numel(s1Grid)
%         e1 = E1(i1,:);                    % 1 x Nx
%         for i2 = 1:numel(s2Grid)
%             if s2Grid(i2) < s1Grid(i1), continue; end
%             e2 = E2(i2,:);
%             % mix across f (vectorized)
%             % model per f: m(f) = rho * ( f*e1 + (1-f)*e2 )
%             for jf = 1:numel(fGrid)
%                 f = fGrid(jf);
%                 m = rho * ( f*e1 + (1-f)*e2 );
%                 if useLog
%                     yM = log(m);
%                 else
%                     yM = m;
%                 end                
%                 L  = weighted_loss(yData - yM.', w, lossName, hDelta);
%                 if L < best.L
%                     best.L = L; best.s1 = s1Grid(i1); best.s2 = s2Grid(i2);
%                     best.f = f; best.rho = rho;
%                 end
%             end
%         end
%     end
% end
% 
% out = struct('sigma1',best.s1,'sigma2',best.s2,'frac',best.f,'rho',best.rho, ...
%              'bestLoss',best.L,'nUsed',numel(yData));
% end
% 
% function S = rmfield_if(S,fields)
% for i=1:numel(fields)
%     if isfield(S,fields{i}), S = rmfield(S,fields{i}); end
% end
% end
% 
% function L = weighted_loss(res, w, name, delta)
% res = res(:); w = w(:);
% switch lower(name)
%     case 'rmse'
%         L = sqrt( sum(w .* (res.^2)) / max(sum(w),eps) );
%     case 'huber'
%         a = abs(res);
%         Lh = (a<=delta) .* (0.5 * res.^2) + (a>delta) .* (delta*(a - 0.5*delta));
%         L = sum(w .* Lh) / max(sum(w),eps);
%     otherwise
%         error('Unknown loss: %s', name);
% end
% end

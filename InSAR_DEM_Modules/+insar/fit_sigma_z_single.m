function out = fit_sigma_z_single(C, kz, opts)
% Fit sigma_z from coherence vs. vertical wavenumber (single-component volume model)
% Model: gamma = exp(-0.5 * (sigma_z * kz).^2)
%
% Inputs:
%   C    : coherence map (0..1)
%   kz   : vertical wavenumber map (rad/m)
%   opts : struct (all optional)
%       .mask           : logical mask of valid pixels
%       .gammaSNR       : map of gamma_snr to divide out (same size as C)
%       .minC           : floor coherence for log (default 0.05)
%       .robust         : true -> robustfit, false -> ordinary (default true)
%       .plot           : true -> quick diagnostic plot
%
% Output:
%   out.sigma_z        : fitted sigma_z (m)
%   out.slope,out.intercept : fit params for ln(Ccorr) = intercept + slope * kz.^2
%   out.r2, out.n, out.resid : diagnostics
%   out.gamma_model    : modeled gamma map
%   out.mask_used      : final sample mask

if nargin<3, opts=struct; end
if ~isfield(opts, 'mask'), opts.mask = true(size(C)); end
if ~isfield(opts,'gammaSNR'), opts.gammaSNR = []; end
if ~isfield(opts,'minC'), opts.minC = 0.5; end
if ~isfield(opts,'robust'), opts.robust = true; end
if ~isfield(opts,'doplot'), opts.doplot = false; end

% d = @(s,def) (isfield(opts,s) && ~isempty(opts.(s))) * opts.(s) + (~(isfield(opts,s)&&~isempty(opts.(s))))*def;
% mask   = d('mask', true(size(C)));
% gammaS = d('gammaSNR', []);
% minC   = d('minC', 0.05);
% useRob = d('robust', true);
% doplot = d('plot', false);

mask   = opts.mask;
gammaS = opts.gammaSNR;
minC   = opts.minC;
useRob = opts.robust;
doplot = opts.doplot;

% valid samples
M = isfinite(C) & isfinite(kz) & mask;
if ~isempty(gammaS), M = M & isfinite(gammaS) & gammaS>0; end
if ~any(M(:)), error('No valid samples for fitting.'); end

% SNR correction (optional)
Cc = C;
if ~isempty(gammaS)
    Cc = C ./ max(gammaS,1e-3);
    Cc = max(min(Cc,0.999), minC); % bound for log
end

x = (kz(M)).^2;
y = log(max(Cc(M), minC));

% Robust/ordinary linear fit: y = a + b*x
if useRob && exist('robustfit','file')
    [b, stats] = robustfit(x, y);
    a0 = b(1); b1 = b(2);
    resid = stats.resid;
else
    X = [ones(numel(x),1) x(:)];
    b = X \ y(:);
    a0 = b(1); b1 = b(2);
    yhat = X*b; resid = y(:) - yhat;
end

% From model: ln(gamma) = -0.5*sigma^2*kz^2  => slope = -0.5*sigma^2
sigma_z = sqrt(max(0, -2*b1));

% Diagnostics
sst = sum((y - mean(y)).^2);
ssr = sum(resid.^2);
R2  = max(0, 1 - ssr/max(sst,eps));

% Modeled gamma
gamma_model = exp(-0.5*(sigma_z.*kz).^2);
if ~isempty(gammaS), gamma_model = gamma_model .* max(gammaS,1e-3); end

out = struct('sigma_z',sigma_z,'slope',b1,'intercept',a0, ...
             'r2',R2,'n',numel(x),'resid',resid, ...
             'gamma_model',gamma_model,'mask_used',M);

if doplot
    figure(); binscatter(x, y); hold on; grid on
    xx = linspace(0, quantile(x,0.95), 200)'; yy = a0 + b1*xx;
    plot(xx, yy, 'r','LineWidth',2)
    xlabel('k_z^2'); ylabel('ln(\gamma_{corr})');
    title(sprintf('Single-comp fit: \\sigma_z = %.3f m, R^2=%.2f, n=%d', sigma_z, R2, out.n))
end
end

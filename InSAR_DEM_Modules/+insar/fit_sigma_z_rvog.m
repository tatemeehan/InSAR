function out = fit_sigma_z_rvog(C, kz, opts)
% Magnitude-only RVoG fit:
%   gamma_model = | alpha * exp(-0.5*(sigma_z*kz).^2) + (1-alpha) * exp(-1i * kz * z0) |
%
% Inputs:
%   C, kz : maps (same size)
%   opts (all optional):
%     .mask        : logical mask
%     .gammaSNR    : map to divide out
%     .bounds      : struct with fields:
%                       sigma_z [min max] (m), alpha [0 1], z0 [min max] (m)
%                    default: sigma_z [0 2], alpha [0 1], z0 [-1 1]
%     .init        : struct with fields: sigma_z, alpha, z0 (initial guesses)
%     .robust      : true to use robust loss in objective (Huber-like), default true
%     .plot        : diagnostics
%
% Output:
%   out.sigma_z, out.alpha, out.z0, out.gamma_model, out.rss, out.n, out.exitflag

if nargin<3, opts=struct; end
d = @(s,def) (isfield(opts,s) && ~isempty(opts.(s))) * opts.(s) + (~(isfield(opts,s)&&~isempty(opts.(s))))*def;

mask   = d('mask', true(size(C)));
gammaS = d('gammaSNR', []);
doplot = d('plot', false);
useRob = d('robust', true);

B = struct('sigma_z',[0 2],'alpha',[0 1],'z0',[-1 1]);  % meters for z0
if isfield(opts,'bounds'), B = utils_merge(B, opts.bounds); end

% valid samples
M = isfinite(C) & isfinite(kz) & mask;
if ~isempty(gammaS), M = M & isfinite(gammaS) & gammaS>0; end
if ~any(M(:)), error('No valid samples for fitting.'); end

% SNR correction (optional)
Cc = C;
if ~isempty(gammaS), Cc = C ./ max(gammaS,1e-3); Cc = min(max(Cc,0), 0.999); end

k = kz(M);  g = Cc(M);
n = numel(g);

% Initial guess from single-component fit
sc = fit_sigma_z_single(C, kz, struct('mask',mask,'gammaSNR',gammaS,'robust',true));
x0 = [ min(max(sc.sigma_z, B.sigma_z(1)), B.sigma_z(2)),  ... % sigma_z
       0.7, ...                                              % alpha
       0.0 ];                                                % z0 (m)
if isfield(opts,'init')
    if isfield(opts.init,'sigma_z'), x0(1)=opts.init.sigma_z; end
    if isfield(opts.init,'alpha'),   x0(2)=opts.init.alpha;   end
    if isfield(opts.init,'z0'),      x0(3)=opts.init.z0;      end
end

lb = [B.sigma_z(1) B.alpha(1) B.z0(1)];
ub = [B.sigma_z(2) B.alpha(2) B.z0(2)];

% Objective (magnitude residuals)
    function [rss, r] = obj(x)
        sig = x(1); a = x(2); z0 = x(3);
        vol = exp(-0.5*(sig.*k).^2);
        grd = exp(-1i*k*z0);
        gm  = abs(a*vol + (1-a)*grd);            % model magnitude
        r   = gm - g;                             % residuals
        if useRob
            % Huber-like soft robustification
            c  = 0.05;
            r  = sqrt(r.^2 + c^2) - c;
        end
        rss = sum(r.^2);
    end

% Optimizer: prefer fmincon if available, else fminsearch with boxing
useFmincon = exist('fmincon','file')==2;
if useFmincon
    optsf = optimoptions('fmincon','Display','none','Algorithm','interior-point');
    [x, fval, exitflag] = fmincon(@(x)obj(x), x0, [],[],[],[], lb, ub, [], optsf);
else
    % Simple bounded transform for fminsearch
    tf  = @(x) [logit((x(1)-lb(1))/(ub(1)-lb(1))), logit((x(2)-lb(2))/(ub(2)-lb(2))), ...
                logit((x(3)-lb(3))/(ub(3)-lb(3)))];
    itf = @(y) [lb(1)+(ub(1)-lb(1))*logistic(y(1)), lb(2)+(ub(2)-lb(2))*logistic(y(2)), ...
                lb(3)+(ub(3)-lb(3))*logistic(y(3))];
    y0  = tf(x0);
    optsnm = optimset('Display','off');
    [y, fval, exitflag] = fminsearch(@(y)obj(itf(y)), y0, optsnm);
    x = itf(y);
end

sig = x(1); a = x(2); z0 = x(3);
vol = exp(-0.5*(sig.*kz).^2);
grd = exp(-1i*kz*z0);
gamma_model = abs(a*vol + (1-a)*grd);
if ~isempty(gammaS), gamma_model = gamma_model .* max(gammaS,1e-3); end

out = struct('sigma_z',sig,'alpha',a,'z0',z0, ...
             'rss',fval,'n',n,'exitflag',exitflag, ...
             'gamma_model',gamma_model,'mask_used',M);

if doplot
    figure; tiledlayout(1,2,'Padding','compact','TileSpacing','compact')
    nexttile
    scatter(k.^2, g, 3,'.'); hold on; grid on
    km = linspace(0, max(k), 300)'; gm = abs(a*exp(-0.5*(sig.*km).^2) + (1-a)*exp(-1i*km*z0));
    plot(km.^2, gm, 'r','LineWidth',2)
    xlabel('k_z^2'); ylabel('|\gamma|'); title(sprintf('RVoG fit: \\sigma_z=%.3f m, \\alpha=%.2f, z_0=%.2f m', sig,a,z0))
    nexttile
    imagesc(gamma_model); axis image xy; colorbar
    title('Modeled |\gamma| (RVoG)')
end
end

function s = utils_merge(a,b)
s = a; f = fieldnames(b);
for i=1:numel(f), s.(f{i}) = b.(f{i}); end
end
function y = logit(x),    y = log(x./(1-x)); end
function x = logistic(y), x = 1./(1+exp(-y)); end

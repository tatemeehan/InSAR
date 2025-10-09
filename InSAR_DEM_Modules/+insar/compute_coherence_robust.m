function [coh, ccoh, meta] = compute_coherence_robust( ...
    slc1, slc2, kernelType, winSize, sigma, opts)

% Defaults / options
if nargin < 3 || isempty(kernelType), kernelType = 'gaussian'; end
if nargin < 4 || isempty(winSize),    winSize    = 7; end
if nargin < 5 || isempty(sigma),      sigma      = winSize/6; end
if nargin < 6, opts = struct(); end
if ~isfield(opts,'alpha'),  opts.alpha  = 2;      end     % base coherence -> weights
if ~isfield(opts,'iters'),  opts.iters  = 2;      end     % robust iters (1–3)
if ~isfield(opts,'loss'),   opts.loss   = 'huber';end     % 'huber'|'cauchy'|'tukey'
if ~isfield(opts,'tau'),    opts.tau    = 0.7;    end     % rad (≈40°)

% Ensure odd window
if mod(winSize,2)==0, winSize = winSize+1; end

% Kernel
switch lower(kernelType)
  case 'gaussian', G = fspecial('gaussian', winSize, sigma);
  case 'boxcar',   G = ones(winSize)/winSize^2;
  otherwise, error('Unknown kernelType');
end

% Ingredients
I  = slc1 .* conj(slc2);     % interferogram
P1 = abs(slc1).^2;           % power 1
P2 = abs(slc2).^2;           % power 2
U  = I ./ max(abs(I),eps);   % unit phasors e^{iφ}

% NaN-safe mask
M = isfinite(I) & isfinite(P1) & isfinite(P2);
I(~M)=0; P1(~M)=0; P2(~M)=0;

% ---- Pass 0: uniform ML coherence (for base weights) ----
num0 = imfilter(I,  G, 'symmetric');
d10  = imfilter(P1, G, 'symmetric');
d20  = imfilter(P2, G, 'symmetric');
cgamma0 = num0 ./ sqrt(d10.*d20 + eps);
coh0    = abs(cgamma0);

% Base spatial/adaptive weights
W = max(coh0,0).^opts.alpha;    % nonnegative
W(~M) = 0;

% ---- Robust reweighting on phase residuals (IRLS) ----
phi = angle(cgamma0);           % init with ML phase
for it = 1:max(0,opts.iters)
    % residual angles in (-pi,pi]
    e = angle( U .* exp(-1i*phi) );
    a = abs(e); t = opts.tau;

    switch lower(opts.loss)
      case 'huber'   % w = min(1, t/|e|)
        wrob = ones(size(a));
        m = a > t; wrob(m) = t ./ (a(m)+eps);

      case 'cauchy'  % w = 1 / (1 + (e/t)^2)
        wrob = 1 ./ (1 + (a./t).^2);

      case 'tukey'   % w = (1-(e/t)^2)^2 for |e|<t, else 0
        wrob = zeros(size(a));
        m = a < t;
        r = (a(m)./t);
        wrob(m) = (1 - r.^2).^2;

      otherwise
        error('Unknown opts.loss');
    end

    % combine with base weights
    Wtot = W .* wrob; Wtot(~M) = 0;

    % update complex circular mean phase (wrap-safe)
    numP = imfilter(Wtot .* U, G, 'symmetric');
    denP = imfilter(Wtot,     G, 'symmetric') + eps;
    Cphase = numP ./ denP;
    phi = angle(Cphase);

    % (optional) early stop if tiny change:
    % if it>1 && median(abs(angle(exp(1i*(phi - phi_prev)))), 'all') < 1e-3, break; end
    % phi_prev = phi;
    W = Wtot;  % carry robust weights forward
end

% ---- Final ML complex coherence with robust weights ----
num = imfilter(W .* I,  G, 'symmetric');
d1  = imfilter(W .* P1, G, 'symmetric');
d2  = imfilter(W .* P2, G, 'symmetric');

ccoh = num ./ sqrt(d1.*d2 + eps);    % complex ML (≤1)
coh  = abs(ccoh);

% Meta
meta.kernelType = kernelType;
meta.winSize    = winSize;
meta.sigma      = sigma;
meta.kernel     = G;
meta.alpha      = opts.alpha;
meta.robust     = rmfield(opts, setdiff(fieldnames(opts), {'alpha','iters','loss','tau'}));
meta.phase      = phi;          % robust phase used for residuals
meta.weights    = W;            % final weights used in ML

end

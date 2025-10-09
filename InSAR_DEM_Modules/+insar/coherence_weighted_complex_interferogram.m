function [Ifilt, phz, q] = coherence_weighted_complex_interferogram( ...
        ccor, winSize, sigma, alpha)
% Coherence-weighted complex filtering that mirrors your scalar
% coherence_weighted_vector_filter(), but operates on the complex field.
%
% INPUTS
%   ccor    : complex coherence (a.k.a. complex cross-correlation), same size as scene
%             angle(ccor) = filtered phase estimate, abs(ccor) = coherence gamma in [0,1]
%   winSize : odd integer, e.g., 5 or 7
%   sigma   : Gaussian std-dev (pixels)
%   alpha   : coherence exponent (>=0), e.g., 2
%
% OUTPUTS
%   Ifilt : complex, coherence-weighted & spatially-filtered interferogram (phasor-like)
%   phz   : angle(Ifilt)  (filtered wrapped phase)
%   q     : abs(Ifilt)    (vector strength; tracks phase consistency)

    % --- input checks mirroring your scalar filter ---
    assert(isnumeric(ccor) && ~isreal(ccor), 'ccor must be complex');
    if nargin < 4 || isempty(alpha),   alpha = 2; end
    if nargin < 3 || isempty(sigma),   sigma = winSize/2; end
    if nargin < 2 || isempty(winSize), winSize = 5; end
    assert(alpha >= 0, 'alpha must be non-negative');
    assert(mod(winSize,2)==1, 'winSize must be an odd integer');

    % --- kernel: SAME normalized Gaussian as your code ---
    halfWin = floor(winSize/2);
    [xg, yg] = meshgrid(-halfWin:halfWin, -halfWin:halfWin);
    G = exp(-(xg.^2 + yg.^2) / (2 * sigma^2));
    G = G / sum(G(:));

    % --- unit phasors from complex coherence (avoid amplitude bias) ---
    phi = angle(ccor);
    Z   = exp(1i * phi);        % unit phasor carrying only phase
    gam = abs(ccor);            % coherence magnitude in [0,1]
    gam = min(max(gam,0),1);    % clip

    % NaN handling (exclude invalids from contributions)
    msk = isfinite(phi) & isfinite(gam);
    Z(~msk)   = 0;
    gam(~msk) = 0;

    % --- weights: gamma^alpha, identical to your scalar filter ---
    W = gam.^alpha;

    % --- vectorized equivalent of per-pixel weighted vector sum ---
    %      numerator = sum( G .* (Z .* W) )    (complex)
    %      denom     = sum( G .*  W      )     (real)
    num = imfilter(Z .* W, G, 'replicate');   % complex
    den = imfilter(W,     G, 'replicate');    % real
    denSafe = den; denSafe(denSafe==0) = eps;

    Ifilt = num ./ denSafe;                   % complex filtered interferogram

    % Optional exact fallback where window had zero total weight (matches your loop)
    zeroWin = (den == 0);
    if any(zeroWin,'all')
        Ifilt(zeroWin) = exp(1i * phi(zeroWin));
    end

    phz = angle(Ifilt);
    q   = abs(Ifilt);                          % vector strength / quality proxy
end

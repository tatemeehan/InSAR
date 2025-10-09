function [Cphase, phi_phaseonly, meta] = insar_phase_complex(slc1, slc2, winSize, sigma)
% Phase-only complex field via circular mean of unit phasors.
% No amplitude weighting; optional lookMask (logical) prevents NaN bleed.

if nargin < 3 || isempty(winSize), winSize = 5; end
if nargin < 4 || isempty(sigma),   sigma   = winSize/2; end
if mod(winSize,2)==0, winSize = winSize+1; end

G = fspecial('gaussian', winSize, sigma);

I  = slc1 .* conj(slc2);
U  = I ./ max(abs(I), eps);                         % unit phasors e^{iφ}

% validity mask (and optional lookMask)
M = isfinite(real(U)) & isfinite(imag(U));
if nargin >= 5 && ~isempty(lookMask), M = M & lookMask; end

% zero-out invalids; normalize by local support to avoid border shrinkage
U(~M) = 0;
num   = imfilter(U,            G, 'symmetric');      % complex
den   = imfilter(double(M),    G, 'symmetric') + eps; % real support count
Cphase = num ./ den;                                  % complex phase-only
phi_phaseonly = angle(Cphase);

meta.kernel = G; meta.support = den; meta.mask = M;
end

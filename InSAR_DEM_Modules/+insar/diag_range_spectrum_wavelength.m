function [lambda_px, Pmean] = diag_range_spectrum_wavelength(Ifilt, coh, pxSize, minSamples)
% pxSize is range pixel spacing in meters (use 1 if 1 m/px). Only used to label in meters if you want.
if nargin<3 || isempty(pxSize), pxSize = 1; end
if nargin<4, minSamples = 100; end

[M,N] = size(Ifilt);
phW   = angle(Ifilt);
mask  = isfinite(coh) & coh>=0.35 & isfinite(phW);

% Work out a common FFT length from the longest valid column
Lmax = 0;
Lcols = zeros(1,N);
for j=1:N
    v = mask(:,j);
    Lcols(j) = nnz(v);
    Lmax = max(Lmax, Lcols(j));
end
if Lmax==0, lambda_px = []; Pmean = []; return; end
Lfft = 2^nextpow2(Lmax);             % common length
K    = Lfft/2 + 1;                   % one-sided bins
Pacc = zeros(K,1);                   % accumulator
nacc = 0;

for j=1:N
    v = mask(:,j);
    L = nnz(v);
    if L < minSamples, continue; end

    ph = unwrap(phW(v,j));
    w  = coh(v,j);
    ph = ph - median(ph,'omitnan');
    % taper and weight
    win = hann(L);
    z   = (ph .* win) .* (w / max(median(w),eps));

    % zero-pad to Lfft and one-sided spectrum
    Z = fft(z, Lfft);
    P = (abs(Z(1:K)).^2) / L;        % power ~ rad^2
    Pacc = Pacc + P;
    nacc = nacc + 1;
end

if nacc==0, lambda_px = []; Pmean=[]; return; end
Pmean = Pacc / nacc;

% frequency axis (cyc/px) for one-sided spectrum
f = (0:K-1).' / Lfft;       % 0 .. 0.5 cyc/px

% convert to wavelength (px); drop DC bin
f = f(2:end);
Pmean = Pmean(2:end);
lambda_px = 1 ./ f;         % pixels per cycle

% Example plots:
figure(); 
subplot(1,2,1);
plot(f, Pmean); xlabel('Spatial frequency (cycles/pixel)'); ylabel('Power (rad^2)');
xlim([0 0.5]);

subplot(1,2,2);
plot(lambda_px*pxSize, Pmean); set(gca,'XDir','reverse');  % long λ on the left
xlabel('Wavelength (m)'); ylabel('Power (rad^2)');
end

function [phaseOut, zOut, Rlocal, nValid, fracValid] = ...
    nan_component_median_filter_phase(phase, window_size, varargin)
%NAN_COMPONENT_MEDIAN_FILTER_PHASE NaN-aware component-wise median phase filter.
%
% Converts wrapped phase to phasors, then applies local nanmedian to real
% and imaginary parts separately.
%
% This is not a true vector median, but is much cheaper and wrap-safe.

p = inputParser;
p.addParameter('MinValidFrac', 0.10, @isnumeric);
p.addParameter('PreserveInputMask', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

minValidFrac = p.Results.MinValidFrac;
preserveInputMask = logical(p.Results.PreserveInputMask);

if mod(window_size, 2) == 0
    error('window_size must be odd.');
end

if isreal(phase)
    valid0 = isfinite(phase);
    z = complex(nan(size(phase)), nan(size(phase)));
    z(valid0) = exp(1i .* double(phase(valid0)));
else
    z = complex(double(real(phase)), double(imag(phase)));
    valid0 = isfinite(real(z)) & isfinite(imag(z));
end

realZ = real(z);
imagZ = imag(z);

realMed = nlfilter(realZ, [window_size window_size], @(x) median(x(isfinite(x)), 'omitnan'));
imagMed = nlfilter(imagZ, [window_size window_size], @(x) median(x(isfinite(x)), 'omitnan'));

validCountFcn = @(x) sum(isfinite(x(:)));
nValid = nlfilter(realZ, [window_size window_size], validCountFcn);
fracValid = nValid ./ (window_size * window_size);

zOut = complex(realMed, imagMed);

bad = nValid == 0 | fracValid < minValidFrac;
zOut(bad) = NaN;

Rlocal = abs(zOut);
phaseOut = angle(zOut);

if preserveInputMask
    phaseOut(~valid0) = NaN;
    zOut(~valid0) = NaN;
    Rlocal(~valid0) = NaN;
end
end
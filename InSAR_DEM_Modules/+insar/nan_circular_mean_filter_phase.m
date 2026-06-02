function [phaseOut, zOut, Rlocal, nValid, fracValid] = ...
    nan_circular_mean_filter_phase(phase, window_size, varargin)
%NAN_CIRCULAR_MEAN_FILTER_PHASE Fast NaN-aware circular phase filter.
%
% This filters wrapped phase by averaging unit phasors:
%   z = exp(1i*phase)
% while ignoring NaNs.
%
% Outputs:
%   phaseOut  : filtered wrapped phase [rad]
%   zOut      : filtered complex phasor
%   Rlocal    : local phase concentration, 0 to 1
%   nValid    : number of finite samples in window
%   fracValid : finite fraction in window

p = inputParser;
p.addParameter('MinValidFrac', 0.10, @isnumeric);
p.addParameter('PreserveInputMask', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

minValidFrac = p.Results.MinValidFrac;
preserveInputMask = logical(p.Results.PreserveInputMask);

if mod(window_size, 2) == 0
    error('window_size must be odd.');
end

phase = double(phase);
valid0 = isfinite(phase);

z = complex(zeros(size(phase)), zeros(size(phase)));
z(valid0) = exp(1i .* phase(valid0));

ker = ones(window_size, window_size);

sumZ = conv2(z, ker, 'same');
nValid = conv2(double(valid0), ker, 'same');
fracValid = nValid ./ numel(ker);

zOut = complex(nan(size(phase)), nan(size(phase)));

good = nValid > 0 & fracValid >= minValidFrac;
zOut(good) = sumZ(good) ./ nValid(good);

Rlocal = abs(zOut);
phaseOut = angle(zOut);

if preserveInputMask
    phaseOut(~valid0) = NaN;
    zOut(~valid0) = NaN;
    Rlocal(~valid0) = NaN;
end
end
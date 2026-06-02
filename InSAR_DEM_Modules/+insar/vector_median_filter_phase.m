function [phaseOut, zOut, R] = vector_median_filter_phase(phase, window_size)
%VECTOR_MEDIAN_FILTER_PHASE Component-wise complex median filter for phase.
%
% Inputs
%   phase       : wrapped phase image [rad], or complex phasor/image
%   window_size : odd integer window size
%
% Outputs
%   phaseOut : filtered wrapped phase [rad]
%   zOut     : filtered complex phasor/image
%   R        : local filtered phasor magnitude, useful as reliability

if mod(window_size, 2) == 0
    error('window_size must be odd.');
end

if isreal(phase)
    z = exp(1i .* phase);
else
    z = phase;
end

real_med = medfilt2(real(z), [window_size window_size], 'symmetric');
imag_med = medfilt2(imag(z), [window_size window_size], 'symmetric');

zOut = complex(real_med, imag_med);

R = abs(zOut);

phaseOut = angle(zOut);

end
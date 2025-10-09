function out = vector_median_filter(phase, window_size)
if isreal(phase)
    z = exp(1i * phase);
    real_med = medfilt2(real(z), [window_size window_size], 'symmetric');
    imag_med = medfilt2(imag(z), [window_size window_size], 'symmetric');
    out = atan2(imag_med, real_med);
else
    real_med = medfilt2(real(phase), [window_size window_size], 'symmetric');
    imag_med = medfilt2(imag(phase), [window_size window_size], 'symmetric');
    out = real_med+1i.*imag_med;
end

end
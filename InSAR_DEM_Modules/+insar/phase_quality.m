function quality = phase_quality(wrapped_phase, window_size)
    z = exp(1i * wrapped_phase);
    mask = ~isnan(z);
    kernel = fspecial('average', window_size);

    % Normalize kernel to valid entries only
    norm_factor = imfilter(double(mask), kernel, 'replicate');
    norm_factor(norm_factor == 0) = NaN;

    % Compute filtered real/imag parts
    real_z = real(z); real_z(~mask) = 0;
    imag_z = imag(z); imag_z(~mask) = 0;
    mean_real = imfilter(real_z, kernel, 'replicate') ./ norm_factor;
    mean_imag = imfilter(imag_z, kernel, 'replicate') ./ norm_factor;

    % Compute complex coherence ratio
    local_mag = imfilter(abs(z), kernel, 'replicate') ./ norm_factor;
    local_mag(local_mag == 0) = NaN;  % avoid divide-by-zero
    quality = sqrt(mean_real.^2 + mean_imag.^2) ./ local_mag;
end
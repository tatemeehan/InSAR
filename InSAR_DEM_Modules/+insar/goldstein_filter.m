function out = goldstein_filter(phase, alpha)
    F = fft2(exp(1i * phase));
    S = fftshift(F);
    [rows, cols] = size(S);
    [x, y] = meshgrid(-floor(cols/2):floor((cols-1)/2), -floor(rows/2):floor((rows-1)/2));
    R = sqrt(x.^2 + y.^2);
    R = R / max(R(:));
    H = R .^ alpha;
    H(isnan(H)) = 0;
    S_filtered = S .* H;
    F_filtered = ifftshift(S_filtered);
    out = angle(ifft2(F_filtered));
end
function gamma_snr = coherence_ceiling_from_snr(SNR1, SNR2)
valid = isfinite(SNR1) & isfinite(SNR2) & SNR1>0 & SNR2>0;
gamma_snr = nan(size(SNR1));
gamma_snr(valid) = 1 ./ sqrt( (1 + 1./SNR1(valid)) .* (1 + 1./SNR2(valid)) );
end

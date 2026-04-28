function eps_eff = eps_eff_energy(diag)
% Computes the Energy Weighted Integral of the 1-D (z) Complex Permittivity
    z = diag.z;
    W = diag.W_energy;
    eps_z = diag.eps;

    den = trapz(z, W);
    if den <= 0 || ~isfinite(den)
        eps_eff = NaN + 1i*NaN;
        return;
    end

    eps_eff = trapz(z, eps_z .* W) / den;
end

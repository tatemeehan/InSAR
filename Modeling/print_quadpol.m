function print_quadpol(out)
%PRINT_QUADPOL Pretty-print 2x2 quad-pol matrices for the benchmark model.
% Works whether out.Ag/Eg/E/Es are scalar or 2x2. Assumes:
%   [HH HV; VH VV] ordering.

    fprintf('\n==== QUAD-POL SUMMARY ====\n');

    print_mat('Ag', out.Ag);
    print_mat('Es', out.Es);
    print_mat('Eg', out.Eg);
    print_mat('E ', out.E);

    fprintf('--------------------------\n');
    if isfield(out,'profile') && isfield(out.profile,'enable') && out.profile.enable
        if isfield(out.profile,'diag') && isfield(out.profile.diag,'sigma_eff')
            fprintf('sigma_eff (energy-weighted) = %.4g S/m\n', out.profile.diag.sigma_eff);
        end
        if isfield(out,'Dg')
            fprintf('Dg (effective) = %.4g%+.4gi\n', real(out.Dg), imag(out.Dg));
        end
    end
    fprintf('==========================\n\n');
end

function print_mat(name, X)
    X = ensure_2x2(X);

    HH = X(1,1); HV = X(1,2);
    VH = X(2,1); VV = X(2,2);

    fprintf('%s HH: |%.3g|, ang=%7.2f deg\n', name, abs(HH), rad2deg(angle(HH)));
    fprintf('%s VV: |%.3g|, ang=%7.2f deg\n', name, abs(VV), rad2deg(angle(VV)));
    fprintf('%s HV: |%.3g|, ang=%7.2f deg\n', name, abs(HV), rad2deg(angle(HV)));
    fprintf('%s VH: |%.3g|, ang=%7.2f deg\n', name, abs(VH), rad2deg(angle(VH)));
end

function X2 = ensure_2x2(X)
    if isempty(X)
        X2 = complex(zeros(2,2));
    elseif isscalar(X)
        X2 = [X 0; 0 X];
    elseif isequal(size(X),[2 2])
        X2 = X;
    else
        error('Expected scalar or 2x2, got %s', mat2str(size(X)));
    end
end

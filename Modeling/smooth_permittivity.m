function eps_smooth = smooth_permittivity(freq, eps_raw, method, param, makePlot)
%SMOOTH_PERMITTIVITY Smooth complex permittivity using flexible methods
%
% eps_smooth = smooth_permittivity(freq, eps_raw, method, param, makePlot)
%
% Inputs:
%   freq      : frequency vector [Hz]
%   eps_raw   : complex permittivity (vector)
%   method    : 'sgolay' | 'spline' | 'lowess' | 'none'
%   param     : struct of method-specific parameters
%   makePlot  : true/false
%
% Output:
%   eps_smooth : smoothed complex permittivity (same size)

% Ensure column
freq = freq(:);
eps_raw = eps_raw(:);

switch lower(method)
    case 'sgolay'
        order    = getOpt(param, 'order', 3);
        framelen = getOpt(param, 'framelen', 9);
        real_s = sgolayfilt(real(eps_raw), order, framelen);
        imag_s = sgolayfilt(imag(eps_raw), order, framelen);
        eps_smooth = real_s + 1i * imag_s;

    case 'spline'
        p = getOpt(param, 'smoothness', 0.99);
        sp_r = csaps(freq, real(eps_raw), p);
        sp_i = csaps(freq, imag(eps_raw), p);
        real_s = fnval(sp_r, freq);
        imag_s = fnval(sp_i, freq);
        eps_smooth = real_s + 1i * imag_s;

    case 'lowess'
        span = getOpt(param, 'span', 0.05);
        real_s = smooth(freq, real(eps_raw), span, 'lowess');
        imag_s = smooth(freq, imag(eps_raw), span, 'lowess');
        eps_smooth = real_s + 1i * imag_s;

    case 'none'
        eps_smooth = eps_raw;

    otherwise
        error('Unknown smoothing method: %s', method);
end

% Optional plot
if makePlot
    figure; 
    subplot(2,1,1); hold on;
    plot(freq, real(eps_raw), 'ko-', 'DisplayName', 'Raw');
    plot(freq, real(eps_smooth), 'r-', 'DisplayName', 'Smoothed');
    ylabel('\epsilon'''); legend; grid on

    subplot(2,1,2); hold on;
    plot(freq, imag(eps_raw), 'ko-', 'DisplayName', 'Raw');
    plot(freq, imag(eps_smooth), 'b-', 'DisplayName', 'Smoothed');
    ylabel('\epsilon'''''); xlabel('Frequency [Hz]'); legend; grid on

    sgtitle(sprintf('Permittivity Smoothing (%s)', method));
end

end

function val = getOpt(s, f, default)
if isfield(s, f)
    val = s.(f);
else
    val = default;
end
end

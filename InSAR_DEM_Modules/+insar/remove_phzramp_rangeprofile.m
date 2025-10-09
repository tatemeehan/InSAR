function [I_corr, prof, validR] = remove_phzramp_rangeprofile(I_filt, coh, orient, sigmaR, cohThresh)
% orient: 'range' (default) or 'azimuth'
% sigmaR: Gaussian sigma (px) to smooth the 1-D profile (e.g., 40–80)
% cohThresh: mask threshold (default 0.3)

if nargin<3 || isempty(orient),    orient   = 'range'; end
if nargin<4 || isempty(sigmaR),    sigmaR   = 60;      end
if nargin<5 || isempty(cohThresh), cohThresh= 0.30;    end

[M,N] = size(I_filt);
valid = isfinite(I_filt) & isfinite(coh) & (coh>=cohThresh);
U = I_filt ./ max(abs(I_filt), eps);

switch lower(orient)
  case 'range'   % profile vs columns (range bins)
    prof = NaN(1,N); validR = false(1,N);
    for j = 1:N
        v = valid(:,j);
        if any(v)
            w = coh(v,j);
            z = U(v,j);
            zbar = sum(w .* z) / max(sum(w), eps);   % coherence-weighted circular mean
            prof(j) = angle(zbar);
            validR(j) = true;
        end
    end
    prof = smooth1Dnan(prof, sigmaR);
    phi = repmat(prof, M, 1);
  case 'azimuth' % profile vs rows
    prof = NaN(M,1); validR = false(M,1);
    for i = 1:M
        v = valid(i,:)';
        if any(v)
            w = coh(i,v);
            z = U(i,v);
            zbar = sum(w .* z) / max(sum(w), eps);
            prof(i) = angle(zbar);
            validR(i) = true;
        end
    end
    prof = smooth1Dnan(prof, sigmaR);
    phi = repmat(prof, 1, N);
  otherwise
    error('orient must be ''range'' or ''azimuth''');
end

% remove in complex domain
I_corr = I_filt .* exp(-1i * phi);

end

function y = smooth1Dnan(x, sigma)
    % NaN-aware 1D Gaussian smoother (wrap-safe profile)
    if sigma<=0, y=x; return; end
    L = max(3, ceil(5*sigma)); g = exp((-( -L:L ).^2)/(2*sigma^2)); g = g/sum(g);
    m = ~isnan(x); x2 = x; x2(~m) = 0; m2 = double(m);
    num = conv(x2, g, 'same'); den = conv(m2, g, 'same'); 
    y = num ./ max(den, eps); y(den==0) = NaN;
    % center profile to zero mean so gauge-free
    y = y - median(y(isfinite(y)), 'omitnan');
end

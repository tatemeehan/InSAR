function [sigma_hat, sigma_ci] = sigma1_bootstrap(kz, gamma_meas, gamma_snr, nBoot,obj)
% model: gamma_model = exp(-0.5 * kz2 * sigma^2) .* gamma_snr
kz2 = kz.^2;
valid = isfinite(kz2) & isfinite(gamma_meas) & isfinite(gamma_snr);
x = kz2(valid); y = gamma_meas(valid); s = gamma_snr(valid);
if nargin<5
    obj = 'huber';
end
if strcmp(obj,'huber')
    isMedianAbs = 0;
    isHuber = 1;
elseif strcmp(obj,'medianAbs')
    isMedianAbs = 1;
    isHuber = 0;
else
    error('Uknown Objective Function')
end
if isMedianAbs
    obj = @(sig) median(abs(y - exp(-0.5*x*(sig.^2)).*s),'omitnan');  % Huber/median L1-ish
    sigma_hat = fminbnd(obj, 0, 1);  % adjust upper bound if needed
elseif isHuber
    % Weighted single-parameter fit (L2 with Huber-like weights)
    w = sqrt(max(s,1e-3));
    f = @(sig) (w .* (y - exp(-0.5*x*(sig.^2)).*s));
    sigma_hat = fminbnd(@(sig) huberNorm(f(sig), 0.05), 0, 1);
end

if nargin<4, nBoot = 250; end
N = numel(x);
boot = zeros(nBoot,1);
for b=1:nBoot
    % idx = randi(N,N,1);
    idx = datasample(1:N,round(0.05.*N),Replace=false); idx = idx(:);
    if isMedianAbs
        objb = @(sig) median(abs(y(idx) - exp(-0.5.*x(idx).*(sig.^2)).*s(idx)),'omitnan');
        boot(b) = fminbnd(objb, 0, 1);
    elseif isHuber
        f = @(sig) (w(idx) .* (y(idx) - exp(-0.5.*x(idx).*(sig.^2)).*s(idx)));
        boot(b) = fminbnd(@(sig) huberNorm(f(sig), 0.05), 0, 1);
    end
    
end
sigma_ci = prctile(boot,[16 84]);  % ~68% interval
end


function v = huberNorm(r,delta), a = abs(r); v = sum((a<=delta).*(0.5*r.^2) + (a>delta).*(delta*(a-delta/2)),'omitnan'); end

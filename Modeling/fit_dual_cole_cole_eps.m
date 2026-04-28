% function [params_fit, eps_model, residual, exitflag] = fit_dual_cole_cole_eps(freq, eps_meas)
function [params_fit, eps_model, residual] = fit_dual_cole_cole_eps(freq, eps_meas)

% FIT_DUAL_COLE_COLE_EPS fits the dual Cole-Cole model to measured permittivity.
% freq      : frequency vector (Hz)
% eps_meas  : measured complex permittivity
% params_fit: fitted parameters
% eps_model : modeled permittivity at fitted parameters
% residual  : residual vector
% exitflag  : exit flag from lsqnonlin

eps0 = 8.854187817e-12;

% --- Initial parameter guess: [eps_inf, eps_s1, eps_s2, tau1, alpha1, tau2, alpha2, sigma]
p0 = [5,     20,     10,     1e-10,   0.1,    1e-11,   0.1,    0.01];

% --- Parameter bounds
lb = [1,     5,      1,      5e-13,   0.0,    5e-13,   0.0,    0.0];
ub = [10,    40,     30,     5e-8,    1.0,    5e-8,    1.0,    2.0];

% % --- Cost function
% cost_fun = @(p) cole_cole_error(freq, eps_meas, p);
% 
% % --- Optimization options
% opts = optimoptions('lsqnonlin', 'Display', 'iter', ...
%     'MaxFunctionEvaluations', 5000, ...
%     'FunctionTolerance', 1e-9);
% 
% % --- Run optimization
% [params_fit, ~, residual, exitflag] = lsqnonlin(cost_fun, p0, lb, ub, opts);
% 

optsCost = struct( ...
    'useRelativeError', false, ...
    'useRegularization', true, ...
    'lambda', 0.01, ...
    'plotModel', false);  % You can toggle on for inspection

costFcn = @(p) compute_dual_cole_cole_cost(p, freq, eps_meas, optsCost);

% Run optimization
[params_fit, residual] = fmincon(costFcn, p0, [], [], [], [], lb, ub, []);

% % --- Evaluate fitted model
eps_model = compute_dual_cole_cole(freq, params_fit);


end

function err = cole_cole_error(freq, eps_meas, p)

eps_model = compute_dual_cole_cole(freq, p);

% Compute relative residuals (optional improvement)
% rel_err_real = (real(eps_model) - real(eps_meas)) ./ max(abs(real(eps_meas)), 1e-3);
% rel_err_imag = (imag(eps_model) - imag(eps_meas)) ./ max(abs(imag(eps_meas)), 1e-3);

rel_err_real = (real(eps_model) - real(eps_meas)).*0.375; %./ max(abs(real(eps_meas)), 1e-3);
rel_err_imag = (imag(eps_model) - imag(eps_meas)).*0.625; %./ max(abs(imag(eps_meas)), 1e-3);

% Stack errors for real + imaginary parts
err = [rel_err_real; rel_err_imag];
end

function eps_model = compute_dual_cole_cole(freq, p)
% Unpack parameters
eps_inf  = p(1);
eps_s1   = p(2);
eps_s2   = p(3);
tau1     = p(4);
alpha1   = p(5);
tau2     = p(6);
alpha2   = p(7);
sigma    = p(8);

eps0 = 8.854187817e-12;
omega = 2 * pi * freq(:);

% Relaxation terms
delta_eps1 = eps_s1 - eps_inf;
delta_eps2 = eps_s2 - eps_s1;

term1 = delta_eps1 ./ (1 + (1i * omega * tau1).^(1 - alpha1));
term2 = delta_eps2 ./ (1 + (1i * omega * tau2).^(1 - alpha2));
cond  = 1i * sigma ./ (omega * eps0);

% Final permittivity model
eps_model = eps_inf + term1 + term2 + cond;
end

function cost = compute_dual_cole_cole_cost(params, freq, eps_meas, opts)
% compute_dual_cole_cole_cost: Evaluate cost between dual Cole–Cole model and measured permittivity
%
% Inputs:
%   - params: [eps_inf, eps_s1, eps_s2, tau1, alpha1, tau2, alpha2, sigma]
%   - freq:   frequency vector (Hz)
%   - eps_meas: complex measured permittivity (Nx1)
%   - opts: struct with fields:
%       * useRelativeError (default: false)
%       * useRegularization (default: false)
%       * lambda (default: 0.01)
%       * plotModel (default: false)
%
% Output:
%   - cost: scalar cost value

% --- 0) Unpack parameters ---
eps_inf = params(1);
eps_s1  = params(2);
eps_s2  = params(3);
tau1    = params(4);
alpha1  = params(5);
tau2    = params(6);
alpha2  = params(7);
sigma   = params(8);

% --- 1) Set defaults ---
if ~isfield(opts, 'useRelativeError'); opts.useRelativeError = false; end
if ~isfield(opts, 'useRegularization'); opts.useRegularization = false; end
if ~isfield(opts, 'lambda'); opts.lambda = 0.01; end
if ~isfield(opts, 'plotModel'); opts.plotModel = false; end

% --- 2) Physical constants ---
eps0 = 8.854187817e-12;
omega = 2 * pi * freq;

% --- 3) Compute Cole-Cole model ---
delta_eps1 = eps_s1 - eps_inf;
delta_eps2 = eps_s2 - eps_s1;

term1 = delta_eps1 ./ (1 + (1i * omega * tau1).^(1 - alpha1));
term2 = delta_eps2 ./ (1 + (1i * omega * tau2).^(1 - alpha2));
cond  = 1i * sigma ./ (omega * eps0);

eps_model = eps_inf + term1 + term2 + cond;

% --- 4) Compute Error ---
real_err = real(eps_model) - real(eps_meas);
imag_err = imag(eps_model) - imag(eps_meas);

if opts.useRelativeError
    eps_real_safe = max(abs(real(eps_meas)), 1e-3);
    eps_imag_safe = max(abs(imag(eps_meas)), 1e-3);

    real_err = real_err ./ eps_real_safe;
    imag_err = imag_err ./ eps_imag_safe;
end

% --- 5) Combine cost ---
cost = sqrt(mean(real_err.^2 + imag_err.^2));

% --- 6) Optional regularization ---
if opts.useRegularization
    cost = cost + opts.lambda * (sigma^2 + tau1^2 + tau2^2);
end

% --- 7) Optional plot ---
if opts.plotModel
    figure(99); clf
    subplot(2,1,1); hold on;
    plot(freq, real(eps_meas), 'ko-');
    plot(freq, real(eps_model), 'r.-');
    ylabel('\epsilon'''); grid on;

    subplot(2,1,2); hold on;
    plot(freq, imag(eps_meas), 'ko-');
    plot(freq, imag(eps_model), 'b.-');
    ylabel('\epsilon'''''); xlabel('Freq [Hz]'); grid on;
    sgtitle(sprintf('Cost: %.4f', cost));
end

end

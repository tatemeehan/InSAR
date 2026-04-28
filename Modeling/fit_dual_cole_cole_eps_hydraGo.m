function [params_fit, eps_model, residual] = fit_dual_cole_cole_eps_hydraGo( ...
    freq, eps_meas, freq_hydra, eps_hydra, weight_hydra, optsCost, eps_std)

% FIT_DUAL_COLE_COLE_EPS fits a dual Cole–Cole model to complex permittivity data,
% including low-frequency anchor from HydraGO Flex.
%
% Inputs:
%   freq         : frequency vector for VNA measurements (Hz)
%   eps_meas     : complex permittivity at freq
%   freq_hydra   : single anchor frequency from HydraGO (e.g., 50e6 Hz)
%   eps_hydra    : complex permittivity from HydraGO
%   weight_hydra : scalar weighting factor (e.g., 5 or 10)
%   optsCost     : struct with fields:
%                   .useRelativeError (bool)
%                   .useRegularization (bool)
%                   .lambda (e.g., 0.01)
%                   .plotModel (true/false)
%                   .sigma_hydra (optional: trust HydraGO conductivity)

if isfield(optsCost, 'sigma_hydra')
    sigma_min = 0.5 * optsCost.sigma_hydra;
    sigma_max = 1.1 * optsCost.sigma_hydra;

    % --- Initial parameter guess: [eps_inf, eps_s1, eps_s2, tau1, alpha1, tau2, alpha2, sigma]
p0 = [5,     20,     10,     1e-10,   0.1,    1e-11,   0.1,   optsCost.sigma_hydra];

% --- Parameter bounds
lb = [1,     5,      1,      5e-13,   0.0,    5e-13,   0.0,    sigma_min];
ub = [20,    40,     30,     5e-8,    1.0,    5e-8,    1.0,    sigma_max];
else
% --- Initial parameter guess: [eps_inf, eps_s1, eps_s2, tau1, alpha1, tau2, alpha2, sigma]
p0 = [5,     20,     10,     1e-10,   0.1,    1e-11,   0.1,    0.01];

% --- Parameter bounds
lb = [1,     5,      1,      5e-13,   0.0,    5e-13,   0.0,    0.0];
ub = [20,    40,     30,     5e-8,    1.0,    5e-8,    1.0,    0.1];
end

% Frequencies to enforce relaxation positivity
freq_constraint = linspace(50e6, 4e9, 50);  % Can adjust resolution

% --- Cost function
costFcn = @(p) compute_dual_cole_cole_cost( ...
    p, freq, eps_meas, freq_hydra, eps_hydra, weight_hydra, optsCost,eps_std);

% Constraint function
nonlincon = @(p) dual_cole_cole_relaxation_constraint(p, freq_constraint);

% Fit with nonlinear constraint
params_fit = fmincon(costFcn, p0, [], [], [], [], lb, ub, nonlincon);

% --- Run optimization
% params_fit = fmincon(costFcn, p0, [], [], [], [], lb, ub, []);

% --- Evaluate final model
eps_model = compute_dual_cole_cole(freq, params_fit);
residual  = costFcn(params_fit);

end

function cost = compute_dual_cole_cole_cost(params, freq, eps_meas, ...
    freq_hydra, eps_hydra, weight_hydra, opts, eps_std)

% Unpack parameters
eps_inf = params(1); eps_s1 = params(2); eps_s2 = params(3);
tau1 = params(4);    alpha1 = params(5);
tau2 = params(6);    alpha2 = params(7);
sigma = params(8);

% Physical constants
eps0 = 8.854187817e-12;
omega = 2 * pi * freq;
omega_hydra = 2 * pi * freq_hydra;

% --- VNA Model ---
delta_eps1 = eps_s1 - eps_inf;
delta_eps2 = eps_s2 - eps_s1;

term1 = delta_eps1 ./ (1 + (1i * omega * tau1).^(1 - alpha1));
term2 = delta_eps2 ./ (1 + (1i * omega * tau2).^(1 - alpha2));
cond  = 1i * sigma ./ (omega * eps0);
eps_model = eps_inf + term1 + term2 + cond;

% --- Penalty for nonphysical negative imaginary permittivity
% imag_eps = imag(eps_model);
% neg_imag = imag_eps(imag_eps <= 0);
% penalty = 0;
% if ~isempty(neg_imag)
%     penalty = 1e2 * sum(neg_imag.^2);  % Soft penalty for ε'' ≤ 0
% end

% % Soft barrier penalty
% eps_floor = 0.02;
% barrier_strength = 2;
% imag_eps = imag(eps_model);
% penalty = barrier_strength * sum(exp(-imag_eps / eps_floor));

% --- VNA Error ---
real_err = (real(eps_model) - real(eps_meas))./real(eps_std);
imag_err = (imag(eps_model) - imag(eps_meas))./imag(eps_std);

if opts.useRelativeError
    real_err = real_err ./ max(abs(real(eps_meas)), 1e-3);
    imag_err = imag_err ./ max(abs(imag(eps_meas)), 1e-3);
end

vna_cost = sqrt(mean(real_err.^2 + imag_err.^2));

% --- HydraGO Error (50 MHz) ---
if weight_hydra > 0
    term1_h = delta_eps1 / (1 + (1i * omega_hydra * tau1)^(1 - alpha1));
    term2_h = delta_eps2 / (1 + (1i * omega_hydra * tau2)^(1 - alpha2));
    cond_h  = 1i * sigma / (omega_hydra * eps0);
    eps_model_h = eps_inf + term1_h + term2_h + cond_h;

    hydra_cost = sqrt((real(eps_model_h) - real(eps_hydra))^2 + ...
        (imag(eps_model_h) - imag(eps_hydra))^2);
else
    hydra_cost = 0;
end
% % --- Final Cost ---
% cost = vna_cost + weight_hydra * hydra_cost;
% 
% % --- Regularization (optional) ---
% if opts.useRegularization
%     if isfield(opts, 'sigma_hydra') && weight_hydra > 0
%         cost = cost + opts.lambda * ((sigma - opts.sigma_hydra)^2 + tau1^2 + tau2^2);
%     else
%         cost = cost + opts.lambda * (sigma^2 + tau1^2 + tau2^2);
%     end
% end

% --- Penalty for nonphysical dipolar relaxation loss (imag < 0)
relax_imag = imag(eps_model) - sigma ./ (omega * eps0);  % dipolar loss only
neg_relax = max(0, -relax_imag);  % zero or positive values where imag < 0
if isfield(opts, 'relaxPenaltyWeight')
    relax_penalty_weight = opts.relaxPenaltyWeight;
else
    relax_penalty_weight = 10;  % default value
end
penalty_relax = relax_penalty_weight * sum(neg_relax.^2);

% --- Regularization (optional)
reg_penalty = 0;
if opts.useRegularization
    if isfield(opts, 'sigma_hydra') && weight_hydra > 0
        reg_penalty = opts.lambda * ((sigma - opts.sigma_hydra)^2 + tau1^2 + tau2^2);
    else
        reg_penalty = opts.lambda * (sigma^2 + tau1^2 + tau2^2);
    end
end

% --- Final Cost ---
cost = vna_cost + weight_hydra * hydra_cost + penalty_relax + reg_penalty;


% --- Plot (optional) ---
if opts.plotModel
    figure(99); clf
    subplot(2,1,1); hold on;
    plot(freq, real(eps_meas), 'ko-', 'DisplayName', 'Real(\epsilon_{meas})');
    plot(freq, real(eps_model), 'r.-', 'DisplayName', 'Real(\epsilon_{model})');
    yline(real(eps_hydra), 'k--', 'DisplayName', 'HydraGO Real');
    ylabel('\epsilon'''); grid on; legend show

    subplot(2,1,2); hold on;
    plot(freq, imag(eps_meas), 'ko-', 'DisplayName', 'Imag(\epsilon_{meas})');
    plot(freq, imag(eps_model), 'b.-', 'DisplayName', 'Imag(\epsilon_{model})');
    yline(imag(eps_hydra), 'k--', 'DisplayName', 'HydraGO Imag');
    ylabel('\epsilon'''''); xlabel('Freq [Hz]'); grid on; legend show
    sgtitle(sprintf('Total Cost: %.4f', cost));
end

end

function eps_model = compute_dual_cole_cole(freq, p)
% Returns the complex permittivity for a dual Cole–Cole model at freq.

eps_inf  = p(1); eps_s1 = p(2); eps_s2 = p(3);
tau1     = p(4); alpha1 = p(5);
tau2     = p(6); alpha2 = p(7);
sigma    = p(8);

eps0 = 8.854187817e-12;
omega = 2 * pi * freq(:);

delta_eps1 = eps_s1 - eps_inf;
delta_eps2 = eps_s2 - eps_s1;

term1 = delta_eps1 ./ (1 + (1i * omega * tau1).^(1 - alpha1));
term2 = delta_eps2 ./ (1 + (1i * omega * tau2).^(1 - alpha2));
cond  = 1i * sigma ./ (omega * eps0);

eps_model = eps_inf + term1 + term2 + cond;
end

function [c, ceq] = nonlin_constraints_eps(params, freq)

% Compute imaginary part of permittivity at all frequencies
eps_model = compute_dual_cole_cole(freq, params);
imag_eps = imag(eps_model);

% Inequality constraint: imag_eps >= 0 → -imag_eps <= 0
c = -imag_eps;

% No equality constraints
ceq = [];
end

function [c, ceq] = dual_cole_cole_relaxation_constraint(params, freq)
% Nonlinear inequality constraint: relaxation loss must be ≥ 0

% Unpack parameters
eps_inf = params(1); eps_s1 = params(2); eps_s2 = params(3);
tau1 = params(4);    alpha1 = params(5);
tau2 = params(6);    alpha2 = params(7);
sigma = params(8);

eps0 = 8.854187817e-12;
omega = 2 * pi * freq;

% Dual Cole–Cole relaxation (no conductivity term)
delta_eps1 = eps_s1 - eps_inf;
delta_eps2 = eps_s2 - eps_s1;
term1 = delta_eps1 ./ (1 + (1i * omega * tau1).^(1 - alpha1));
term2 = delta_eps2 ./ (1 + (1i * omega * tau2).^(1 - alpha2));

eps_relax = imag(term1 + term2);

% Inequality constraint: relax loss ≥ 0
c = -eps_relax;  % all values must be ≤ 0 for constraint to pass
ceq = [];        % no equality constraints
end

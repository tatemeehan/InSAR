function result = fit_cole_cole_eps(freq, eps_complex)
%FIT_COLE_COLE_EPS Fits a Cole-Cole model with conductivity to complex permittivity data
%
% result = fit_cole_cole_eps(freq, eps_complex)
%
% Inputs:
%   freq         - frequency vector [Hz]
%   eps_complex  - complex permittivity vector (same size as freq)
%
% Output:
%   result: struct with fields
%       .params     - fitted parameters: [eps_s, eps_inf, tau, alpha, sigma]
%       .model_eps  - modeled complex permittivity
%       .rmse       - root-mean-square error

    % Constants
    eps0 = 8.854187817e-12; % vacuum permittivity [F/m]

    % Initial guess: [eps_s, eps_inf, tau, alpha, sigma]
    eps_real = real(eps_complex);
    eps_imag = imag(eps_complex);
    % eps_s_guess = max(eps_real);
    % eps_inf_guess = min(eps_real);
    % tau_guess = 1 / (2*pi*mean(freq));
    % alpha_guess = 0.05;
    % sigma_guess = 0.01;
    eps_s_guess = 1;
    eps_inf_guess = 5;
    tau_guess = 2 / (2*pi*mean(freq));
    alpha_guess = 0.0;
    sigma_guess = 0.01;
    p0 = [eps_s_guess, eps_inf_guess, tau_guess, alpha_guess, sigma_guess];

    % Bounds
    % lb = [1, 1, 1e-12, 0,     0];
    % ub = [100, 100, 1e-3,  0.5,   10];
    lb = [5, 1, 1e-12, 0, 0];
    ub = [50, 20, 1e-7, 1, 5];

    % Cost function for complex fit
    cost_fun = @(p) compute_error(p, freq, eps_complex, eps0);

    % Fit using nonlinear least squares
    % options = optimoptions('lsqnonlin','Display','iter','TolFun',1e-9,'TolX',1e-9);
    options = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 2000, ...
    'TolFun', 1e-9, ...
    'TolX', 1e-9);
    [p_fit, ~, residual] = lsqnonlin(cost_fun, p0, lb, ub, options);

    % Reconstruct fitted model
    eps_fit = cole_cole_model(p_fit, freq, eps0);
    rmse = sqrt(mean(abs(eps_fit - eps_complex).^2));

    % Output structure
    result.params = p_fit;
    result.model_eps = eps_fit;
    result.rmse = rmse;

end

function eps_model = cole_cole_model(p, freq, eps0)
    eps_s = p(1);
    eps_inf = p(2);
    tau = p(3);
    alpha = p(4);
    sigma = p(5);
    omega = 2 * pi * freq;
    jw_tau_alpha = (1 - (1i * omega * tau).^(1 - alpha));
    eps_model = eps_inf + (eps_s - eps_inf) ./ jw_tau_alpha + 1i.*sigma ./ (omega * eps0);
end

function err = compute_error(p, freq, eps_meas, eps0)
    eps_model = cole_cole_model(p, freq, eps0);
    err = [real(eps_model - eps_meas); imag(eps_model - eps_meas)];
    % err = sqrt(mean(real(eps_model - eps_meas).^2 + imag(eps_model - eps11_meas).^2));
    % err = [ ...
    % (real(eps_model) - real(eps_meas)) ./ abs(real(eps_meas)); 
    % (imag(eps_model) - imag(eps_meas)) ./ abs(imag(eps_meas))];

end

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
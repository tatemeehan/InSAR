function eps_clay = clay_permittivity_colecole(f, eps_real_clay, clayParams)
% CLAY_PERMITTIVITY_COLECOLE
% Single Cole–Cole term for clay.
%
% clayParams: struct
%   .eps_inf_clay : high-freq permittivity of clay
%   .tau_clay     : relaxation time (s)
%   .alpha_clay   : broadening parameter

    omega = 2*pi*f;

    eps_inf_clay = clayParams.eps_inf_clay;
    tau_clay     = clayParams.tau_clay;
    alpha_clay   = clayParams.alpha_clay;

    delta_eps_clay = eps_real_clay - eps_inf_clay;

    % jomega_tau_clay = (1i * omega * tau_clay).^(1 - alpha_clay);
    % eps_clay = eps_inf_clay + delta_eps_clay ./ (1 + jomega_tau_clay);
    eps_clay = eps_inf_clay + delta_eps_clay ./ (1 - (1i*omega*tau_clay).^(1-alpha_clay));

end

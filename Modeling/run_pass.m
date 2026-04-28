function out = run_pass(theta_i, f, snowState, soil_eps_g0, pars)
%RUN_PASS Build eps_snow from snowState then call benchmark.

% --- Build wet snow permittivity (you already have this pipeline) ---
eps_snow = wetsnow_permittivity_tinga73_extended(f, snowState.Tk, snowState.rho, snowState.lwc);

% --- Call benchmark ---
% out = insar_forward_snow_soil_benchmark(theta_i, f, snowState.Hs, eps_snow, soil_eps_g0, pars);
% Call 4 component Model
out = insar_forward_snow_soil_4component(theta_i, f, snowState.Hs, eps_snow, soil_eps_g0, pars);

% --- Attach for logging ---
out.snowState = snowState;
out.eps_snow  = eps_snow;
end

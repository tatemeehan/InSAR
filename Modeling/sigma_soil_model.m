function sigma_eff = sigma_soil_model(theta_liq, phi, frac_clay, sP)
% theta_liq: liquid water content (m3/m3)
% phi      : porosity (SWC)
% frac_clay: 0..1
% sP fields: sigma0, k_clay, sigma_pore, m, sigma_floor

if ~isfield(sP,'sigma_floor'), sP.sigma_floor = 0; end
if ~isfield(sP,'sigma0'),      sP.sigma0      = 0.03; end % dS/m
if ~isfield(sP,'k_clay'),      sP.k_clay      = 0.25; end % dS/m per clay frac
if ~isfield(sP,'sigma_pore'),  sP.sigma_pore  = 0.35; end % dS/m
if ~isfield(sP,'m'),           sP.m           = 1.5;  end

Se = max(theta_liq ./ max(phi, realmin), 0);  % effective saturation
sigma_surf = sP.sigma0 + sP.k_clay .* frac_clay;
sigma_conn = sP.sigma_pore .* (Se .^ sP.m);

sigma_eff = max(sP.sigma_floor + sigma_surf + sigma_conn, 0); % dS/m
end

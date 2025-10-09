%% BiStatic InSAR SWE Retrieval
addpath(genpath('E:\MCS\MCS_LiDAR_GPR'))

LWC = 0;
density = 300;
% 1.3 GHz
tauParams = [2.70604364852145e-10, -8.95922225296886e-13, 1.10036955414893e-15 ,-4.79169239741482e-19];
tauFun = @(p,rho) p(1) + p(2).*rho + p(3).*rho.^2 + p(4).*rho.^3;
tau = tauFun(tauParams,density);      % Relaxation Time

tingaParams = struct( ...
    'eps_s', 87.5, ...
    'eps_inf', 4.9, ...
    'tau', tau, ...
    'alpha', 0, ...
    'sigma', 0.001);


% Loop Over InSAR Penetration Pairs
for kk = 2:numel(insarData)
    if insarData(kk).penetrationValid
        % Inputs you already have:
        dz        = penetrationDeramp;%-insarData(kk).penetration;             % vertical height change [m]
        theta_inc = deg2rad(geomData.sarGeometry{insarData(kk).geomIdx}.incidence);  % local incidence (LOS vs normal)
        slope      = deg2rad(demData.slope);               % slope from vertical
        % Permittivity Modeling
        fHz = params.f.*10^9;
        eps_w = water_permittivity_colecole(fHz, tingaParams);
        depol = 0.2262;
        Perm = real(wetsnow_permittivity_tinga73(fHz,273.15,density,LWC, depol, @ice_permittivity_maetzler06, eps_w));
        V = (params.c./sqrt(Perm));
        n = sqrt(Perm); % snow refractive index

        % Normal thickness (Incidence geometry)
        hs_inc = dz ./ max(cos(slope), 1e-6);

        % Refraction correction
        theta_s = asin( sin(theta_inc) ./ n );
        K = (n.*cos(theta_s)) ./ max(cos(theta_inc), 1e-6);

        % Physical snow thickness (normal direction)
        snowThickness = hs_inc ./ K;

        % vertical depth (z):
        snowDepth = snowThickness .* cos(slope);
        swe = (snowDepth.*density)./1000;
    end
end
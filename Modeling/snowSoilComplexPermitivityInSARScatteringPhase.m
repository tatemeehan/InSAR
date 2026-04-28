% Water parameters (Kaatze-polynomial based)
epssCoefs = [88.119019607843190,-0.484234778121778,0.004464654282766,-2.393188854489220e-05]';
tauCoefs = [1.736936274509802e-11,-4.823117475060188e-13,5.488493292053644e-15,-2.366013071895416e-17]';
wparams.epssCoefs   = epssCoefs;   % your 4-coef vec
wparams.tauCoefs    = tauCoefs;
wparams.eps_inf     = 4.8;
wparams.alpha_free  = 0.01;
wparams.alpha_bound = 0.2;
wparams.frac_bound  = 0.3;    % 30% of Δε in bound mode (tune later)

tauBW.free_scale  = 1.0;
tauBW.bound_scale = 5.0;      % bound water ~5× slower than free
tauBW.B_bound     = 0.02;     % VWC sensitivity (tune)
tauBW.max_tau     = 5e-10;

% Simple soil defaults, or pass a struct
soilParams.eps_real_soil = 4.0;
soilParams.eps_real_clay = 25;
soilParams.frac_clay     = 0.1;

% Clay Parameters
clayParams.eps_inf_clay = 5;
clayParams.tau_clay     = 5e-10;
clayParams.alpha_clay   = 0.15;

% freqs = [50e6, 1.3e9, 3.2e9];  % 50 MHz, L-band, S-band
frequency = 1.3e9;
% Soil data
VWC    = 0.202666666666667;      % N×1
SWC    = 0.556753;               % scalar porosity
EC_dSm = 0.2933;       % N×1, dS/m
Tc     = 17.33;      % N×1

[eps_soil, components] = soilPermittivity_advanced( ...
    VWC, SWC, EC_dSm.*0.1, soilParams, Tc, 1.3e9, wparams, tauBW, clayParams)
%% Snow 
addpath(genpath('E:\MCS\Tinga73\tinga_extended_model_bundle'))
rho = 350;
lwc = 0.015;
tempSnow = 273.15;

   eps_w = water_permittivity_maetzler87(frequency, 273.15);

    eps_snow = wetsnow_permittivity_tinga73_extended(frequency, tempSnow, rho, lwc, ...
        @ice_permittivity_maetzler06, eps_w)

%% InSAR Forward Model
% Given by you:
f = 1.3e9;
eps_soil = 9.3143 + 1.9261i;
eps_snow = 1.8003 + 0.0183i;

theta_i = 40;   % deg
Hs = 1.0;       % m (example; change as desired)

pars = struct();
pars.Ls = 5.0;      % m
pars.Lg = 0.05;     % m
pars.As = 0.05;     % snow weight (Level 1)
pars.Ag = 1.00;     % soil weight (Level 1)

out = insar_forward_snow_soil_benchmark(theta_i, f, Hs, eps_snow, eps_soil, pars);

fprintf('kz_snow = %.4f + %.4fi  [rad/m]\n', real(out.kzs), imag(out.kzs));
fprintf('kz_soil = %.4f + %.4fi  [rad/m]\n', real(out.kzg), imag(out.kzg));
fprintf('alpha_s = %.4f 1/m, alpha_g = %.4f 1/m\n', out.alpha_s, out.alpha_g);
fprintf('phi = %.4f rad (%.2f deg)\n', out.phi, rad2deg(out.phi));
fprintf('|Es|/|Eg| = %.3f\n', abs(out.Es)/max(abs(out.Eg),eps));

%% Full Test Run
% --- Fixed settings ---
f = 1.3e9;
theta_i = 40;
Hs = 1.0;

% Fixed snow state (your existing call)
rho = 350;
lwc = 0.015;
tempSnow = 273.15;
eps_w = water_permittivity_maetzler87(f, 273.15);
eps_snow = wetsnow_permittivity_tinga73_extended(f, tempSnow, rho, lwc, ...
    @ice_permittivity_maetzler06, eps_w);

% Benchmark scattering params
pars = struct();
pars.Ls = 5.0;
pars.Lg = 0.05;
pars.As = 0.05;

% Turn on roughness kernel
pars.use_roughness = true;
pars.rough = struct();
pars.rough.sigma_h = 0.01;   % 1 cm RMS
pars.rough.ell     = 0.10;   % 10 cm corr length
pars.rough.A0      = 1.0;    % overall soil scale
pars.rough.Cdiff   = 0;    % diffuse fraction strength
pars.rough.psi     = 2*pi*rand; % Level-2 random phase

% If roughness is enabled, pars.Ag is ignored (or can remain as fallback)
pars.Ag = 1.0;

% Soil fixed state inputs (from your snippet)
SWC    = 0.556753;
EC_dSm = 0.2933;
Tc     = 17.33;

sigma_bulk = EC_dSm * 0.1;  % S/m

% --- Sweep VWC ---
VWC_grid = linspace(0.01, 0.30, 60).';

phi = zeros(size(VWC_grid));
alpha_g = zeros(size(VWC_grid));
beta_g  = zeros(size(VWC_grid));
ratio_EsEg = zeros(size(VWC_grid));

Es = zeros(size(VWC_grid));
Eg = zeros(size(VWC_grid));
E  = zeros(size(VWC_grid));
R12 = zeros(size(VWC_grid));   % optional diagnostic
Ag = R12; Ts = R12; Dg = Ag;

for i = 1:numel(VWC_grid)
    VWC = VWC_grid(i);

    [eps_soil, comps] = soilPermittivity_advanced( ...
        VWC, SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams);

    out = insar_forward_snow_soil_benchmark(theta_i, f, Hs, eps_snow, eps_soil, pars);

    phi(i) = out.phi;
    alpha_g(i) = out.alpha_g;
    beta_g(i)  = out.beta_g;
    ratio_EsEg(i) = abs(out.Es) / max(abs(out.Eg), eps);

    Es(i) = out.Es;
    Eg(i) = out.Eg;
    E(i)  = out.E;

    % Optional: reflection coefficient phase/magnitude diagnostic
    R12(i) = (out.kzs - out.kzg) / (out.kzs + out.kzg);
    Ag(i) = out.Ag;
    Ts(i) = out.Ts;
    Dg(i) = out.Dg;

    % out_store(i) = out;
end
% Amplitude Centroid Depth and Phase Center Depth
% a = 1./pars.Lg+2.*alpha_g; b = 2.*beta_g;
a = real(Dg); b = imag(Dg);
zphi_g = a./(a.^2+b.^2); zamp = 1./a;

phi_pre = angle(Ts .* Ag);     % wrapped
phi_int = angle(1 ./ Dg);      % wrapped

phi_sum_wrapped = wrapToPi(phi_pre + phi_int);
phi_num_wrapped = angle(Eg);

phase_err = wrapToPi(phi_num_wrapped - phi_sum_wrapped);
fprintf('max wrapped phase error = %.3g rad\n', max(abs(phase_err)));
phi_rec = unwrap(phi_pre + phi_int);

% For plotting relative to first VWC
phi_num_u = unwrap(phi_num_wrapped);
phi_pre_u = unwrap(phi_pre);                 % unwrap each separately for smooth curves
phi_int_u = unwrap(phi_int);

% then compare differences (not absolute branches)
dphi_num = rad2deg(phi_num_u - phi_num_u(1));
dphi_pre = rad2deg(phi_pre_u - phi_pre_u(1));
dphi_int = rad2deg(phi_int_u - phi_int_u(1));

% Reference at driest point
phi_ref = phi(1);
dphi = angle(exp(1i*(phi - phi_ref)));   % wrapped relative phase

% Useful diagnostic depth scale
delta_g = 1 ./ alpha_g;

% --- Quick prints ---
fprintf('VWC=%.3f eps_g=%s\n', VWC_grid(1), num2str(soilPermittivity_advanced(VWC_grid(1), SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams)));
fprintf('VWC=%.3f eps_g=%s\n', VWC_grid(end), num2str(soilPermittivity_advanced(VWC_grid(end), SWC, sigma_bulk, soilParams, Tc, f, wparams, tauBW, clayParams)));

%% --- Plots ---
figure; 
subplot(1,3,1)
plot(VWC_grid, rad2deg(dphi)); grid on;
xlabel('VWC (m^3/m^3)'); ylabel('\Delta\phi relative to VWC=0.01 (deg)');
title('Soil-moisture-driven phase change (benchmark model)');

% figure; 
subplot(1,3,2)
plot(VWC_grid, delta_g); grid on;
xlabel('VWC (m^3/m^3)'); ylabel('\delta_g = 1/\alpha_g (m)');
title('Soil effective depth scale vs VWC');

% figure; 
subplot(1,3,3)
plot(VWC_grid, ratio_EsEg); grid on;
xlabel('VWC (m^3/m^3)'); ylabel('|E_s|/|E_g|');
title('Relative snow-vs-soil contribution vs VWC');


% Complex-plane trajectories
figure; 
subplot(1,3,1)
plot(real(Eg), imag(Eg), '-o'); grid on; axis equal;
xlabel('Re(E_g)'); ylabel('Im(E_g)'); title('Soil phasor E_g as VWC varies');

% figure;
subplot(1,3,2)
plot(real(Es), imag(Es), '-o'); grid on; axis equal;
xlabel('Re(E_s)'); ylabel('Im(E_s)'); title('Snow phasor E_s (should be ~constant)');

% figure; 
subplot(1,3,3)
plot(real(E), imag(E), '-o'); grid on; axis equal;
xlabel('Re(E)'); ylabel('Im(E)'); title('Total phasor E = E_s + E_g as VWC varies');

% Label a few points so you can connect VWC to location
hold on;
ix = round(linspace(1,numel(VWC_grid),6));
text(real(E(ix)), imag(E(ix)), compose('%.2f', VWC_grid(ix)));


phi_s = unwrap(angle(Es));
phi_g = unwrap(angle(Eg));
phi_E = unwrap(angle(E));

figure; plot(VWC_grid, rad2deg(phi_s - phi_s(1)), ...
             VWC_grid, rad2deg(phi_g - phi_g(1)), ...
             VWC_grid, rad2deg(phi_E - phi_E(1)));
grid on; xlabel('VWC'); ylabel('\Delta phase (deg)');
legend('\Delta\phi_s (snow)','\Delta\phi_g (soil)','\Delta\phi_{total}');
title('Who is causing the phase change?');

figure; hold on; axis equal; grid on;
quiver(zeros(size(Es)),zeros(size(Es)),real(Eg),imag(Eg),0,'r','LineWidth',2)
quiver(zeros(size(Es)),zeros(size(Es)),real(Es),imag(Es),0,'b','LineWidth',2)
quiver(zeros(size(Es)),zeros(size(Es)),real(E), imag(E), 0,'k','LineWidth',2)
legend('E_g','E_s','E_{tot}')
xlabel('Re'); ylabel('Im')

% Plot amplitude and phase depth centroid
figure();plot(VWC_grid,zamp,'k');hold on; plot(VWC_grid,zphi_g,'r'); plot(VWC_grid, delta_g./2,'b')
grid on; grid minor;
xlabel('VWC (m^3/m^3)')
ylabel('Depth (m)')
legend('z_{amp}', 'z_{\phi}','\delta_z')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)

%% Soil Contribution Phase Plots
figure; plot(VWC_grid,dphi_num,'k',VWC_grid,dphi_pre,'b--',VWC_grid,dphi_int,'r-.');
grid on; xlabel('VWC'); ylabel('\Delta phase (deg)');
legend('\Delta\phi_g total','\Delta\phi_{pre}=\angle(T_sA_g)','\Delta\phi_{int}=\angle(1/D_g)');
title('Soil phase decomposition');

dphi_R12 = rad2deg(unwrap(angle(R12)) - angle(R12(1)));
figure; plot(VWC_grid,dphi_R12); grid on;
xlabel('VWC'); ylabel('\Delta \angle R_{12} (deg)');
title('Reflection coefficient phase change vs VWC');

b = imag(Dg);             % this is 2*beta_g
phi_int = unwrap(angle(1./Dg));

% finite difference equivalent depth (meters)
z_eq = -gradient(phi_int, b);   % dphi/db

figure; plot(VWC_grid, z_eq); grid on;
xlabel('VWC'); ylabel('z_{eq,int} (m)');
title('Equivalent depth from integral phase term');
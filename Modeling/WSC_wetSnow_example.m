%% Weastern Snow Conference Abstract Example
clear; close all; clc;
addpath('E:\MCS\Tinga73\tinga_extended_model_bundle')
%% Model of Wet Snow Perturbations over Wet Soil
rng(1);   % make MC repeatable on demand

% --- Radar settings ---
f = 1.3e9;
theta_i = 35;

% --- Soil baseline permittivity input (only used if profile disabled);
% if profile enabled, benchmark uses eps_g_surf from soil_profile_integral anyway
eps_g0 = 5 + 1i*0.5;  % placeholder; your pipeline already builds this

% --- Pars baseline (your existing pars) ---
pars = struct();
pars.output_phase_convention = "data_matched";      % or "native"
pars.Ls = 5.0;
pars.Lg = 0.05;
pars.As = 0.05;

pars.snow_surface.enable = true;
pars.snow_surface.sigma_h = 0.01;  % m
pars.snow_surface.ell     = 0.10;  % m
pars.snow_surface.A0      = 1.0;
pars.snow_surface.Cdiff   = 0.2;
pars.snow_surface.psi     = 0;     % diffuse phase
pars.snow_surface.psi0    = 0;     % coherent phase knob (optional)
pars.snow_surface.Xpol    = 0.02;  % optional
pars.snow_surface.psi_x   = 0;

pars.soil_surface.enable = true;   % uses your existing soil_roughness_kernel knobs

pars.soil_volume.enable = true;
pars.soil_volume.Ag     = 1.0;     % soil volume strength (not interface Fresnel)

% pars.apply_TasTsa_to_snowvol = true;  % default ON in code above

pars.use_roughness = true;
pars.rough.sigma_h = 0.02;     % try 0.01–0.10
pars.rough.ell     = 0.10;
pars.rough.A0      = 1.0;
pars.Ag = 1.0; % ignored when use_roughness=true
pars.rough.Cdiff   = 0.2;
pars.rough.psi     = 2*pi*rand;
pars.rough.Xpol    = 0.1;      % cross-pol leakage knob
pars.rough.psi_x   = 2*pi*rand;

pars.soil_surface = pars.rough;       % copy kernel params
pars.soil_surface.enable = true;      % ensure enable survives overwrite

pars.pol_mode = "matrix";
% Polarization
if strcmp(pars.pol_mode,"single")
pars.pol = 'HH';  % uses TE
% or
% pars.pol = 'VV';  % uses TM
% HV or VH
% pars.pol = 'HV';
end

% --- Soil profile (use what you already have) ---
pars.soil_profile = struct();
pars.soil_profile.enable = true;
pars.soil_profile.zmax = 0.50;
pars.soil_profile.dz   = 0.001;

pars.soil_profile.type   = 'exp';
pars.soil_profile.VWC0   = 0.20;
pars.soil_profile.VWCinf = 0.10;
pars.soil_profile.z0     = 0.07;

pars.soil_profile.SWC = 0.55;
pars.soil_profile.sigma_bulk = 0.03;   % S/m (only used if your frozen physical model overwrites sigma)
% Exponential Temperture Model
pars.soil_profile.Tc_profile = struct('type','exp','Tc0',-1.0,'Tc_inf',2.0,'zT',0.10);

pars.soil_profile.freeze.enable = false; % start unfrozen for wet snow showcase
% Uncomment if True
% --- Enable freezing wrapper ---
% pars.soil_profile.freeze = struct( ...
%     'enable',  true, ...
%     'T0',      0.0, ...
%     'dT',      0.5, ...
%     'p_sigma', 2, ...
%     'mode', 'physical',...
%     'eta_sigma', 2,...
%     'sigma_floor', 0.0);
% % 1) Logistic
% % pars.soil_profile.freeze.uwf = struct('mode','logistic','T0',0,'dT',0.5,'fmin',0.02);
% 
% % 2) Temperature-only (fewer knobs)
% % pars.soil_profile.freeze.uwf = struct('mode','temp_only','Tm',0,'gamma',1.2,'fmin',0.02);
% 
% % 3) Texture-based (uses soilParams.frac_clay by default)
% pars.soil_profile.freeze.uwf = struct('mode','texture','fmin',0.02);

% --- Temporal decorrelation knobs ---
pars.temporal.rho_soil = 0.90;
pars.temporal.rho_snow = 0.80;

% --- Snow cross-pol knobs (your new snow matrix behavior) ---
pars.snow.Xpol = 0.05;
pars.snow.psi_x = 2*pi*rand;
pars.snow.VV_scale = 1.0;

% Soil Parameters
% Water parameters (Kaatze-polynomial based)
epssCoefs = [88.119019607843190,-0.484234778121778,0.004464654282766,-2.393188854489220e-05]';
tauCoefs = [1.736936274509802e-11,-4.823117475060188e-13,5.488493292053644e-15,-2.366013071895416e-17]';
wparams.epssCoefs   = epssCoefs;   % your 4-coef vec
wparams.tauCoefs    = tauCoefs;
wparams.eps_inf     = 4.8;
wparams.alpha_free  = 0.01;
wparams.alpha_bound = 0.2;
wparams.frac_bound  = 0.3;  

% Bound Water Relaxation Parameters
tauBW.free_scale  = 1.0;
tauBW.bound_scale = 5.0;      % bound water ~5× slower than free
tauBW.B_bound     = 0.02;     % VWC sensitivity (tune)
tauBW.max_tau     = 5e-10;

% Simple soil defaults
soilParams.eps_real_soil = 4.0;
soilParams.eps_real_clay = 25;
soilParams.frac_clay     = 0.1;

% Clay Parameters
clayParams.eps_inf_clay = 5;
clayParams.tau_clay     = 5e-10;
clayParams.alpha_clay   = 0.15;

pars.soil_profile.soilParams = soilParams;
pars.soil_profile.wparams    = wparams;
pars.soil_profile.tauBW      = tauBW;
pars.soil_profile.clayParams = clayParams;

% Snow Elevation Change
% Default uses Snow Depth Differnce
% Optionally a z_rel for each pass can be used e.g.
% pars.snow_surface.z_rel = zSnowSurf23 - zSnowFree;
pars.snow_surface.include_height_phase = true;
pars.snow_surface.height_phase_sign = +1;

% ---------------- Snow definitions ----------------
snow1 = struct('Hs',1.00,'rho',350,'lwc',0.015,'Tk',273.15);
snow2 = struct('Hs',1.00,'rho',350,'lwc',0.060,'Tk',273.15);  % wetter pass

pars1 = pars; pars2 = pars;
pars1.snowState = snow1;
pars2.snowState = snow2;

% --- Deterministic (mean-field) outputs ---
o1 = run_pass(theta_i, f, snow1, eps_g0, pars1);
o2 = run_pass(theta_i, f, snow2, eps_g0, pars2);

% dphi_det = angle(ensure_2x2(o1.E) .* conj(ensure_2x2(o2.E)));
dphi_det= angle((o1.E) .* conj((o2.E)));

% --- Monte Carlo coherence ---
Nmc = 250;
% coh = coherence_mc(theta_i,f, snow1.Hs, o1.eps_snow, eps_g0, pars1, pars2, Nmc); 
coh = coherence_mc_4component(theta_i,f, snow1.Hs, o1.eps_snow, eps_g0, pars1, pars2, Nmc); 

% NOTE: if your coherence_mc signature still expects Hs, eps_snow, update it to only use pars.snowState

% --- Print summary ---
lab = ["HH","HV";"VH","VV"];
fprintf('\n--- Deterministic Δphi (deg) ---\n');
disp(array2table(rad2deg(dphi_det), 'VariableNames', {'H*','V*'}, 'RowNames', {'H','V'}));

fprintf('--- Coherence ---\n');
fprintf('|gamma_HH|=%.3f, ang=%.1f deg\n', abs(coh.HH), rad2deg(angle(coh.HH)));
fprintf('|gamma_VV|=%.3f, ang=%.1f deg\n', abs(coh.VV), rad2deg(angle(coh.VV)));
fprintf('|gamma_HV|=%.3f, ang=%.1f deg\n', abs(coh.HV), rad2deg(angle(coh.HV)));
fprintf('|gamma_VH|=%.3f, ang=%.1f deg\n', abs(coh.VH), rad2deg(angle(coh.VH)));

fprintf('--- Magnitude partition (HH) ---\n');
% fprintf('|Es1|=%.3g  |Eg1|=%.3g  (|Es|/|Eg|=%.3g)\n', abs(o1.Es(1,1)), abs(o1.Eg(1,1)), abs(o1.Es(1,1))/max(abs(o1.Eg(1,1)),realmin));
% fprintf('|Es2|=%.3g  |Eg2|=%.3g  (|Es|/|Eg|=%.3g)\n', abs(o2.Es(1,1)), abs(o2.Eg(1,1)), abs(o2.Es(1,1))/max(abs(o2.Eg(1,1)),realmin));

fprintf('|Eas1|=%.3g |Es1|=%.3g |Esg1|=%.3g |Egv1|=%.3g\n', ...
    abs(o1.Eas(1,1)), abs(o1.Es(1,1)), abs(o1.Esg(1,1)), abs(o1.Egv(1,1)));

Eg1 = o1.Esg + o1.Egv;
fprintf('|Es1|/|Eg1|=%.3g   |Eas1|/|E|=%.3g\n', ...
    abs(o1.Es(1,1))/max(abs(Eg1(1,1)),realmin), abs(o1.Eas(1,1))/max(abs(o1.E(1,1)),realmin));

%% Case Snow SWE Change
lwc = linspace(0,0.1,1000);
hs = linspace(1,2,1000);
dphi_det = zeros(2,2,numel(lwc,1));
dphi_det_snow = zeros(2,2,numel(lwc,1));
dphi_det_soil = zeros(2,2,numel(lwc,1));
dphi_Ts_snow = zeros(numel(lwc),1);
dphi_soil_res = dphi_det_soil;
dphi_det_air = dphi_det_soil;
deltaSWE = zeros(numel(lwc,1),1);
TempK = 270;
for kk = 1:numel(lwc)
% ---------------- Snow definitions ----------------
snow1 = struct('Hs',hs(1),'rho',300,'lwc',lwc(1),'Tk',TempK);
snow2 = struct('Hs',hs(kk),'rho',300,'lwc',lwc(1),'Tk',TempK);  % wetter pass
deltaSWE(kk) = snow2.Hs.*snow2.rho-snow1.Hs.*snow1.rho;
pars1 = pars; pars2 = pars;
pars1.snowState = snow1;
pars2.snowState = snow2;

% --- Deterministic (mean-field) outputs ---
o2 = run_pass(theta_i, f, snow1, eps_g0, pars1);
o1 = run_pass(theta_i, f, snow2, eps_g0, pars2);
% 
% dphi_det(kk) = angle(ensure_2x2(o1.E) .* conj(ensure_2x2(o2.E)));
% dphi_det_snow(kk) = angle(ensure_2x2(o1.Es) .* conj(ensure_2x2(o2.Es)));
% dphi_det_soil(kk) = angle(ensure_2x2(o1.Eg) .* conj(ensure_2x2(o2.Eg)));

% dphi_det(:,:,kk) = angle((o1.E) .* conj((o2.E)));
% dphi_det_snow(:,:,kk) = angle((o1.Es) .* conj((o2.Es)));
% dphi_Ts_snow(kk) = angle(o1.Ts.*conj(o2.Ts));
% dphi_det_soil(:,:,kk) = angle((o1.Eg) .* conj((o2.Eg)));

E1 = o1.E;   E2 = o2.E;
Es1 = o1.Es; Es2 = o2.Es;
Eas1 = o1.Eas; Eas2 = o2.Eas;

% NEW: ground/soil total in 4-component model
Eg1 = o1.Esg + o1.Egv;
Eg2 = o2.Esg + o2.Egv;

Eg1_res = Eg1 ./ o1.Ts_eff;
Eg2_res = Eg2 ./ o2.Ts_eff;

dphi_soil_res(:,:,kk) = angle( Eg1_res .* conj(Eg2_res) );

dphi_det(:,:,kk)      = angle(E1  .* conj(E2));
dphi_det_air(:,:,kk) = angle(Eas1 .* conj(Eas2));
dphi_det_snow(:,:,kk) = angle(Es1 .* conj(Es2));
dphi_det_soil(:,:,kk) = angle(Eg1 .* conj(Eg2));

dphi_Ts_snow(kk) = angle(o1.Ts_excess .* conj(o2.Ts_excess));
end
phz = squeeze(dphi_det(2,2,:));
phzAir = squeeze(dphi_det_air(2,2,:));
phzSnow = squeeze(dphi_det_snow(2,2,:));
phzSoil = squeeze(dphi_det_soil(2,2,:));
phzSoilTs = squeeze(dphi_soil_res(2,2,:));

% Oviesgarhan SWE
oviDSWE = 0:0.001:100;
oviDPhi = oveis_phi_wrapped(oviDSWE, theta_i, f, 300);
% Make a Conspiratory Figure
% figure(); 
% hold on; 
% plot(lwc.*100,phzSnow,'linewidth', 2, 'color',[0.5,0.5,0.5],"DisplayName","Snow \phi"); 
% plot(lwc.*100,phzSoil,'linewidth', 1.5, 'color',[0.5,0.5,0.5],'lineStyle',"--","DisplayName","Soil \phi");
% plot(lwc.*100,phz,'k','linewidth', 1.5, "DisplayName","Measured \phi");
% xlabel('Liquid Water Content Change (%)')
% ylabel('\Delta \phi (rad)')
% ylim([-pi pi])
% legend('Location','southwest')
% set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
% grid on; grid minor;

figure(); 
hold on; 
plot(deltaSWE,phzAir,'linewidth', 1.5, 'color',[0.25,0.25,0.25],'lineStyle',"--","DisplayName","Snow Interface \phi"); 
plot(deltaSWE,phzSnow,'linewidth', 1.5, 'color',[0.5,0.5,0.5],"DisplayName","Snow Volume \phi"); 
plot(deltaSWE,phzSoil,'linewidth', 1.5, 'color',[0.5,0.5,0.5],'lineStyle',"--","DisplayName","Snow Transmission \phi");
plot(deltaSWE,phz,'k','linewidth', 1.5, "DisplayName","Measured \phi");
% plot(oviDSWE,oviDPhi,'r',"DisplayName","Oveisgharan \phi")
% plot(deltaSWE,phzSoilTs)
xlabel('SWE Change (mm)')
ylabel('\Delta \phi (rad)')
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; grid minor;

figure();
hold on;
plot(oviDSWE,oviDPhi,'k',"DisplayName","Oveisgharan \phi")
plot(deltaSWE,dphi_Ts_snow,'--k',"DisplayName","Snow Transmission \phi")
plot(deltaSWE,phzSoilTs,'--r',"DisplayName","Ground Residual \phi")
xlabel('SWE Change (mm)')
ylabel('\Delta \phi (rad)')
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; grid minor;

figure();plot(deltaSWE,phzSoil,'k','linewidth',1.5,"DisplayName","dphi det soil \phi")
hold on; plot(deltaSWE,dphi_Ts_snow,'--r',"DisplayName","dphi Ts snow \phi")

xlabel('SWE Change (mm)')
ylabel('\Delta \phi (rad)')
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest')
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; grid minor;

figure(); subplot(1,3,1)
hold on;
plot(oviDSWE,oviDPhi,'k','linewidth', 2,"DisplayName","\Delta\phi Oveisgharan et al. (2024)")
plot(deltaSWE,dphi_Ts_snow,'linewidth', 2, 'color',[0.5,0.5,0.5],"DisplayName","\Delta\phi Snow Propagation")
% plot(deltaSWE,phzSoilTs,'--r',"DisplayName","Ground Residual \phi")
% title('Canonical Refraction Limit')
title('Refraction-Based Phase Limit')
xlabel('\DeltaSWE (mm)')
ylabel('\Delta\phi (rad)')
xticks([0,25,50,75,100])
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.75)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; %grid minor;
axis square
text(-0.1,1.0675,'(a)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');
subplot(1,3,2)
hold on;
plot(deltaSWE,phz,'k','linewidth', 2, "DisplayName","\Delta\phi Measured");
plot(deltaSWE,dphi_Ts_snow,'linewidth', 2, 'color',[0.5,0.5,0.5],"DisplayName","\Delta\phi Snow Propagation");
% plot(deltaSWE,phzAir,'linewidth', 1.5, 'color',[0.25,0.25,0.25],'lineStyle',"--","DisplayName","Snow Interface \phi"); 
plot(deltaSWE,phzSnow,'linewidth', 1.25, 'color',[0.5,0.5,0.5],'lineStyle',"--","DisplayName","\Delta\phi Snow Volume "); 
% plot(oviDSWE,oviDPhi,'r',"DisplayName","Oveisgharan \phi")
% plot(deltaSWE,phzSoilTs)
% title('Coherent Phasor Mixing')
title('Complex Phasor Summation')
xlabel('\DeltaSWE (mm)')
ylabel('\Delta\phi (rad)')
xticks([0,25,50,75,100])
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.75)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on;%grid minor;
axis square
text(-0.1,1.0675,'(b)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');


%% Constant SWE Change with Soil Change
test_dphi_det = zeros(2,2,3);
test_dphi_det_snow = zeros(2,2,3);
test_dphi_det_soil = zeros(2,2,3);
test_dphi_Ts_snow = zeros(3,1);
test_dphi_soil_res = dphi_det_soil;
test_dphi_det_air = dphi_det_soil;
% Recall Pars
for kk = 1:3
    pars1 = pars; pars2 = pars;

% ---------------- Snow definitions ----------------
snow1 = struct('Hs',1,'rho',300,'lwc',lwc(1),'Tk',TempK);
snow2 = struct('Hs',1.0833,'rho',300,'lwc',lwc(1),'Tk',TempK);  % wetter pass
deltaSWEsoil = snow2.Hs.*snow2.rho-snow1.Hs.*snow1.rho;

% Snow Pars for run
pars1.snowState = snow1;
pars2.snowState = snow2;

soilprof = struct();
soilprof.enable = true;
soilprof.zmax = 0.50;
soilprof.dz   = 0.001;
soilprof.type   = 'exp';
soilprof.SWC = 0.55;
soilprof.sigma_bulk = 0.03;
soilprof.z0 = 0.07;

soilprof.Tc_profile = struct('type','exp','Tc0',-2.0,'Tc_inf',2.0,'zT',0.10);

soilprof.soilParams = soilParams;
soilprof.wparams    = wparams;
soilprof.tauBW      = tauBW;
soilprof.clayParams = clayParams;
soilprof.freeze = struct('enable',false); % default unfrozen
soilprof.freeze.uwf = struct('mode','texture','fmin',0.02);

switch kk
  case 1 % Frozen
    soilprof.VWC0 = 0.20; soilprof.VWCinf = 0.10;
    soilprof.freeze.enable = true;
    soilprof.freeze.mode = 'physical';
    soilprof.freeze.T0 = 0.0; soilprof.freeze.dT = 0.5;
    soilprof.freeze.p_sigma = 2;
    soilprof.freeze.eta_sigma = 2;
    soilprof.freeze.sigma_floor = 0.0;

  case 2 % Moist (unfrozen)
    soilprof.VWC0 = 0.20; soilprof.VWCinf = 0.10;
    soilprof.freeze.enable = false;

  case 3 % Wet (unfrozen, higher VWC)
    soilprof.VWC0 = 0.30; soilprof.VWCinf = 0.15;
    soilprof.freeze.enable = false;
end

pars1.soil_profile = soilprof;
pars2.soil_profile = soilprof;  % identical soil within this kk
assert(isequaln(pars1.soil_profile, pars2.soil_profile), 'Soil differs between passes!');

% --- Deterministic (mean-field) outputs ---
o1 = run_pass(theta_i, f, snow1, eps_g0, pars1);
o2 = run_pass(theta_i, f, snow2, eps_g0, pars2);

% Compute Enery Weighted Effective Complex Permittivity
switch kk
    case 1 % Frozen
        eps_eff_frozen = eps_eff_energy(o1.profile.diag);
    case 2 % Moist (unfrozen)
        eps_eff_moist = eps_eff_energy(o1.profile.diag);
    case 3 % Wet (unfrozen, higher VWC)
        eps_eff_wet= eps_eff_energy(o1.profile.diag);
end
% 
% dphi_det(kk) = angle(ensure_2x2(o1.E) .* conj(ensure_2x2(o2.E)));
% dphi_det_snow(kk) = angle(ensure_2x2(o1.Es) .* conj(ensure_2x2(o2.Es)));
% dphi_det_soil(kk) = angle(ensure_2x2(o1.Eg) .* conj(ensure_2x2(o2.Eg)));

% dphi_det(:,:,kk) = angle((o1.E) .* conj((o2.E)));
% dphi_det_snow(:,:,kk) = angle((o1.Es) .* conj((o2.Es)));
% dphi_Ts_snow(kk) = angle(o1.Ts.*conj(o2.Ts));
% dphi_det_soil(:,:,kk) = angle((o1.Eg) .* conj((o2.Eg)));

E1 = o1.E;   E2 = o2.E;
Es1 = o1.Es; Es2 = o2.Es;
Eas1 = o1.Eas; Eas2 = o2.Eas;

% NEW: ground/soil total in 4-component model
Eg1 = o1.Esg + o1.Egv;
Eg2 = o2.Esg + o2.Egv;

Eg1_res = Eg1 ./ o1.Ts_eff;
Eg2_res = Eg2 ./ o2.Ts_eff;

dphi_soil_res = angle( Eg1_res .* conj(Eg2_res) );

test_dphi_det(:,:,kk)       = angle(E1  .* conj(E2));
test_dphi_det_air(:,:,kk)  = angle(Eas1 .* conj(Eas2));
test_dphi_det_snow(:,:,kk)  = angle(Es1 .* conj(Es2));
test_dphi_det_soil(:,:,kk)  = angle(Eg1 .* conj(Eg2));

test_dphi_Ts_snow(kk) = angle(o1.Ts_excess .* conj(o2.Ts_excess));
end

soilPhz = squeeze(test_dphi_det(2,2,:));
eps_eff = [eps_eff_frozen; eps_eff_moist; eps_eff_wet];

eps_r = real(eps_eff);
eps_i = imag(eps_eff);

subplot(1,3,3)
hold on;
yline(0,'k-','LineWidth',1);
bar(1,soilPhz(1),0.8, 'FaceColor',[0,0,0],"DisplayName","Frozen");
bar(2,soilPhz(2), 0.8,'FaceColor',[0.35,0.35,0.35],"DisplayName","Moist");
bar(3,soilPhz(3),0.8,'FaceColor',[0.65,0.65,0.65], "DisplayName","Wet");

% plot(deltaSWE,phzAir,'linewidth', 1.5, 'color',[0.25,0.25,0.25],'lineStyle',"--","DisplayName","Snow Interface \phi"); 
% plot(deltaSWE,phzSnow,'linewidth', 1.5, 'color',[0.5,0.5,0.5],'lineStyle',"--","DisplayName","\Delta \phi Snow Volume "); 
% plot(oviDSWE,oviDPhi,'r',"DisplayName","Oveisgharan \phi")
% plot(deltaSWE,phzSoilTs)
% title('Soil State Varies \Delta\phi Total')
title([{'Sensitivity of \Delta\phi Measured'},{'to Soil Dielectric State'}])
xlabel(['\Delta','SWE = 25 mm']);
% xlh = xlabel([{'Soil State Varies \Delta \phi Total'},{'(\Delta SWE = 30 mm)'}]);
ylabel('\Delta\phi (rad)')
xticks([1 2 3])
% xticklabels({" ", " ", " " })
xticklabels({'Frozen','Moist','Wet'});
yticks([-pi,-pi/2,0,pi/2,pi])
% yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
yticklabels({'', '','0','\pi/2','\pi'})

ylim([-pi pi])
xlim([0,4])
yl = ylim;
% xlh.Position(2) = yl(1)-0.125;
% legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.65)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; %grid minor;
axis square
text(-0.1,1.0675,'(c)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');
% set(gcf,'Position',[226.600000000000	242	1137.60000000000	420])
set(gcf,'Position',[226.600000000000	242	1116.80000000000	420])


% --- In Panel C, after making the main axes ---
axC = gca;

% Place inset axes (tune position numbers)
axInset1 = axes('Position',[0.69125 0.355 0.2135 0.1]); % eps'
axInset2 = axes('Position',[0.69125 0.225 0.2135 0.1]); % eps''

% Example values: eps_soil = [eps_frozen; eps_moist; eps_wet] (complex)
% eps' inset
axes(axInset1);
hold on;
% bar(0.7,eps_r(1),1,'FaceColor',[0,0,0]);
% bar(2,eps_r(2),1,'FaceColor',[0.35,0.35,0.35]);
% bar(3.3,eps_r(3),1,'FaceColor',[0.65,0.65,0.65]);
bar(0.825,eps_r(1),0.75,'FaceColor',[0,0,0]);
bar(2,eps_r(2),0.75,'FaceColor',[0.35,0.35,0.35]);
bar(3.175,eps_r(3),0.75,'FaceColor',[0.65,0.65,0.65]);

ylim([0,15])
% xticks([]); box on;
% xticks([0.7,2,3.3]); box on;
xticks([0.825,2,3.175]); box on;

xticklabels({" ", " ", " " })
% xticks([1 2 3]);
ylabel('\epsilon^\prime','FontSize',9);
set(gca,'FontSize',9,'FontWeight','bold','FontName','Serif');
grid on;
% text(.2,1.25,'Soil Complex Permittivity','Units','normalized','FontWeight','bold','FontSize',10,'FontName','Serif');
text(.1,1.25,'Effective Soil Permittivity (\epsilon^\prime, \epsilon^{\prime\prime})','Units','normalized','FontWeight','bold','FontSize',10,'FontName','Serif');



% eps'' inset
axes(axInset2);
% bar(eps_i,0.75);
hold on;
% bar(0.7,eps_i(1),1,'FaceColor',[0,0,0]);
% bar(2,eps_i(2),1,'FaceColor',[0.35,0.35,0.35]);
% bar(3.3,eps_i(3),1,'FaceColor',[0.65,0.65,0.65]);
bar(0.825,eps_i(1),0.75,'FaceColor',[0,0,0]);
bar(2,eps_i(2),0.75,'FaceColor',[0.35,0.35,0.35]);
bar(3.175,eps_i(3),0.75,'FaceColor',[0.65,0.65,0.65]);
% xticks([]); box on;
% xticks([0.7,2,3.3]); box on;
xticks([0.825,2,3.175]); box on;

xticklabels({" ", " ", " " })
% xticks(1:3); xticklabels({'F','M','W'}); box on;
ylabel('\epsilon^{\prime\prime}','FontSize',9);
set(gca,'FontSize',9,'FontWeight','bold','FontName','Serif');
grid on;
% Save this beauty
% exportgraphics(gcf,['E:\WSC24\AbstrctFigure2.png'],Resolution=300)
%% 0 -100 mm SWE Change with Soil Change
hs2 = linspace(1,4./3,100);
test_dphi_det_frz = zeros(2,2,numel(hs2));
test_dphi_det_moist = zeros(2,2,numel(hs2));
test_dphi_det_wet = zeros(2,2,numel(hs2));
deltaSWEsoilHs = zeros(numel(hs2),1);
% test_dphi_det_snow = zeros(2,2,3);
% test_dphi_det_soil = zeros(2,2,3);
% test_dphi_Ts_snow = zeros(3,1);
% test_dphi_soil_res = dphi_det_soil;
% test_dphi_det_air = dphi_det_soil;
% Recall Pars
for kk = 1:3
    for jj = 1:numel(hs2)
    pars1 = pars; pars2 = pars;

% ---------------- Snow definitions ----------------
snow1 = struct('Hs',1,'rho',300,'lwc',lwc(1),'Tk',TempK);
snow2 = struct('Hs',hs2(jj),'rho',300,'lwc',lwc(1),'Tk',TempK);  % wetter pass
deltaSWEsoilHs(jj) = snow2.Hs.*snow2.rho-snow1.Hs.*snow1.rho;

% Snow Pars for run
pars1.snowState = snow1;
pars2.snowState = snow2;

soilprof = struct();
soilprof.enable = true;
soilprof.zmax = 0.50;
soilprof.dz   = 0.001;
soilprof.type   = 'exp';
soilprof.SWC = 0.55;
soilprof.sigma_bulk = 0.03;
soilprof.z0 = 0.07;

soilprof.Tc_profile = struct('type','exp','Tc0',-2.0,'Tc_inf',2.0,'zT',0.10);

soilprof.soilParams = soilParams;
soilprof.wparams    = wparams;
soilprof.tauBW      = tauBW;
soilprof.clayParams = clayParams;
soilprof.freeze = struct('enable',false); % default unfrozen
soilprof.freeze.uwf = struct('mode','texture','fmin',0.02);

switch kk
  case 1 % Frozen
    soilprof.VWC0 = 0.20; soilprof.VWCinf = 0.10;
    soilprof.freeze.enable = true;
    soilprof.freeze.mode = 'physical';
    soilprof.freeze.T0 = 0.0; soilprof.freeze.dT = 0.5;
    soilprof.freeze.p_sigma = 2;
    soilprof.freeze.eta_sigma = 2;
    soilprof.freeze.sigma_floor = 0.0;

  case 2 % Moist (unfrozen)
    soilprof.VWC0 = 0.20; soilprof.VWCinf = 0.10;
    soilprof.freeze.enable = false;

  case 3 % Wet (unfrozen, higher VWC)
    soilprof.VWC0 = 0.30; soilprof.VWCinf = 0.15;
    soilprof.freeze.enable = false;
end

pars1.soil_profile = soilprof;
pars2.soil_profile = soilprof;  % identical soil within this kk
assert(isequaln(pars1.soil_profile, pars2.soil_profile), 'Soil differs between passes!');

% --- Deterministic (mean-field) outputs ---
o2 = run_pass(theta_i, f, snow1, eps_g0, pars1);
o1 = run_pass(theta_i, f, snow2, eps_g0, pars2);

% Compute Enery Weighted Effective Complex Permittivity
switch kk
    case 1 % Frozen
        eps_eff_frozen = eps_eff_energy(o1.profile.diag);
        E1 = o1.E;   E2 = o2.E;
        test_dphi_det_frz(:,:,jj)       = angle(E1  .* conj(E2));
    case 2 % Moist (unfrozen)
        eps_eff_moist = eps_eff_energy(o1.profile.diag);
        E1 = o1.E;   E2 = o2.E;
        test_dphi_det_moist(:,:,jj)       = angle(E1  .* conj(E2));

    case 3 % Wet (unfrozen, higher VWC)
        eps_eff_wet= eps_eff_energy(o1.profile.diag);
        E1 = o1.E;   E2 = o2.E;
        test_dphi_det_wet(:,:,jj)       = angle(E1  .* conj(E2));
end
% 
% dphi_det(kk) = angle(ensure_2x2(o1.E) .* conj(ensure_2x2(o2.E)));
% dphi_det_snow(kk) = angle(ensure_2x2(o1.Es) .* conj(ensure_2x2(o2.Es)));
% dphi_det_soil(kk) = angle(ensure_2x2(o1.Eg) .* conj(ensure_2x2(o2.Eg)));

% dphi_det(:,:,kk) = angle((o1.E) .* conj((o2.E)));
% dphi_det_snow(:,:,kk) = angle((o1.Es) .* conj((o2.Es)));
% dphi_Ts_snow(kk) = angle(o1.Ts.*conj(o2.Ts));
% dphi_det_soil(:,:,kk) = angle((o1.Eg) .* conj((o2.Eg)));

% E1 = o1.E;   E2 = o2.E;
% Es1 = o1.Es; Es2 = o2.Es;
% Eas1 = o1.Eas; Eas2 = o2.Eas;
% 
% % NEW: ground/soil total in 4-component model
% Eg1 = o1.Esg + o1.Egv;
% Eg2 = o2.Esg + o2.Egv;
% 
% Eg1_res = Eg1 ./ o1.Ts_eff;
% Eg2_res = Eg2 ./ o2.Ts_eff;
% 
% dphi_soil_res = angle( Eg1_res .* conj(Eg2_res) );




% test_dphi_det_air(:,:,kk)  = angle(Eas1 .* conj(Eas2));
% test_dphi_det_snow(:,:,kk)  = angle(Es1 .* conj(Es2));
% test_dphi_det_soil(:,:,kk)  = angle(Eg1 .* conj(Eg2));
% 
% test_dphi_Ts_snow(kk) = angle(o1.Ts_excess .* conj(o2.Ts_excess));
    end
end

soilPhzFrz = squeeze(test_dphi_det_frz(2,2,:));
soilPhzMoist = squeeze(test_dphi_det_moist(2,2,:));
soilPhzWet = squeeze(test_dphi_det_wet(2,2,:));

eps_eff = [eps_eff_frozen; eps_eff_moist; eps_eff_wet];

eps_r = real(eps_eff);
eps_i = imag(eps_eff);

%% Figure!
figure(); subplot(1,3,1)
box on;
set(gcf,'Position',[226.600000000000	242	1116.80000000000	420])
hold on;
plot(oviDSWE,oviDPhi,'k','linewidth', 2,"DisplayName","Refraction Limit")
plot(deltaSWE,dphi_Ts_snow,'linewidth', 2, 'color',[0.5,0.5,0.5],"DisplayName","Snow Propagation")
% plot(deltaSWE,phzSoilTs,'--r',"DisplayName","Ground Residual \phi")
% title('Canonical Refraction Limit')
title('Refraction-Based Phase Limit')
xlabel('\DeltaSWE (mm)')
ylabel('\Delta\phi (rad)')
xticks([0,25,50,75,100])
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.75)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; %grid minor;
axis square
text(-0.1,1.0675,'(a)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');
subplot(1,3,2)
box on;
hold on;
plot(deltaSWE,phz,'k','linewidth', 2, "DisplayName","Modeled Total");
plot(deltaSWE,dphi_Ts_snow,'linewidth', 2, 'color',[0.5,0.5,0.5],"DisplayName","Snow Propagation");
% plot(deltaSWE,phzAir,'linewidth', 1.5, 'color',[0.25,0.25,0.25],'lineStyle',"--","DisplayName","Snow Interface \phi"); 
plot(deltaSWE,phzSnow,'linewidth', 1.25, 'color',[0.65,0.65,0.65],'lineStyle',"--","DisplayName","Snow Volume "); 
% plot(oviDSWE,oviDPhi,'r',"DisplayName","Oveisgharan \phi")
% plot(deltaSWE,phzSoilTs)
title('Coherent Phasor Mixing')
% title('Complex Phasor Summation')
xlabel('\DeltaSWE (mm)')
ylabel('\Delta\phi (rad)')
xticks([0,25,50,75,100])
yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
ylim([-pi pi])
xlim([0,100])
legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.75)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on;%grid minor;
axis square
text(-0.1,1.0675,'(b)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');


% subplot(1,3,3)
% Top axis (phase curves)
% posSWE = [0.69125 0.5525  0.2135 0.25];  % adjust to your layout
posSWE = [0.69125 0.5881  0.2135 0.2145];  % adjust to your layout

axSWE = axes('Position',[posSWE]);
axes(axSWE);
box on;

hold on;
plot(deltaSWEsoilHs,soilPhzFrz,"DisplayName",'Frozen','LineWidth',2,'Color',[0,0,0]);
plot(deltaSWEsoilHs,soilPhzMoist,"DisplayName","Moist",'LineWidth',2,'Color',[0.35,0.35,0.35]);
plot(deltaSWEsoilHs,soilPhzWet,'DisplayName','Wet','LineWidth',2,'Color',[0.65,0.65,0.65])
ylim([-pi pi])
yticks([-pi -pi/2 0 pi/2 pi])
xticks([0,25,50,75,100])

% xticklabels([])  % hide top x labels
% title('Soil State Varies \Delta\phi Total')
% title([{'Sensitivity of \Delta\phi Modeled'},{'to Soil Dielectric State'}])
set(gca,'FontSize',10,'FontWeight','bold','FontName','Serif');
title([{'Soil State Modulates \Delta\phi'}],'FontSize',13)
xlabel(['\Delta','SWE (mm)'],"FontSize",12);

% xlh = xlabel([{'Soil State Varies \Delta \phi Total'},{'(\Delta SWE = 30 mm)'}]);
ylabel('\Delta\phi (rad)','FontSize',12)
% xticks([1 2 3])
% % xticklabels({" ", " ", " " })
% xticklabels({'Frozen','Moist','Wet'});
% yticks([-pi,-pi/2,0,pi/2,pi])
yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
% yticklabels({'', '','0','\pi/2','\pi'})

ylim([-pi pi])
xlim([0,100])
yl = ylim;
% xlh.Position(2) = yl(1)-0.125;
% legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.65)
% set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; %grid minor;
% axis square
% legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.75)

text(-0.1,1.16,'(c)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');
% set(gcf,'Position',[226.600000000000	242	1137.60000000000	420])
set(gcf,'Position',[226.600000000000	242	1116.80000000000	420])


% --- In Panel C, after making the main axes ---
% axC = gca;

% Place inset axes (tune position numbers)
axInset1 = axes('Position',[0.69125 0.36 0.2135 0.1]); % eps'
axInset2 = axes('Position',[0.69125 0.225 0.2135 0.1]); % eps''

% Example values: eps_soil = [eps_frozen; eps_moist; eps_wet] (complex)
% eps' inset
axes(axInset1);
hold on;
% bar(0.7,eps_r(1),1,'FaceColor',[0,0,0]);
% bar(2,eps_r(2),1,'FaceColor',[0.35,0.35,0.35]);
% bar(3.3,eps_r(3),1,'FaceColor',[0.65,0.65,0.65]);
bar(0.825,eps_r(1),0.75,'FaceColor',[0,0,0]);
bar(2,eps_r(2),0.75,'FaceColor',[0.35,0.35,0.35]);
bar(3.175,eps_r(3),0.75,'FaceColor',[0.65,0.65,0.65]);

ylim([0,15])
% xticks([]); box on;
% xticks([0.7,2,3.3]); box on;
xticks([0.825,2,3.175]); box on;

xticklabels({" ", " ", " " })
% xticks([1 2 3]);
set(gca,'FontSize',10,'FontWeight','bold','FontName','Serif');
ylabel('\epsilon^\prime','FontSize',12);

grid on;
% text(.2,1.25,'Soil Complex Permittivity','Units','normalized','FontWeight','bold','FontSize',10,'FontName','Serif');
% text(.1,1.25,'Effective Soil Permittivity (\epsilon^\prime, \epsilon^{\prime\prime})','Units','normalized','FontWeight','bold','FontSize',10,'FontName','Serif');



% eps'' inset
axes(axInset2);
% bar(eps_i,0.75);
hold on;
% bar(0.7,eps_i(1),1,'FaceColor',[0,0,0]);
% bar(2,eps_i(2),1,'FaceColor',[0.35,0.35,0.35]);
% bar(3.3,eps_i(3),1,'FaceColor',[0.65,0.65,0.65]);
bar(0.825,eps_i(1),0.75,'FaceColor',[0,0,0]);
bar(2,eps_i(2),0.75,'FaceColor',[0.35,0.35,0.35]);
bar(3.175,eps_i(3),0.75,'FaceColor',[0.65,0.65,0.65]);
% xticks([]); box on;
% xticks([0.7,2,3.3]); box on;
xticks([0.825,2,3.175]); box on;
% ylim([0 5])
% xticklabels({" ", " ", " " })
% xticks(1:3); xticklabels({'F','M','W'}); box on;
set(gca,'FontSize',10,'FontWeight','bold','FontName','Serif');
ylabel('\epsilon^{\prime\prime}','FontSize',12);
xticklabels({'Frozen','Moist','Wet'});
axInset2.XAxis.FontSize = 12;
xlabel('Effective Soil Permittivity','FontSize',13)

grid on;
% Save this beauty
% exportgraphics(gcf,['E:\WSC24\AbstrctFigureTomO.png'],Resolution=300)

%% MCS Parameterization 02/23/26 - 02/26/26
test_dphi_det = zeros(2,2,3);
test_dphi_det_snow = zeros(2,2,3);
test_dphi_det_soil = zeros(2,2,3);
test_dphi_Ts_snow = zeros(3,1);
test_dphi_soil_res = zeros(2,2,numel(1,1));
test_dphi_det_air = zeros(2,2,numel(1,1));
% Recall Pars
for kk = 1
    pars1 = pars; pars2 = pars;

% ---------------- Snow definitions ----------------
% snow1 = struct('Hs',1,'rho',300,'lwc',lwc(1),'Tk',TempK);
% snow2 = struct('Hs',1.0833,'rho',300,'lwc',lwc(1),'Tk',TempK);  % wetter pass
TempK = 273.15;
snow1 = struct('Hs',0.862,'rho',311.5,'lwc',0.0025,'Tk',TempK);
snow2 = struct('Hs',0.941,'rho',312.5,'lwc',0.015,'Tk',TempK);  % wetter pass
deltaSWEsoil = snow2.Hs.*snow2.rho-snow1.Hs.*snow1.rho;

% Snow Pars for run
pars1.snowState = snow1;
pars2.snowState = snow2;

soilprof = struct();
soilprof.enable = true;
soilprof.zmax = 0.50;
soilprof.dz   = 0.001;
soilprof.type   = 'exp';
soilprof.SWC = 0.55;
soilprof.sigma_bulk = 0.03;
soilprof.z0 = 0.07;

soilprof.Tc_profile = struct('type','exp','Tc0',-.25,'Tc_inf',-1.5,'zT',.05);

soilprof.soilParams = soilParams;
soilprof.wparams    = wparams;
soilprof.tauBW      = tauBW;
soilprof.clayParams = clayParams;
soilprof.freeze = struct('enable',false); % default unfrozen
% soilprof.freeze.uwf = struct('mode','texture','fmin',0.25);
soilprof.freeze.uwf = struct('mode','logistic','T0',-0.35,'dT',0.2, 'fmin',0.02);

switch kk
  case 1 % Frozen
    soilprof.VWC0 = 0.45; soilprof.VWCinf = 0.30;
    soilprof.freeze.enable = true;
    soilprof.freeze.mode = 'physical';
    soilprof.freeze.T0 = -1.0; soilprof.freeze.dT = 0.25;
    soilprof.freeze.p_sigma = 2;
    soilprof.freeze.eta_sigma = 2;
    soilprof.freeze.sigma_floor = 0.0;

  case 2 % Moist (unfrozen)
    soilprof.VWC0 = 0.20; soilprof.VWCinf = 0.10;
    soilprof.freeze.enable = false;

  case 3 % Wet (unfrozen, higher VWC)
    soilprof.VWC0 = 0.30; soilprof.VWCinf = 0.15;
    soilprof.freeze.enable = false;
end

pars1.soil_profile = soilprof;
pars2.soil_profile = soilprof;  % identical soil within this kk
assert(isequaln(pars1.soil_profile, pars2.soil_profile), 'Soil differs between passes!');

% --- Deterministic (mean-field) outputs ---
o1 = run_pass(theta_i, f, snow1, eps_g0, pars1);
o2 = run_pass(theta_i, f, snow2, eps_g0, pars2);

% Compute Enery Weighted Effective Complex Permittivity
switch kk
    case 1 % Frozen
        eps_eff_frozen = eps_eff_energy(o1.profile.diag);
    case 2 % Moist (unfrozen)
        eps_eff_moist = eps_eff_energy(o1.profile.diag);
    case 3 % Wet (unfrozen, higher VWC)
        eps_eff_wet= eps_eff_energy(o1.profile.diag);
end
% 
% dphi_det(kk) = angle(ensure_2x2(o1.E) .* conj(ensure_2x2(o2.E)));
% dphi_det_snow(kk) = angle(ensure_2x2(o1.Es) .* conj(ensure_2x2(o2.Es)));
% dphi_det_soil(kk) = angle(ensure_2x2(o1.Eg) .* conj(ensure_2x2(o2.Eg)));

% dphi_det(:,:,kk) = angle((o1.E) .* conj((o2.E)));
% dphi_det_snow(:,:,kk) = angle((o1.Es) .* conj((o2.Es)));
% dphi_Ts_snow(kk) = angle(o1.Ts.*conj(o2.Ts));
% dphi_det_soil(:,:,kk) = angle((o1.Eg) .* conj((o2.Eg)));

E1 = o1.E;   E2 = o2.E;
Es1 = o1.Es; Es2 = o2.Es;
Eas1 = o1.Eas; Eas2 = o2.Eas;

% NEW: ground/soil total in 4-component model
Eg1 = o1.Esg + o1.Egv;
Eg2 = o2.Esg + o2.Egv;

Eg1_res = Eg1 ./ o1.Ts_eff;
Eg2_res = Eg2 ./ o2.Ts_eff;

dphi_soil_res = angle( Eg1_res .* conj(Eg2_res) );

test_dphi_det(:,:,kk)       = angle(E1  .* conj(E2));
test_dphi_det_air(:,:,kk)  = angle(Eas1 .* conj(Eas2));
test_dphi_det_snow(:,:,kk)  = angle(Es1 .* conj(Es2));
test_dphi_det_soil(:,:,kk)  = angle(Eg1 .* conj(Eg2));

test_dphi_Ts_snow(kk) = angle(o1.Ts_excess .* conj(o2.Ts_excess));
end

soilPhz = squeeze(test_dphi_det(2,2,:));
% eps_eff = [eps_eff_frozen; eps_eff_moist; eps_eff_wet];
eps_eff = [eps_eff_frozen; eps_eff_frozen; eps_eff_frozen];


eps_r = real(eps_eff);
eps_i = imag(eps_eff);
figure();
subplot(1,3,3)
hold on;
yline(0,'k-','LineWidth',1);
bar(1,soilPhz(1),0.8, 'FaceColor',[0,0,0],"DisplayName","Frozen");
bar(2,soilPhz(2), 0.8,'FaceColor',[0.35,0.35,0.35],"DisplayName","Moist");
bar(3,soilPhz(3),0.8,'FaceColor',[0.65,0.65,0.65], "DisplayName","Wet");

% plot(deltaSWE,phzAir,'linewidth', 1.5, 'color',[0.25,0.25,0.25],'lineStyle',"--","DisplayName","Snow Interface \phi"); 
% plot(deltaSWE,phzSnow,'linewidth', 1.5, 'color',[0.5,0.5,0.5],'lineStyle',"--","DisplayName","\Delta \phi Snow Volume "); 
% plot(oviDSWE,oviDPhi,'r',"DisplayName","Oveisgharan \phi")
% plot(deltaSWE,phzSoilTs)
% title('Soil State Varies \Delta\phi Total')
title([{'Sensitivity of \Delta\phi Measured'},{'to Soil Dielectric State'}])
xlabel(['\Delta','SWE = 25 mm']);
% xlh = xlabel([{'Soil State Varies \Delta \phi Total'},{'(\Delta SWE = 30 mm)'}]);
ylabel('\Delta\phi (rad)')
xticks([1 2 3])
% xticklabels({" ", " ", " " })
xticklabels({'Frozen','Moist','Wet'});
yticks([-pi,-pi/2,0,pi/2,pi])
% yticklabels({'-\pi', '-\pi/2','0','\pi/2','\pi'})
yticklabels({'', '','0','\pi/2','\pi'})

ylim([-pi pi])
xlim([0,4])
yl = ylim;
% xlh.Position(2) = yl(1)-0.125;
% legend('Location','southwest','FontSize',10,'BackgroundAlpha',0.65)
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)
grid on; %grid minor;
axis square
text(-0.1,1.0675,'(c)','Units','normalized','FontWeight','bold','FontSize',14,'FontName','Serif');
% set(gcf,'Position',[226.600000000000	242	1137.60000000000	420])
set(gcf,'Position',[226.600000000000	242	1116.80000000000	420])


% --- In Panel C, after making the main axes ---
axC = gca;

% Place inset axes (tune position numbers)
axInset1 = axes('Position',[0.69125 0.355 0.2135 0.1]); % eps'
axInset2 = axes('Position',[0.69125 0.225 0.2135 0.1]); % eps''

% Example values: eps_soil = [eps_frozen; eps_moist; eps_wet] (complex)
% eps' inset
axes(axInset1);
hold on;
% bar(0.7,eps_r(1),1,'FaceColor',[0,0,0]);
% bar(2,eps_r(2),1,'FaceColor',[0.35,0.35,0.35]);
% bar(3.3,eps_r(3),1,'FaceColor',[0.65,0.65,0.65]);
bar(0.825,eps_r(1),0.75,'FaceColor',[0,0,0]);
bar(2,eps_r(2),0.75,'FaceColor',[0.35,0.35,0.35]);
bar(3.175,eps_r(3),0.75,'FaceColor',[0.65,0.65,0.65]);

ylim([0,15])
% xticks([]); box on;
% xticks([0.7,2,3.3]); box on;
xticks([0.825,2,3.175]); box on;

xticklabels({" ", " ", " " })
% xticks([1 2 3]);
ylabel('\epsilon^\prime','FontSize',9);
set(gca,'FontSize',9,'FontWeight','bold','FontName','Serif');
grid on;
% text(.2,1.25,'Soil Complex Permittivity','Units','normalized','FontWeight','bold','FontSize',10,'FontName','Serif');
text(.1,1.25,'Effective Soil Permittivity (\epsilon^\prime, \epsilon^{\prime\prime})','Units','normalized','FontWeight','bold','FontSize',10,'FontName','Serif');



% eps'' inset
axes(axInset2);
% bar(eps_i,0.75);
hold on;
% bar(0.7,eps_i(1),1,'FaceColor',[0,0,0]);
% bar(2,eps_i(2),1,'FaceColor',[0.35,0.35,0.35]);
% bar(3.3,eps_i(3),1,'FaceColor',[0.65,0.65,0.65]);
bar(0.825,eps_i(1),0.75,'FaceColor',[0,0,0]);
bar(2,eps_i(2),0.75,'FaceColor',[0.35,0.35,0.35]);
bar(3.175,eps_i(3),0.75,'FaceColor',[0.65,0.65,0.65]);
% xticks([]); box on;
% xticks([0.7,2,3.3]); box on;
xticks([0.825,2,3.175]); box on;

xticklabels({" ", " ", " " })
% xticks(1:3); xticklabels({'F','M','W'}); box on;
ylabel('\epsilon^{\prime\prime}','FontSize',9);
set(gca,'FontSize',9,'FontWeight','bold','FontName','Serif');
grid on;
% Save this beauty
% exportgraphics(gcf,['E:\WSC24\AbstrctFigure2.png'],Resolution=300)

%% Sweep of lwc
lwc1_vec = linspace(0.0, 0.02, 401);   % 0 to 2%
theta_i = 60;
out = sweep_initial_lwc_phase(theta_i, f, snow1, snow2, eps_g0, pars1, pars2, lwc1_vec, 'VV');

figure('Color','w'); hold on;
plot(100*out.lwc1, out.ph.tot, 'k', 'LineWidth', 1.5, 'DisplayName', 'Measured \phi (wrapped)');
plot(100*out.lwc1, out.ph.g, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5, 'DisplayName', 'Ground \phi');
plot(100*out.lwc1, out.ph.ts, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5, 'DisplayName', 'T_s \phi');

yl = ylim;
for x = out.wrap.lwc_tot(:).'
    xline(100*x, ':k', 'LineWidth', 1.0);
end

xlabel('Initial LWC in pass 1 (%)');
ylabel('\Delta \phi (rad)');
ylim([-pi pi]);
grid on; grid minor;
legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
title('Wrapped phase vs. initial LWC');

figure('Color','w'); hold on;
plot(100*out.lwc1, out.phu.tot, 'k', 'LineWidth', 1.8, 'DisplayName', 'Measured \phi');
plot(100*out.lwc1, out.phu.g, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5, 'DisplayName', 'Ground \phi');
plot(100*out.lwc1, out.phu.ts, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5, 'DisplayName', 'T_s \phi');
plot(100*out.lwc1, out.phu.as, ':', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'DisplayName', 'Air-snow \phi');
plot(100*out.lwc1, out.phu.sv, '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'DisplayName', 'Snow volume \phi');

xlabel('Initial LWC in pass 1 (%)');
ylabel('Unwrapped \Delta \phi (rad)');
grid on; grid minor;
legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
title('Unwrapped phase vs. initial LWC');

figure('Color','w'); hold on;
plot(100*out.lwc1, out.mag.E1, 'k', 'LineWidth', 1.5, 'DisplayName', '|E_1| total');
plot(100*out.lwc1, out.mag.Eg1, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.5, 'DisplayName', '|E_g|');
plot(100*out.lwc1, out.mag.Eas1, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5, 'DisplayName', '|E_{as}|');
plot(100*out.lwc1, out.mag.Es1, '-.', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'DisplayName', '|E_s|');

xline(100*out.minmag.lwc_E1, ':k', 'LineWidth', 1.0, 'DisplayName', 'min |E_1|');

xlabel('Initial LWC in pass 1 (%)');
ylabel('Magnitude');
grid on; grid minor;
legend('Location','best');
set(gca,'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
title('Magnitude balance vs. initial LWC');

%% Sweep lwc and incidence angle

theta_vec = 35:1:75;
lwc1_vec  = 0:0.00025:0.02;   % 0 to 2%

snow1 = struct('Hs',0.862,'rho',311.5,'lwc',0.0025,'Tk',273.15); % early template
snow2 = struct('Hs',0.941,'rho',312.5,'lwc',0.0150,'Tk',273.15); % late fixed

% Example observed referenced phase to compare against
obs_phase_ref = 1.3415;   % replace with representative measured value if desired

out2D = sweep_lwc_theta_phase(theta_vec, lwc1_vec, f, snow1, snow2, ...
    eps_g0, pars1, pars2, obs_phase_ref, 'VV');

figure('Color','w');
imagesc(100*lwc1_vec, theta_vec, out2D.ph_wrap.tot);
set(gca,'YDir','normal');
xlabel('Initial LWC in pass 1 (%)');
ylabel('Incidence angle (deg)');
title('Wrapped modeled phase, \angle(E_{late}E_{early}^*)');
colorbar;
colormap(parula);
caxis([-pi pi]);

figure('Color','w');
imagesc(100*lwc1_vec, theta_vec, out2D.ph_unwrap.tot);
set(gca,'YDir','normal');
xlabel('Initial LWC in pass 1 (%)');
ylabel('Incidence angle (deg)');
title('Unwrapped modeled phase');
colorbar;
colormap(parula);

figure('Color','w');
imagesc(100*lwc1_vec, theta_vec, out2D.misfit_wrap);
set(gca,'YDir','normal');
xlabel('Initial LWC in pass 1 (%)');
ylabel('Incidence angle (deg)');
title('Wrapped phase misfit to observed referenced phase');
colorbar;
colormap(parula);
caxis([-pi pi]);

figure('Color','w');
imagesc(100*lwc1_vec, theta_vec, out2D.misfit_wrap);
set(gca,'YDir','normal'); hold on;
contour(100*lwc1_vec, theta_vec, out2D.misfit_wrap, [0 0], 'k', 'LineWidth', 2);
xlabel('Initial LWC in pass 1 (%)');
ylabel('Incidence angle (deg)');
title('Wrapped phase misfit with zero-misfit contour');
colorbar;
caxis([-pi pi]);

figure('Color','w');
imagesc(100*lwc1_vec, theta_vec, out2D.wrap_index);
set(gca,'YDir','normal');
xlabel('Initial LWC in pass 1 (%)');
ylabel('Incidence angle (deg)');
title('Wrap index of modeled phase');
colorbar;

%% Let's synthesize the MCS interfeogram

%% -------- helpers --------
function X2 = ensure_2x2(X)
if isempty(X), X2 = complex(zeros(2,2));
elseif isscalar(X), X2 = [X 0; 0 X];
else, X2 = X;
end
end

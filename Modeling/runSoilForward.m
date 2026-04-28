% Example inputs from HydraGO/site
theta_corr = d1cors;
VWC      = theta_corr./100;                    % m3/m3
SWC      = 0.5171;                          % porosity (best guess)
EC_dSm   = hp1sigmatc;                     % dS/m from HydraGO
Tc       = T1.HP_1_SoilTemperature_C_;                    % deg C
f_L      = 1.3e9;                         % L-band Hz

% Water parameters (Kaatze-polynomial based)
epssCoefs = [88.119019607843190,-0.484234778121778,0.004464654282766,-2.393188854489220e-05]';
tauCoefs = [1.736936274509802e-11,-4.823117475060188e-13,5.488493292053644e-15,-2.366013071895416e-17]';
wparams.epssCoefs = epssCoefs;   % your 4-coef vec
wparams.tauCoefs  = tauCoefs;
wparams.eps_inf   = 4.8;
wparams.alpha     = 0.001;

% Simple soil defaults, or pass a struct
soilParams.eps_real_soil = 4.0;
soilParams.eps_real_clay = 25;
soilParams.frac_clay     = 0.1;

outL = soil_eps_from_hydra_forward( ...
    VWC, SWC, EC_dSm, Tc, f_L, ...
    soilParams, wparams, []);   % [] → default tauBW
%% Multifrequency
% HydraGo Weighted Inversion
freq_hydra = 50e6;
eps_hydra = 12.2703 + 1i.*10.6567;
weight_hydra = 2.5;
sigma_hydra = 0.02933;

VWC    = 0.202666666666667;      % N×1
SWC    = 0.556753;               % scalar porosity
EC_dSm = 0.2933;       % N×1, dS/m
Tc     = 17.33;      % N×1
freqs  = [50e6, 1.3e9, 3.2e9];   % [50 MHz, L, S]
nF     = numel(freqs);
nN     = numel(VWC);

% Preallocate arrays: [nN x nF]
eps_complex_all     = complex(zeros(nN, nF));
eps_real_all        = zeros(nN, nF);
eps_im_total_all    = zeros(nN, nF);
eps_im_cond_all     = zeros(nN, nF);
eps_im_dip_all      = zeros(nN, nF);
sigma_total_all     = zeros(nN, nF);
sigma_cond_all      = zeros(nN, nF);
sigma_dip_all       = zeros(nN, nF);

for k = 1:nF
    f = freqs(k);

    out = soil_eps_from_hydra_forward( ...
        VWC, SWC, EC_dSm, Tc, f, ...
        soilParams, wparams, tauBW);

    eps_complex_all(:,k)  = out.eps_complex;
    eps_real_all(:,k)     = out.eps_real;
    eps_im_total_all(:,k) = out.eps_im_total;
    eps_im_cond_all(:,k)  = out.eps_im_cond;
    eps_im_dip_all(:,k)   = out.eps_im_dip;

    sigma_total_all(:,k)  = out.sigma_eff_total;
    sigma_cond_all(:,k)   = out.sigma_eff_cond;
    sigma_dip_all(:,k)    = out.sigma_eff_dip;
end

%% Compact loop version
freqs = [50e6, 1.3e9, 3.2e9];  % 50 MHz, L-band, S-band
results = soil_spectrum(VWC, SWC, EC_dSm, Tc, freqs, soilParams, wparams, []);
%% Forward Modeling Inversion
% Load Complex S11 
dataDir = 'E:\WCP\0925campaign\fieldfox\data';
% Load Calibration Data
airReal = readtable([dataDir,'\S11_AIR_REAL.csv']);
airImag = readtable([dataDir,'\S11_AIR_IMAG.csv']);
shortReal = readtable([dataDir,'\S11_SHORT_REAL.csv']);
shortImag = readtable([dataDir,'\S11_SHORT_IMAG.csv']);
waterReal = readtable([dataDir,'\S11_WATER_REAL.csv']);
waterImag = readtable([dataDir,'\S11_WATER_IMAG.csv']);
% Additional Calibration files for Short Synthesis
airRealCal = readtable([dataDir,'\S11_OPEN_CAL_AIR_REAL.csv']);
airImagCal = readtable([dataDir,'\S11_OPEN_CAL_AIR_IMAG.csv']);
openRealCal = readtable([dataDir,'\S11_OPEN_CAL_REAL.csv']);
openImagCal = readtable([dataDir,'\S11_OPEN_CAL_IMAG.csv']);
shortRealCal = readtable([dataDir,'\S11_SHORT_CAL_REAL.csv']);
shortImagCal = readtable([dataDir,'\S11_SHORT_CAL_IMAG.csv']);


% Load Soil Measurements
wsReal = readtable([dataDir,'\S11_WS4_REAL.csv']);
wsImag = readtable([dataDir,'\S11_WS4_IMAG.csv']);

% Form Complex Number
freq_vna = airReal.x_DATAFreq; % in Hz
s_air    = airReal.S11 - 1j*airImag.S11;
s_short  = shortReal.S11 - 1j*shortImag.S11;
s_water  = waterReal.S11 - 1j*waterImag.S11;
s_sample = wsReal.S11 - 1j*wsImag.S11;
% Synthesize Short
s_air_cal = airRealCal.S11 - 1j*airImagCal.S11;
s_open_cal = openRealCal.S11 - 1j*openImagCal.S11;
s_short_cal = shortRealCal.S11 - 1j*shortImagCal.S11;


s_air = smooth_permittivity(freq_vna, s_air, 'lowess', struct('span', 0.175), false);
s_short = smooth_permittivity(freq_vna, s_short, 'lowess', struct('span', 0.175), false);
s_water = smooth_permittivity(freq_vna, s_water, 'lowess', struct('span', 0.175), false);
s_sample = smooth_permittivity(freq_vna, s_sample, 'lowess', struct('span', 0.175), false);

% S11 Data Inversion
eps0 = 8.854187817e-12; % vacuum permittivity [F/m]
% 3pt
% [eps_complex3, freq] = invert_stuchly3pt(s_air, s_short, s_water, s_sample, freq);
% 2Pt
% [eps_complex, freq] = invert_stuchly2pt(s_air, s_water, s_sample, freq);
s_std = 0.005;
[eps_vna, eps_std] = invert_stuchly2pt_uncertainty(s_air, s_water, s_sample, freq_vna, s_std);


% Remove NaN
nanIx = find(isnan(freq_vna)|isnan(eps_vna)|freq_vna<.5e9);
freq_vna(nanIx) = []; eps_vna(nanIx) = []; eps_std(nanIx) = [];
minEpsImag = min((imag(eps_vna)));

if minEpsImag < 0
eps_vna = real(eps_vna)+1i.*(imag(eps_vna)-(1.25.*minEpsImag));
end

% Pack data
VWC    = 0.202666666666667;      % N×1
SWC    = 0.556753;               % scalar porosity
EC_dSm = 0.2933;       % N×1, dS/m
Tc     = 17.33;      % N×1
freqs  = [50e6, 1.3e9, 3.2e9];   % [50 MHz, L, S]
nF     = numel(freqs);
nN     = numel(VWC);
data.VWC        = VWC;
data.SWC        = SWC;
data.EC_dSm     = EC_dSm;
data.Tc         = Tc;
data.freq_vna   = freq_vna;
data.eps_vna_meas = eps_vna;
data.freq_hydra = 50e6;
data.eps_hydra  = eps_hydra;
data.weight_hydra = 5;
data.tau_max    = 5e-10;
data.lambda_reg = 1e-2;

% Initial guesses for p
p0 = [4.0, 25.0, 0.15, 1.0, 0.01, 1.0];   % [eps_soil, eps_clay, frac_clay, A, B, tau_scale]

lb = [2.0, 10.0, 0.001, 0.5, 0.0, 0.5];
ub = [8.0, 40.0, 0.5 2.0, 0.1, 2.0];

costFcn = @(p) soil_forward_cost(p, data, wparams);

opts = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',3e3);
p_fit = fmincon(costFcn, p0, [], [], [], [], lb, ub, [], opts);

% Build fitted soilParams/tauBW
soilParams_fit.eps_real_soil = p_fit(1);
soilParams_fit.eps_real_clay = p_fit(2);
soilParams_fit.frac_clay     = p_fit(3);
tauBW_fit.A         = p_fit(4);
tauBW_fit.B         = p_fit(5);
tauBW_fit.tau_scale = p_fit(6);
tauBW_fit.max_tau   = data.tau_max;

%% Model Output
resultsInv = soil_spectrum(VWC, SWC, EC_dSm, Tc, freqs, soilParams_fit, wparams, tauBW_fit);
figure();plot(freq_vna,real(eps_vna),'k','LineWidth',2);
hold on; plot(cat(1,resultsInv.freq),real(cat(1,resultsInv.eps_complex)),'Color',[0.5,0.5,0.5],'linewidth',2)
plot(freq_vna,imag(eps_vna),'--k','LineWidth',2);
hold on; plot(cat(1,resultsInv.freq),imag(cat(1,resultsInv.eps_complex)),'Color',[0.5,0.5,0.5],'linewidth',2,'LineStyle','--')
grid on; grid minor
xlabel('Frequency (Hz)')
ylabel('Complex Permittivity')
legend('Real VNA','Real Model','Imag VNA', 'Imag Model')
title()

%% Take dual
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

clayParams.eps_inf_clay = 5;
clayParams.tau_clay     = 5e-10;
clayParams.alpha_clay   = 0.15;

freqs = [50e6, 1.3e9, 3.2e9];  % 50 MHz, L-band, S-band
results = soil_spectrum(VWC, SWC, EC_dSm, Tc, freqs, soilParams, wparams, [],[]);

% Test
[eps_complex, components] = soilPermittivity_advanced( ...
    VWC, SWC, EC_dSm.*0.1, soilParams, Tc, 1.3e9, wparams, tauBW, clayParams);
%% Dual Inversion
% --- 1. Load / define your data for one sample ---
% Pack data
VWC    = 0.202666666666667;      % N×1
SWC    = 0.556753;               % scalar porosity
EC_dSm = 0.2933;       % N×1, dS/m
Tc     = 17.33;      % N×1
freqs  = [50e6, 1.3e9, 3.2e9];   % [50 MHz, L, S]
nF     = numel(freqs);
nN     = numel(VWC);

data.VWC       = VWC;           % scalar or N×1 (use same VWC for lab sample)
data.SWC       = SWC;           % porosity
data.EC_dSm    = EC_dSm;        % bulk EC (dS/m)
data.Tc        = Tc;            % °C
data.freq_vna  = freq_vna;          % Nv×1 Hz
data.eps_vna_meas = eps_vna;   % Nv×1 complex

data.freq_hydra = 50e6;             % optional HydraGO anchor
data.eps_hydra  = eps_hydra;        % complex at 50 MHz

data.w_vna_real = 1.0;
data.w_vna_imag = 1.0;
data.w_hydra    = 5.0;              % emphasize HydraGO point

data.tau_max      = 5e-10;
data.clay_eps_inf = 5.0;
data.lambda_reg = 1e-2;


% --- 2. Water parameters (fixed parts from Kaatze) ---

wparams.epssCoefs   = epssCoefs;   % your 4×1 static perm poly
wparams.tauCoefs    = tauCoefs;    % your 4×1 tau poly
wparams.eps_inf     = 4.8;
wparams.alpha_free  = 0.01;
wparams.alpha_bound = 0.2;
% wparams.frac_bound will come from p(4)

% --- 3. Initial guess and bounds for p ---

p0 = [ ...
    4.0;      % eps_real_soil
    20.0;     % eps_real_clay
    0.20;     % frac_clay
    0.30;     % frac_bound
    5.0;      % bound_scale
    0.02;     % B_bound
    5e-10;    % tau_clay
    0.15];    % alpha_clay

lb = [ ...
    2.0;      % eps_real_soil
    8.0;      % eps_real_clay
    0.00;     % frac_clay
    0.00;     % frac_bound
    1.0;      % bound_scale
    0.0;      % B_bound
    5e-13;    % tau_clay
    0.0];     % alpha_clay

ub = [ ...
    8.0;      % eps_real_soil
    40.0;     % eps_real_clay
    0.50;     % frac_clay
    0.80;     % frac_bound
    30.0;     % bound_scale
    0.2;      % B_bound
    5e-9;     % tau_clay
    0.5];     % alpha_clay

% --- 4. Run lsqnonlin ---

resfun = @(p) soil_inversion_residual(p, data, wparams);

opts = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',4000, ...
    'MaxIterations',200);

[p_fit, resnorm, residual, exitflag, output] = lsqnonlin(resfun, p0, lb, ub, opts);

disp('Fitted parameter vector p_fit:');
disp(p_fit);
% Write to structure
soilParams_fit.eps_real_soil = p_fit(1);
soilParams_fit.eps_real_clay = p_fit(2);
soilParams_fit.frac_clay     = p_fit(3);

wparams_fit          = wparams;
wparams_fit.frac_bound = p_fit(4);

tauBW_fit.free_scale  = 1.0;
tauBW_fit.bound_scale = p_fit(5);
tauBW_fit.B_bound     = p_fit(6);
tauBW_fit.max_tau     = data.tau_max;

clayParams_fit.eps_inf_clay = data.clay_eps_inf;
clayParams_fit.tau_clay     = p_fit(7);
clayParams_fit.alpha_clay   = p_fit(8);
%% Forward Model
freqs   = [50e6, 1.3e9, 3.2e9];
sigma_Sm = data.EC_dSm * 0.1;

resultsInv = soil_spectrum( ...
    data.VWC, data.SWC, data.EC_dSm, data.Tc, freqs, ...
    soilParams_fit, wparams_fit, tauBW_fit, clayParams_fit);
figure();plot(freq_vna,real(eps_vna),'k','LineWidth',2);
hold on; plot(cat(1,resultsInv.freq),real(cat(1,resultsInv.eps_complex)),'Color',[0.5,0.5,0.5],'linewidth',2)
plot(freq_vna,imag(eps_vna),'--k','LineWidth',2);
hold on; plot(cat(1,resultsInv.freq),imag(cat(1,resultsInv.eps_complex)),'Color',[0.5,0.5,0.5],'linewidth',2,'LineStyle','--')
grid on; grid minor
xlabel('Frequency (Hz)')
ylabel('Complex Permittivity')
legend('Real VNA','Real Model','Imag VNA', 'Imag Model')
title('dick')

plot_soil_dielectric_decomposition( ...
    soilParams, wparams, tauBW, clayParams, ...
    VWC, SWC, EC_dSm, Tc)
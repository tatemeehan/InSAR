%% Field Fox S11 Measurements to Complex Dielectric Permittivity
% Based on:
%Stuchly, M. A.; Stuchly, S. S. 
% “Coaxial line reflection method for measuring dielectric properties of0
% biological substances at radio and microwave frequencies — a review.”
% IEEE Transactions on Instrumentation and Measurement, vol. IM‑29, 1980, pp. 176‑183.
% https://doi.org/10.1186/1756-6649-14-3

%% Load Complex S11 
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
ws3Real = readtable([dataDir,'\S11_WS4_REAL.csv']);
ws3Imag = readtable([dataDir,'\S11_WS4_IMAG.csv']);

% Form Complex Number
freq = airReal.x_DATAFreq; % in Hz
s_air    = airReal.S11 - 1j*airImag.S11;
s_short  = shortReal.S11 - 1j*shortImag.S11;
s_water  = waterReal.S11 - 1j*waterImag.S11;
s_sample = ws3Real.S11 - 1j*ws3Imag.S11;
% Synthesize Short
s_air_cal = airRealCal.S11 - 1j*airImagCal.S11;
s_open_cal = openRealCal.S11 - 1j*openImagCal.S11;
s_short_cal = shortRealCal.S11 - 1j*shortImagCal.S11;
% % Transfer Function
% H = s_air_cal./s_open_cal;
% % H = movmean(H,101);
% s_short = H.*s_short_cal;
% % s_short = s_short_cal;

% Undo calibration to get "raw" S11, then apply new transfer
% s_short_uncal = s_short_cal ./ s_open_cal;  % Remove open cal (approx)
% s_short = s_short_uncal .* s_air_cal;       % Apply probe open as new ref

% s_short = s_air;
% Smooth Data
% s_air = smoothdata(s_air,2,"gaussian","omitmissing","SamplePoints",1001);
% s_short = smoothdata(s_short,2,"gaussian","omitmissing","SamplePoints",1001);
% s_water = smoothdata(s_water,2,"gaussian","omitmissing","SamplePoints",1001);
% s_sample = smoothdata(s_sample,2,"movmean","omitmissing","SamplePoints",1001);
% s_air = movmean(s_air,101);
% s_short = movmean(s_short,101);
% s_water = movmean(s_water,101);
% s_sample = movmean(s_sample,101);

s_air = smooth_permittivity(freq, s_air, 'lowess', struct('span', 0.175), false);
s_short = smooth_permittivity(freq, s_short, 'lowess', struct('span', 0.175), false);
s_water = smooth_permittivity(freq, s_water, 'lowess', struct('span', 0.175), false);
s_sample = smooth_permittivity(freq, s_sample, 'lowess', struct('span', 0.175), false);

%% Data Inversion
eps0 = 8.854187817e-12; % vacuum permittivity [F/m]
% 3pt
% [eps_complex3, freq] = invert_stuchly3pt(s_air, s_short, s_water, s_sample, freq);
% 2Pt
% [eps_complex, freq] = invert_stuchly2pt(s_air, s_water, s_sample, freq);
s_std = 0.005;
[eps_complex, eps_std] = invert_stuchly2pt_uncertainty(s_air, s_water, s_sample, freq, s_std);


% Remove NaN
nanIx = find(isnan(freq)|isnan(eps_complex)|freq<.5e9);
freq(nanIx) = []; eps_complex(nanIx) = []; eps_std(nanIx) = [];
minEpsImag = min((imag(eps_complex)));

if minEpsImag < 0
eps_complex = real(eps_complex)+1i.*(imag(eps_complex)-(1.25.*minEpsImag));
end
% eps_complex = smooth_permittivity(freq, eps_complex, 'lowess', struct('span', 0.25), false);

% result = fit_cole_cole_eps(freq, eps_complex);
% [params_fit, eps_model, residual] = fit_dual_cole_cole_eps(freq, eps_complex);
% [params_fit, eps_model] = fit_dual_cole_cole_eps(freq, eps_complex);
% freq50 = 50e6;
% eps_model50 = compute_dual_cole_cole(freq50, params_fit);
% freqq = 50e6:1e6:4e9;
% eps_model_q = compute_dual_cole_cole(freqq, params_fit);

% HydraGo Weighted Inversion
freq_hydra = 50e6;
eps_hydra = 12.2703 + 1i.*10.6567;
weight_hydra = 2.5;
sigma_hydra = 0.02933;
optsCost = struct('useRelativeError', true, ...
              'useRegularization', true, ...
              'lambda', 0.01, ...
              'plotModel', false, ...
             'sigma_hydra',sigma_hydra,...
             'relaxPenaltyWeight', 5);   % optional
[params_fit, eps_model, residual] = fit_dual_cole_cole_eps_hydraGo( ...
    freq, eps_complex, freq_hydra, eps_hydra, weight_hydra, optsCost, eps_std);
eps_relaxation = imag(eps_model)-params_fit(8)./(eps0.*2.*pi.*freq);
% eps_relaxation = imag(eps_model)-0.02533./(eps0.*2.*pi.*freq);

% Plot results
figure();
plot(freq/1e9, real(eps_complex), 'b.-', freq/1e9, real(eps_model), 'r--');
hold on;
plot(freq/1e9, imag(eps_complex), 'c.-', freq/1e9, imag(eps_model), 'm--');
xlabel('Frequency (GHz)');
ylabel('Permittivity');
legend('Real Data', 'Real Fit', 'Imag Data', 'Imag Fit');
grid on;

%% Forward Model Inversion Results and Plot
modFreq = [50e6:10e6:4e9]';
modEps = compute_dual_cole_cole(modFreq, params_fit);
imagModEps = imag(modEps);
condModEps = params_fit(8)./(eps0.*2.*pi.*modFreq);
relaxModEps = imagModEps-condModEps;
minRelaxModEps = min(relaxModEps);
if minRelaxModEps < 0
    relaxModEps = relaxModEps-minRelaxModEps;
    condModEps = condModEps+minRelaxModEps;
    imagModEps = relaxModEps+condModEps;
    modEps = complex(real(modEps),imagModEps);
end
% Plot results
figure();
plot(freq/1e9, real(eps_complex), 'k', modFreq/1e9, real(modEps), 'm--');
hold on;
plot(freq/1e9, imag(eps_complex), 'k', modFreq/1e9, imag(modEps), 'm--');
xlabel('Frequency (GHz)');
ylabel('Permittivity');
legend('Real Data', 'Real Fit', 'Imag Data', 'Imag Fit');
grid on;

% Plot results
figure();
plot(modFreq/1e9, condModEps, 'k--');
hold on;
plot(modFreq/1e9, relaxModEps, 'm--');
xlabel('Frequency (GHz)');
ylabel('Imaginary Permittivity');
legend('Conductivity Term', 'Dipolar Relaxation Term');
grid on;

figure(); subplot(1,2,1)
plot(freq/1e9, real(eps_complex), 'k','LineWidth',2)
hold on;
plot(modFreq/1e9, real(modEps), 'color',[0.5,0.5,0.5],'LineWidth',2);
hold on;
plot(freq/1e9, imag(eps_complex), '--k', 'LineWidth',2)
    plot(modFreq/1e9, imag(modEps), 'color',[0.5,0.5,0.5],'LineWidth',2,'LineStyle','--');
ylim([0 14])
xlabel('Frequency (GHz)');
ylabel('Complex Permittivity');
legend('Real Data', 'Real Fit', 'Imag Data', 'Imag Fit');
grid on;
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)

subplot(1,2,2)
plot(modFreq/1e9, condModEps, 'k-','LineWidth',2);
hold on;
plot(modFreq/1e9, relaxModEps, 'color',[0.5,0.5,0.5],'LineWidth',2);
ylim([0 14])

xlabel('Frequency (GHz)');
ylabel('Imaginary Permittivity');
legend('Conductivity Term', 'Dipolar Relaxation Term');
grid on;
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',12)

%% Signal Penetration Modeling
[alpha, delta, beta] = compute_attenuation_penetration(modFreq, modEps);
[dp,xa] = xtinction(modEps,modFreq);
%% Forward Cole-Cole Model
    eps_s_guess = 5:25;
    eps_inf_guess = 4.8;
    tau_guess = 1e-12:0.5e-11:1e-10;%2 / (2*pi*mean(freq));
    alpha_guess = 0.0:0.1:1;
    sigma_guess = 0.0:0.05:0.5;
    N = numel(eps_s_guess).*numel(eps_inf_guess).*numel(tau_guess).*numel(alpha_guess).*numel(sigma_guess);
    eps_model = zeros(numel(freq),N);
    % p = [eps_s_guess(:), eps_inf_guess(:), tau_guess(:), alpha_guess(:), sigma_guess(:)];
    iter = 1;
    for ii = 1:numel(eps_s_guess)
        for jj = 1:numel(eps_inf_guess)
            for kk = 1:numel(tau_guess)
                for ll = 1:numel(alpha_guess)
                    for mm = 1:numel(sigma_guess)
                        eps_s = eps_s_guess(ii);
                        eps_inf = eps_inf_guess(jj);
                        tau = tau_guess(kk);
                        alpha = alpha_guess(ll);
                        sigma = sigma_guess(mm);
                        omega = 2 * pi * freq;
                        jw_tau_alpha = (1 + (1i * omega * tau).^(1 - alpha));
                        eps_model(:,iter) = eps_inf + (eps_s - eps_inf) ./ ...
                            jw_tau_alpha - sigma ./ (1i.*omega * eps0);
                        iter = iter+1;
                    end
                end
            end
        end
    end
    figure();
    for nn = 1000:1000:N
subplot(2,1,1)
 hold on;

plot(freq,real(eps_model(:,nn)))
subplot(2,1,2)
 hold on;

plot(freq,imag(eps_model(:,nn)))
    end

    %% Dual cole-cole
    eps_inf  = 4.8;
    eps_s1   = 10:5:30;
    eps_s2   = 20:5:35;
    tau1     = logspace(-11, -9, 5);
    tau2     = logspace(-10, -8, 5);
    alpha1   = 0:0.3:1;
    alpha2   = 0:0.3:1;
    sigma    = 0:0.25:1.5;

    N = numel(eps_inf).*numel(eps_s1).*numel(eps_s2).*numel(tau1).*numel(tau2).*numel(alpha1).*numel(alpha2).*numel(sigma);
    eps_model = zeros(numel(freq),N);
    eps_double_relax = eps_model;
    eps_cond = eps_model;
    iter = 1;
    for ii = 1:numel(eps_inf)
        for jj = 1:numel(eps_s1)
            for kk = 1:numel(eps_s2)
                for ll = 1:numel(tau1)
                    for mm = 1:numel(tau2)
                        for nn = 1:numel(alpha1)
                            for oo = 1:numel(alpha2)
                                for pp = 1:numel(sigma)
                                    tmp_eps_inf = eps_inf(ii);
                                    tmp_eps_s1 = eps_s1(jj);
                                    tmp_eps_s2 = eps_s2(kk);
                                    tmp_tau1 = tau1(ll);
                                    tmp_tau2 = tau1(mm);
                                    tmp_alpha1 = alpha1(nn);
                                    tmp_alpha2 = alpha2(oo);
                                    tmp_sigma = sigma(pp);
                                    % Compute relaxation contributions
                                    delta_eps1 = tmp_eps_s1 - tmp_eps_inf;
                                    delta_eps2 = tmp_eps_s2 - tmp_eps_s1;

                                    term1 = delta_eps1 ./ (1 - (1i * omega * tmp_tau1).^(1 - tmp_alpha1));
                                    term2 = delta_eps2 ./ (1 - (1i * omega * tmp_tau2).^(1 - tmp_alpha2));
                                    cond  = 1i.*tmp_sigma ./ (omega * eps0);
                                    eps_double_relax(:,iter) = term1+term2;
                                    eps_cond(:,iter) = cond;
                                    eps_model(:,iter) = eps_inf + term1 + term2 + cond;
                                    iter = iter+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
%     figure();
%     for qq = 1000:1000:N
% subplot(2,1,1)
%  hold on;
% 
% plot(freq,real(eps_model(:,qq)))
% subplot(2,1,2)
%  hold on;
% 
% plot(freq,imag(eps_model(:,qq)))
%     end

        figure();
    qq = 10000;
plot(freq, imag(eps_double_relax(:,qq)), 'b--', 'DisplayName', 'Relaxation Only')
hold on
plot(freq, imag(eps_cond(:,qq)), 'r--', 'DisplayName', 'Conductivity Only')
plot(freq, imag(eps_model(:,qq)), 'k-', 'DisplayName', 'Total Model')
% 
    legend

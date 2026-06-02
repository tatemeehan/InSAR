
% SWE23_m = snowDepth23  .* snowDensity23 ./ 1000;   % m water equivalent
% SWE26_m = snowDepth26 .* ML.GPR.Density.mean ./ 1000;   % m water equivalent

SWE23_m = snowDepth23  .* 300 ./ 1000;   % m water equivalent
SWE26_m = snowDepth26 .* 300 ./ 1000;   % m water equivalent

deltaSWE_mm = (SWE26_m - SWE23_m) .* 1000;  % mm water equivalent

% rhoMean = 0.5 .* (snowDensity23 + ML.GPR.Density.mean);
rhoMean = 300.*ones(size(snowDensity23));

phi_oveis = oveis_phi_wrapped( ...
    deltaSWE_mm, ...
    incidenceAngle0, ...
    f, ...
    rhoMean, ...
    'Sign', +1, ...
    'WrapOutput', true);

%%
phi_oveis_pass = oveis_phi_from_passes( ...
    snowDepth23, 300.*ones(size(snowDensity23)), ...
    snowDepth26, ML.GPR.Density.mean, ...
    incidenceAngle0, ...
    f, ...
    'Sign', +1, ...
    'WrapOutput', true);

%%

win = 11;

phi_oveis_vmf = insar.vector_median_filter(phi_oveis_pass, win);

statsOveis = phase_model_bias_stats( ...
    phi_uas_vmf, ...
    phi_oveis_vmf, ...
    referenceMask, ...
    referenceMask);

statsSub = phase_model_bias_stats( ...
    phi_uas_vmf, ...
    phi_modelNoEas_vmf, ...
    referenceMask, ...
    referenceMask);

fprintf('\nOveis/refraction-limit comparison\n');

fprintf('Oveis-style refraction model:\n');
fprintf('  bias = %.4f rad\n', statsOveis.bias);
fprintf('  R = %.3f\n', statsOveis.R);
fprintf('  circStd = %.3f rad\n', statsOveis.circStd);
fprintf('  median |res| = %.3f rad\n', statsOveis.medAbs);

fprintf('Our subsurface model:\n');
fprintf('  bias = %.4f rad\n', statsSub.bias);
fprintf('  R = %.3f\n', statsSub.R);
fprintf('  circStd = %.3f rad\n', statsSub.circStd);
fprintf('  median |res| = %.3f rad\n', statsSub.medAbs);

%%
fitMask  = referenceMask;
evalMask = coh > 0.75 & isfinite(phi_uas_vmf);

statsOveis_eval = phase_model_bias_stats( ...
    phi_uas_vmf, ...
    phi_oveis_vmf, ...
    fitMask, ...
    evalMask);

statsSub_eval = phase_model_bias_stats( ...
    phi_uas_vmf, ...
    phi_modelNoEas_vmf, ...
    fitMask, ...
    evalMask);

fprintf('\nIndependent eval-mask comparison\n');

fprintf('Oveis model:\n');
fprintf('  bias = %.4f rad\n', statsOveis_eval.bias);
fprintf('  R = %.3f\n', statsOveis_eval.R);
fprintf('  circStd = %.3f rad\n', statsOveis_eval.circStd);
fprintf('  median |res| = %.3f rad\n', statsOveis_eval.medAbs);

fprintf('Subsurface model:\n');
fprintf('  bias = %.4f rad\n', statsSub_eval.bias);
fprintf('  R = %.3f\n', statsSub_eval.R);
fprintf('  circStd = %.3f rad\n', statsSub_eval.circStd);
fprintf('  median |res| = %.3f rad\n', statsSub_eval.medAbs);

%%
HS26 = snowDepth26;
HS23 = snowDepth23;
rhoConstVals = 200:10:450;

statsOveisConst = struct([]);

for ir = 1:numel(rhoConstVals)

    rhoConst = rhoConstVals(ir);

    % Self-consistent constant-density SWE change
    deltaSWE_mm_const = (HS26 - HS23) .* rhoConst;

    phi_oveis_const = oveis_phi_wrapped( ...
        deltaSWE_mm_const, ...
        incidenceAngle0, ...
        f, ...
        rhoConst, ...
        'Sign', +1, ...
        'WrapOutput', true);

    phi_oveis_const_vmf = insar.vector_median_filter(phi_oveis_const, 11);

    stats = phase_model_bias_stats( ...
        phi_uas_vmf, ...
        phi_oveis_const_vmf, ...
        referenceMask, ...
        referenceMask);

    statsOveisConst(ir).rhoConst = rhoConst;
    statsOveisConst(ir).bias = stats.bias;
    statsOveisConst(ir).R = stats.R;
    statsOveisConst(ir).circStd = stats.circStd;
    statsOveisConst(ir).medAbs = stats.medAbs;

end

Rvals = [statsOveisConst.R];
Mvals = [statsOveisConst.medAbs];

[~, iBestR] = max(Rvals);
[~, iBestM] = min(Mvals);

fprintf('Best constant rho by R: %.1f kg/m3, R = %.3f\n', ...
    statsOveisConst(iBestR).rhoConst, statsOveisConst(iBestR).R);

fprintf('Best constant rho by medAbs: %.1f kg/m3, medAbs = %.3f rad\n', ...
    statsOveisConst(iBestM).rhoConst, statsOveisConst(iBestM).medAbs);

figure('Color','w');
plot(rhoConstVals, Rvals, 'k', 'LineWidth', 2);
xlabel('Constant density (kg/m^3)');
ylabel('R');
title('Oveis constant-density sensitivity');
grid on;

figure('Color','w');
plot(rhoConstVals, Mvals, 'k', 'LineWidth', 2);
xlabel('Constant density (kg/m^3)');
ylabel('Median |residual| (rad)');
title('Oveis constant-density sensitivity');
grid on;

%%
rho26 = ML.GPR.Density.mean;
rho23 = snowDensity23;
% -------------------------
% Common
% -------------------------
rhoConst = 300;
wrap = @(x) angle(exp(1i .* x));

% SWE change from spatial density states
deltaSWE_mm_spatial = HS26 .* rho26 - HS23 .* rho23;

% -------------------------
% A) Spatial-density Oveis
% -------------------------
rhoMean = 0.5 .* (rho23 + rho26);

phi_oveis_spatial = oveis_phi_wrapped( ...
    deltaSWE_mm_spatial, incidenceAngle0, f, rhoMean, ...
    'Sign', +1, 'WrapOutput', true);

% -------------------------
% B) Hybrid constant-density conversion
%    spatial SWE change, constant rho conversion
% -------------------------
phi_oveis_hybrid300 = oveis_phi_wrapped( ...
    deltaSWE_mm_spatial, incidenceAngle0, f, rhoConst, ...
    'Sign', +1, 'WrapOutput', true);

% -------------------------
% C) Self-consistent constant-density depth-change model
%    constant rho used to build SWE and convert back
% -------------------------
deltaSWE_mm_const300 = (HS26 - HS23) .* rhoConst;

phi_oveis_const300 = oveis_phi_wrapped( ...
    deltaSWE_mm_const300, incidenceAngle0, f, rhoConst, ...
    'Sign', +1, 'WrapOutput', true);

win = 11;

phi_oveis_spatial_vmf   = insar.vector_median_filter(phi_oveis_spatial, win);
phi_oveis_hybrid300_vmf = insar.vector_median_filter(phi_oveis_hybrid300, win);
phi_oveis_const300_vmf  = insar.vector_median_filter(phi_oveis_const300, win);

statsSpatial = phase_model_bias_stats( ...
    phi_uas_vmf, phi_oveis_spatial_vmf, referenceMask, referenceMask);

statsHybrid300 = phase_model_bias_stats( ...
    phi_uas_vmf, phi_oveis_hybrid300_vmf, referenceMask, referenceMask);

statsConst300 = phase_model_bias_stats( ...
    phi_uas_vmf, phi_oveis_const300_vmf, referenceMask, referenceMask);

fprintf('\nOveis model variants\n');

fprintf('Spatial density SWE + spatial mean rho:\n');
fprintf('  bias = %.4f, R = %.3f, medAbs = %.3f\n', ...
    statsSpatial.bias, statsSpatial.R, statsSpatial.medAbs);

fprintf('Hybrid: spatial SWE + constant rho=300:\n');
fprintf('  bias = %.4f, R = %.3f, medAbs = %.3f\n', ...
    statsHybrid300.bias, statsHybrid300.R, statsHybrid300.medAbs);

fprintf('Self-consistent constant rho=300 depth-change:\n');
fprintf('  bias = %.4f, R = %.3f, medAbs = %.3f\n', ...
    statsConst300.bias, statsConst300.R, statsConst300.medAbs);

dphi_const_vs_hybrid = wrap(phi_oveis_const300 - phi_oveis_hybrid300);

fprintf('\nConstant 300 self-consistent vs hybrid 300\n');
fprintf('  median diff = %.4f rad\n', ...
    median(dphi_const_vs_hybrid(:), 'omitnan'));
fprintf('  P05/P95 diff = %.4f / %.4f rad\n', ...
    prctile(dphi_const_vs_hybrid(:), 5), ...
    prctile(dphi_const_vs_hybrid(:), 95));
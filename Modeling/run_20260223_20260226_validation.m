%% run_20260223_20260226_validation
% Clean validation harness for UAS InSAR phase model comparison.
%
% Required workspace variables:
%   phi_uas_vmf
%   insarData.phaseFillTrustedMask
%   demData.R
%   Coords.utmX, Coords.utmY
%   coh
%   snowDepth23, snowDepth26
%   snowDensity23
%   ML.GPR.Density.mean
%   incidenceAngle0
%   f
%   phi_model_vmf          % all components, filtered
%   phi_modelNoEas_vmf     % no-Eas / subsurface, filtered
%
% Optional:
%   Esub23, Esub26, Eas23, Eas26 for surface-scale sweep

clear statsRows

wrap = @(x) angle(exp(1i .* x));
win = 11;

%% ------------------------------------------------------------------------
% 1. Build masks
% -------------------------------------------------------------------------

referenceMaskInterp = mapinterp( ...
    insarData.phaseFillTrustedMask, ...
    demData.R, ...
    Coords.utmX, ...
    Coords.utmY);

referenceMask = isfinite(referenceMaskInterp) & referenceMaskInterp > 0.5;

evalMask_ref = referenceMask;

evalMask_coh075 = ...
    isfinite(phi_uas_vmf) & ...
    isfinite(coh) & coh > 0.75;

fprintf('Validation masks:\n');
fprintf('  referenceMask pixels = %d\n', nnz(evalMask_ref));
fprintf('  coh > 0.75 pixels    = %d\n', nnz(evalMask_coh075));

%% ------------------------------------------------------------------------
% 2. Build Oveis / refraction-limit variants
% -------------------------------------------------------------------------

HS23 = snowDepth23;
HS26 = snowDepth26;

rho23 = snowDensity23;
rho26 = ML.GPR.Density.mean;

rhoConst = 300;

% A) Spatial SWE + spatial mean density
deltaSWE_mm_spatial = HS26 .* rho26 - HS23 .* rho23;
rhoMean = 0.5 .* (rho23 + rho26);

phi_oveis_spatial = oveis_phi_wrapped( ...
    deltaSWE_mm_spatial, ...
    incidenceAngle0, ...
    f, ...
    rhoMean, ...
    'Sign', +1, ...
    'WrapOutput', true);

% B) Hybrid: spatial SWE + constant density conversion
phi_oveis_hybrid300 = oveis_phi_wrapped( ...
    deltaSWE_mm_spatial, ...
    incidenceAngle0, ...
    f, ...
    rhoConst, ...
    'Sign', +1, ...
    'WrapOutput', true);

% C) Self-consistent constant density depth-change model
deltaSWE_mm_const300 = (HS26 - HS23) .* rhoConst;

phi_oveis_const300 = oveis_phi_wrapped( ...
    deltaSWE_mm_const300, ...
    incidenceAngle0, ...
    f, ...
    rhoConst, ...
    'Sign', +1, ...
    'WrapOutput', true);

% Filter Oveis products
phi_oveis_spatial_vmf   = insar.vector_median_filter(phi_oveis_spatial, win);
phi_oveis_hybrid300_vmf = insar.vector_median_filter(phi_oveis_hybrid300, win);
phi_oveis_const300_vmf  = insar.vector_median_filter(phi_oveis_const300, win);

%% ------------------------------------------------------------------------
% 3. Score model set
% -------------------------------------------------------------------------

models = struct([]);

models(1).name = 'Complex_All';
models(1).phi  = phi_model_vmf;

models(2).name = 'Complex_NoEas_EsgEgv';
models(2).phi  = phi_modelNoEas_vmf;

models(3).name = 'Oveis_SpatialDensity';
models(3).phi  = phi_oveis_spatial_vmf;

models(4).name = 'Oveis_HybridRho300';
models(4).phi  = phi_oveis_hybrid300_vmf;

models(5).name = 'Oveis_ConstRho300';
models(5).phi  = phi_oveis_const300_vmf;

statsRows = {};
iRow = 0;

for im = 1:numel(models)

    modelName = models(im).name;
    phiModel = models(im).phi;

    % Reference-mask evaluation
    statsRef = score_phase_model( ...
        phi_uas_vmf, ...
        phiModel, ...
        referenceMask, ...
        evalMask_ref);

    iRow = iRow + 1;
    statsRows{iRow,1} = make_stats_row( ...
        modelName, ...
        'referenceMask', ...
        statsRef);

    % Independent high-coherence evaluation
    statsCoh = score_phase_model( ...
        phi_uas_vmf, ...
        phiModel, ...
        referenceMask, ...
        evalMask_coh075);

    iRow = iRow + 1;
    statsRows{iRow,1} = make_stats_row( ...
        modelName, ...
        'coh075', ...
        statsCoh);

end

statsStruct = vertcat(statsRows{:});
statsTable = struct2table(statsStruct);

disp(statsTable_sorted(:, { ...
    'evalMaskName', ...
    'modelName', ...
    'bias', ...
    'R_bias', ...
    'phaseCorr_complex_bias', ...
    'circStd_bias', ...
    'medAbs_bias', ...
    'nEval'}));

%% ------------------------------------------------------------------------
% 4. Constant-density Oveis sweep
% -------------------------------------------------------------------------

rhoConstVals = 200:10:450;

rhoSweepRows = struct([]);

for ir = 1:numel(rhoConstVals)

    rhoC = rhoConstVals(ir);

    deltaSWE_mm_const = (HS26 - HS23) .* rhoC;

    phi_oveis_const = oveis_phi_wrapped( ...
        deltaSWE_mm_const, ...
        incidenceAngle0, ...
        f, ...
        rhoC, ...
        'Sign', +1, ...
        'WrapOutput', true);

    phi_oveis_const_vmf = insar.vector_median_filter(phi_oveis_const, win);

    stats = score_phase_model( ...
        phi_uas_vmf, ...
        phi_oveis_const_vmf, ...
        referenceMask, ...
        evalMask_ref);

    rhoSweepRows(ir).rhoConst = rhoC;
    rhoSweepRows(ir).bias = stats.bias;
    rhoSweepRows(ir).R_bias = stats.R_bias;
    rhoSweepRows(ir).circStd_bias = stats.circStd_bias;
    rhoSweepRows(ir).medAbs_bias = stats.medAbs_bias;
    rhoSweepRows(ir).nEval = stats.nEval;

end

rhoSweepTable = struct2table(rhoSweepRows);

[~, iBestR] = max(rhoSweepTable.R_bias);
[~, iBestM] = min(rhoSweepTable.medAbs_bias);

fprintf('\nOveis constant-density sweep:\n');
fprintf('  Best rho by R      = %.1f kg/m3, R = %.3f\n', ...
    rhoSweepTable.rhoConst(iBestR), rhoSweepTable.R_bias(iBestR));
fprintf('  Best rho by medAbs = %.1f kg/m3, medAbs = %.3f rad\n', ...
    rhoSweepTable.rhoConst(iBestM), rhoSweepTable.medAbs_bias(iBestM));

figure('Color','w');
plot(rhoSweepTable.rhoConst, rhoSweepTable.R_bias, 'k', 'LineWidth', 2);
xlabel('Constant density (kg/m^3)');
ylabel('Bias-corrected R');
title('Oveis constant-density sensitivity');
grid on;

figure('Color','w');
plot(rhoSweepTable.rhoConst, rhoSweepTable.medAbs_bias, 'k', 'LineWidth', 2);
xlabel('Constant density (kg/m^3)');
ylabel('Median |bias-corrected residual| (rad)');
title('Oveis constant-density sensitivity');
grid on;

%% ------------------------------------------------------------------------
% 5. Optional: surface phasor amplitude sweep
% -------------------------------------------------------------------------

if exist('Esub23','var') && exist('Esub26','var') && ...
        exist('Eas23','var') && exist('Eas26','var')

    aVals = 0:0.02:1;

    surfRows = struct([]);

    for ia = 1:numel(aVals)

        a = aVals(ia);

        E23 = Esub23 + a .* Eas23;
        E26 = Esub26 + a .* Eas26;

        phiSurfMix = angle(E26 .* conj(E23));
        phiSurfMix_vmf = insar.vector_median_filter(phiSurfMix, win);

        stats = score_phase_model( ...
            phi_uas_vmf, ...
            phiSurfMix_vmf, ...
            referenceMask, ...
            evalMask_ref);

        surfRows(ia).surfaceScale = a;
        surfRows(ia).bias = stats.bias;
        surfRows(ia).R_bias = stats.R_bias;
        surfRows(ia).circStd_bias = stats.circStd_bias;
        surfRows(ia).medAbs_bias = stats.medAbs_bias;
        surfRows(ia).nEval = stats.nEval;

    end

    surfaceScaleTable = struct2table(surfRows);

    [~, iBestR] = max(surfaceScaleTable.R_bias);
    [~, iBestM] = min(surfaceScaleTable.medAbs_bias);

    fprintf('\nSurface phasor scale sweep:\n');
    fprintf('  Best a by R      = %.3f, R = %.3f\n', ...
        surfaceScaleTable.surfaceScale(iBestR), ...
        surfaceScaleTable.R_bias(iBestR));
    fprintf('  Best a by medAbs = %.3f, medAbs = %.3f rad\n', ...
        surfaceScaleTable.surfaceScale(iBestM), ...
        surfaceScaleTable.medAbs_bias(iBestM));

    figure('Color','w');
    plot(surfaceScaleTable.surfaceScale, surfaceScaleTable.R_bias, 'k', 'LineWidth', 2);
    xlabel('Surface phasor scale a');
    ylabel('Bias-corrected R');
    title('Surface phasor scale sensitivity');
    grid on;

    figure('Color','w');
    plot(surfaceScaleTable.surfaceScale, surfaceScaleTable.medAbs_bias, 'k', 'LineWidth', 2);
    xlabel('Surface phasor scale a');
    ylabel('Median |bias-corrected residual| (rad)');
    title('Surface phasor scale sensitivity');
    grid on;

end

%% ------------------------------------------------------------------------
% Phase Bias as SWE error
rhoRef = 300;

for ii = 1:height(statsTable)

    err = phase_bias_to_oveis_error( ...
        statsTable.bias(ii), incidenceAngle0, f, rhoRef);

    statsTable.eq_dH_median_m(ii,1) = err.dH_m_median;
    statsTable.eq_abs_dH_median_m(ii,1) = err.abs_dH_m_median;
    statsTable.eq_SWE_median_mm(ii,1) = err.SWE_mm_median;
    statsTable.eq_abs_SWE_median_mm(ii,1) = err.abs_SWE_mm_median;

end

disp(statsTable(:, { ...
    'evalMaskName', ...
    'modelName', ...
    'bias', ...
    'eq_SWE_median_mm', ...
    'R_bias', ...
    'medAbs_bias'}));
%% ------------------------------------------------------------------------
% 6. Save validation products
% -------------------------------------------------------------------------

outDir = 'D:\IDALS\2026\forward_model_runs\validation_clean\';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

save(fullfile(outDir, 'validation_summary_tables.mat'), ...
    'statsTable', ...
    'rhoSweepTable', ...
    '-v7.3');

if exist('surfaceScaleTable','var')
    save(fullfile(outDir, 'surface_scale_table.mat'), ...
        'surfaceScaleTable', ...
        '-v7.3');
end

fprintf('\nSaved validation tables to:\n%s\n', outDir);

%% ------------------------------------------------------------------------
% Local helper
% -------------------------------------------------------------------------
function row = make_stats_row(modelName, evalMaskName, stats)

row = struct( ...
    'modelName', char(modelName), ...
    'evalMaskName', char(evalMaskName), ...
    'bias', stats.bias, ...
    'R_raw', stats.R_raw, ...
    'phaseCorr_complex_raw', stats.phaseCorr_complex_raw, ...
    'circStd_raw', stats.circStd_raw, ...
    'medAbs_raw', stats.medAbs_raw, ...
    'R_bias', stats.R_bias, ...
    'phaseCorr_complex_bias', stats.phaseCorr_complex_bias, ...
    'circStd_bias', stats.circStd_bias, ...
    'medAbs_bias', stats.medAbs_bias, ...
    'nFit', stats.nFit, ...
    'nEval', stats.nEval);

end
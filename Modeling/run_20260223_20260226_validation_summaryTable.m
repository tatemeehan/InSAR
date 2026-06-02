figure();
imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,((cosd(LiDAR.C.aspect+45)+sind(LiDAR.C.aspect+45))).*sind(2.5.*LiDAR.C.slope));colormap(bone)
utils.freezeColors; hold on;
hI = imagesc(demData.X(1,:)./1000,demData.Y(:,1)./1000,phi_oveis_hybrid300_vmf,'AlphaData',0.625);daspect([1,1,1]);colormap([[1 1 1];cmap]);hc=colorbar;
ylabel(hc,'\Delta \phi (rad)','fontname','serif','fontweight','bold','fontsize',14)

% xlabel('Longitude');ylabel('Latitude');
xlabel('Easting (km)');ylabel('Northing (km)');
clim([-pi pi])
set(gca,'YDir','normal','fontname','serif','fontweight','bold','fontsize',14)

ax = ancestor(hI, 'axes');
ax.XAxis.Exponent = 0;
xtickformat('%.1f')
ax.YAxis.Exponent = 0;
ytickformat('%.1f')

%%
rawFinite  = isfinite(phi_uas);
filtFinite = isfinite(phi_uas_vmf);

refMask = logical(referenceMask);
cohMask = isfinite(coh) & coh > 0.75;

fprintf('\nFilter support diagnostics\n');

fprintf('Raw finite pixels:        %d\n', nnz(rawFinite));
fprintf('Filtered finite pixels:   %d\n', nnz(filtFinite));
fprintf('Filtered/raw fraction:    %.3f\n', nnz(filtFinite) / nnz(rawFinite));

fprintf('\nReference mask support\n');
fprintf('referenceMask pixels:             %d\n', nnz(refMask));
fprintf('filtered finite in referenceMask: %.3f\n', ...
    nnz(filtFinite & refMask) / nnz(refMask));

fprintf('\nCoherence mask support\n');
fprintf('coh > 0.75 pixels:                %d\n', nnz(cohMask));
fprintf('filtered finite in coh > 0.75:    %.3f\n', ...
    nnz(filtFinite & cohMask) / nnz(cohMask));

fprintf('\nOverlap of filtered support with referenceMask\n');
fprintf('fraction of filtered pixels inside referenceMask: %.3f\n', ...
    nnz(filtFinite & refMask) / nnz(filtFinite));

%%
figure('Color','w');
imagesc(filtFinite);
axis image;
colorbar;
title('Finite support after phase filtering');

figure('Color','w');
imagesc(refMask);
axis image;
colorbar;
title('Reference mask');

figure('Color','w');
imagesc(cohMask & filtFinite);
axis image;
colorbar;
title('coh > 0.75 and finite after filtering');
%%
wrap = @(x) angle(exp(1i .* x));

modelPhases = struct();
modelPhases.UAS = phi_uas_vmf;
modelPhases.Complex_All = phi_model_vmf;
modelPhases.Complex_NoEas = phi_modelNoEas_vmf;
modelPhases.Oveis_HybridRho300 = phi_oveis_hybrid300_vmf;

maskNames = {'referenceMask','coh075'};
masks = {referenceMask, coh > 0.75 & isfinite(phi_uas_vmf)};

rows = {};
iRow = 0;

names = fieldnames(modelPhases);

for imask = 1:numel(masks)

    thisMask = logical(masks{imask});

    for in = 1:numel(names)

        name = names{in};
        phi = modelPhases.(name);

        valid = thisMask & isfinite(phi);

        phiv = wrap(phi(valid));

        circMean = angle(mean(exp(1i .* phiv), 'omitnan'));
        circR = abs(mean(exp(1i .* phiv), 'omitnan'));

        % Ordinary median of wrapped phase values.
        % Fine for diagnostics, but remember it is not fully circular.
        medWrapped = median(phiv, 'omitnan');

        q = prctile(phiv, [5 95]);

        iRow = iRow + 1;
        rows{iRow,1} = struct( ...
            'maskName', maskNames{imask}, ...
            'phaseName', name, ...
            'circMean', circMean, ...
            'circR', circR, ...
            'medianWrapped', medWrapped, ...
            'p05', q(1), ...
            'p95', q(2), ...
            'n', nnz(valid));
    end
end

phaseSummaryTable = struct2table(vertcat(rows{:}));

disp(phaseSummaryTable);

%%
wrap = @(x) angle(exp(1i .* x));

modelPhases = struct();
modelPhases.UAS = phi_uas_vmf;
modelPhases.Complex_All = phi_model_vmf;
modelPhases.Complex_NoEas = phi_modelNoEas_vmf;
modelPhases.Oveis_HybridRho300 = phi_oveis_hybrid300_vmf;

maskNames = {'referenceMask','coh075'};
masks = {referenceMask, coh > 0.75 & isfinite(phi_uas_vmf)};

rows = {};
iRow = 0;

names = fieldnames(modelPhases);

phiUAS = phi_uas_vmf;

for imask = 1:numel(masks)

    thisMask = logical(masks{imask});

    for in = 1:numel(names)

        name = names{in};
        phi = modelPhases.(name);

        valid = thisMask & isfinite(phi) & isfinite(phiUAS);

        phiv = wrap(phi(valid));

        circMean = angle(mean(exp(1i .* phiv), 'omitnan'));
        circR = abs(mean(exp(1i .* phiv), 'omitnan'));

        medWrapped = median(phiv, 'omitnan');
        q = prctile(phiv, [5 95]);

        % ------------------------------------------------------------
        % Complex phase correlation relative to UAS
        % ------------------------------------------------------------
        zU = exp(1i .* phiUAS(valid));
        zM = exp(1i .* phi(valid));

        phaseCorr_complex = abs(sum(zU .* conj(zM), 'omitnan')) ./ ...
            sqrt(sum(abs(zU).^2, 'omitnan') .* sum(abs(zM).^2, 'omitnan'));

% Wrapped residual relative to UAS
% res = wrap(phiUAS(valid) - phi(valid));
res = wrap(phi(valid)-phiUAS(valid));


residualR = abs(mean(exp(1i .* res), 'omitnan'));
residualMedianAbs = median(abs(res), 'omitnan');
residualCircMean = angle(mean(exp(1i .* res), 'omitnan'));

% Bias-corrected residual
resBias = wrap(res - residualCircMean);

residualR_biasCorrected = abs(mean(exp(1i .* resBias), 'omitnan'));
residualMedianAbs_biasCorrected = median(abs(resBias), 'omitnan');

iRow = iRow + 1;
rows{iRow,1} = struct( ...
    'maskName', maskNames{imask}, ...
    'phaseName', name, ...
    'circMean', circMean, ...
    'circR_self', circR, ...
    'medianWrapped', medWrapped, ...
    'p05', q(1), ...
    'p95', q(2), ...
    'phaseCorr_complex_vs_UAS', phaseCorr_complex, ...
    'residualR_vs_UAS', residualR, ...
    'residualCircMean_vs_UAS', residualCircMean, ...
    'residualMedianAbs_vs_UAS', residualMedianAbs, ...
    'residualR_biasCorrected_vs_UAS', residualR_biasCorrected, ...
    'residualMedianAbs_biasCorrected_vs_UAS', residualMedianAbs_biasCorrected, ...
    'n', nnz(valid));
    end
end

phaseSummaryTable = struct2table(vertcat(rows{:}));

disp(phaseSummaryTable);

%%
models = struct();
models.Complex_All = phi_model_vmf;
models.Complex_NoEas = phi_modelNoEas_vmf;
models.Complex_Esg = phi_Esg_vmf;
models.Oveis_HybridRho300 = phi_oveis_hybrid300_vmf;
models.Oveis_Spatial = phi_oveis_spatial_vmf;
models.Oveis_Rho300 = phi_oveis_const300_vmf;




names = fieldnames(models);

for i = 1:numel(names)

    name = names{i};

    pstats = pearson_phase_bias_corrected( ...
        phi_uas_vmf, ...
        models.(name), ...
        referenceMask, ...
        referenceMask);

    fprintf('\n%s\n', name);
    fprintf('  bias = %.3f rad\n', pstats.bias);
    fprintf('  Pearson wrapped, bias-corrected = %.3f\n', pstats.r_wrapped_bias);
    fprintf('  Pearson cos component = %.3f\n', pstats.r_cos_bias);
    fprintf('  Pearson sin component = %.3f\n', pstats.r_sin_bias);
    fprintf('  Mean component Pearson = %.3f\n', pstats.r_component_mean_bias);

end

%% Surface/subsurface amplitude dominance
Eas23 = out23Cases.all.Eas;
Eas26 = out26Cases.all.Eas;

Esg23 = out23Cases.all.Esg;
Esg26 = out26Cases.all.Esg;

A_surface_ifg = sqrt(abs(Eas23) .* abs(Eas26));
A_esg_ifg     = sqrt(abs(Esg23) .* abs(Esg26));

ratioLog_Eas_Esg = log10(A_surface_ifg ./ A_esg_ifg);
ratioLog_Eas_Esg(~isfinite(ratioLog_Eas_Esg)) = NaN;

valid = referenceMask & isfinite(ratioLog_Eas_Esg);

fprintf('log10(|Eas|/|Esg|) over referenceMask:\n');
fprintf('  median = %.3f\n', median(ratioLog_Eas_Esg(valid), 'omitnan'));
fprintf('  P05/P95 = %.3f / %.3f\n', ...
    prctile(ratioLog_Eas_Esg(valid),5), ...
    prctile(ratioLog_Eas_Esg(valid),95));

fprintf('Fraction Eas < Esg = %.2f %%\n', ...
    100*nnz(valid & ratioLog_Eas_Esg < 0)/nnz(valid));

%% Does adding Eas degrade the model everywhere, or only where it is weak?
wrap = @(x) angle(exp(1i .* x));

phiData = phi_uas_vmf;

phiEsg = phi_Esg_vmf;       % or phi_modelNoEas_vmf
phiAll = phi_model_vmf;

ratio = ratioLog_Eas_Esg;

bins = {
    ratio < -0.5, 'Esg dominates';
    abs(ratio) <= 0.5, 'mixed';
    ratio > 0.5, 'Eas dominates'
};

for ib = 1:size(bins,1)

    binMask = bins{ib,1};
    binName = bins{ib,2};

    evalMask = referenceMask & binMask & ...
        isfinite(phiData) & isfinite(phiEsg) & isfinite(phiAll);

    resEsg = wrap(phiData(evalMask) - phiEsg(evalMask));
    resAll = wrap(phiData(evalMask) - phiAll(evalMask));

    R_Esg = abs(mean(exp(1i .* resEsg), 'omitnan'));
    R_All = abs(mean(exp(1i .* resAll), 'omitnan'));

    med_Esg = median(abs(resEsg), 'omitnan');
    med_All = median(abs(resAll), 'omitnan');

    fprintf('\n%s\n', binName);
    fprintf('  n = %d\n', nnz(evalMask));
    fprintf('  Esg: R = %.3f, medAbs = %.3f\n', R_Esg, med_Esg);
    fprintf('  All: R = %.3f, medAbs = %.3f\n', R_All, med_All);
end

%% Surface coherence / reference suppression test
valid = isfinite(coh) & isfinite(ratioLog_Eas_Esg);

subMask  = valid & ratioLog_Eas_Esg < -0.5;
mixMask  = valid & abs(ratioLog_Eas_Esg) <= 0.5;
surfMask = valid & ratioLog_Eas_Esg > 0.5;

fprintf('\nObserved coherence by predicted return regime:\n');
fprintf('  Esg dominated median coh = %.3f\n', median(coh(subMask), 'omitnan'));
fprintf('  mixed median coh         = %.3f\n', median(coh(mixMask), 'omitnan'));
fprintf('  Eas dominated median coh = %.3f\n', median(coh(surfMask), 'omitnan'));

%% Flip the sign of the surface phase change
wrap = @(x) angle(exp(1i .* x));

aVals = -1:0.02:1;

Esub23 = out23Cases.all.Esg + out23Cases.all.Egv;
Esub26 = out26Cases.all.Esg + out26Cases.all.Egv;

% Or use Esg only, since Egv appears negligible:
% Esub23 = out23Cases.all.Esg;
% Esub26 = out26Cases.all.Esg;

Eas23 = out23Cases.all.Eas;
Eas26 = out26Cases.all.Eas;

fitMask  = referenceMask;
evalMask = referenceMask;

for ia = 1:numel(aVals)

    a = aVals(ia);

    E23 = Esub23 + a .* Eas23;
    E26 = Esub26 + a .* Eas26;

    phiMix = angle(E26 .* conj(E23));
    phiMix_vmf = insar.true_vector_median_filter_phase(phiMix, 11, ...
        'Normalize', true, ...
        'UseParallel', true);

    stats = score_phase_model( ...
        phi_uas_vmf, ...
        phiMix_vmf, ...
        fitMask, ...
        evalMask);

    sweep(ia).a = a;
    sweep(ia).bias = stats.bias;
    sweep(ia).R = stats.R_bias;
    sweep(ia).medAbs = stats.medAbs_bias;
    sweep(ia).circStd = stats.circStd_bias;

end

surfaceSignSweep = struct2table(sweep);

[~, iBestR] = max(surfaceSignSweep.R);
[~, iBestM] = min(surfaceSignSweep.medAbs);

fprintf('\nSigned Eas sweep\n');
fprintf('Best a by R      = %.3f, R = %.3f, medAbs = %.3f\n', ...
    surfaceSignSweep.a(iBestR), surfaceSignSweep.R(iBestR), surfaceSignSweep.medAbs(iBestR));

fprintf('Best a by medAbs = %.3f, R = %.3f, medAbs = %.3f\n', ...
    surfaceSignSweep.a(iBestM), surfaceSignSweep.R(iBestM), surfaceSignSweep.medAbs(iBestM));

figure('Color','w');
yyaxis left
plot(surfaceSignSweep.a, surfaceSignSweep.R, 'k-', 'LineWidth', 2);
ylabel('Bias-corrected R');

yyaxis right
plot(surfaceSignSweep.a, surfaceSignSweep.medAbs, 'k--', 'LineWidth', 2);
ylabel('Median |bias-corrected residual| (rad)');

xlabel('Signed Eas coefficient a');
title('Signed air/snow surface phasor sensitivity');
grid on;
xline(0,'--');
xline(1,':');
xline(-1,':');

%%
%%
wrap = @(x) angle(exp(1i .* x));

% Existing Oveis product from old convention
phi_oveis_old_vmf = phi_oveis_hybrid300_vmf;

% New convention, if old convention was exactly opposite sign
phi_oveis_new_vmf = wrap(-phi_oveis_old_vmf);

statsOld = score_phase_model( ...
    phi_uas_vmf, ...
    phi_oveis_old_vmf, ...
    referenceMask, ...
    referenceMask);

statsNew = score_phase_model( ...
    phi_uas_vmf, ...
    phi_oveis_new_vmf, ...
    referenceMask, ...
    referenceMask);

fprintf('\nOveis sign convention test\n');

fprintf('Old Oveis convention:\n');
fprintf('  correction-to-model bias = %.3f rad\n', statsOld.bias);
fprintf('  R = %.3f\n', statsOld.R_bias);
fprintf('  medAbs = %.3f rad\n', statsOld.medAbs_bias);

fprintf('New Oveis convention:\n');
fprintf('  correction-to-model bias = %.3f rad\n', statsNew.bias);
fprintf('  R = %.3f\n', statsNew.R_bias);
fprintf('  medAbs = %.3f rad\n', statsNew.medAbs_bias);

pOld = pearson_phase_bias_corrected( ...
    phi_uas_vmf, phi_oveis_old_vmf, referenceMask, referenceMask);

pNew = pearson_phase_bias_corrected( ...
    phi_uas_vmf, phi_oveis_new_vmf, referenceMask, referenceMask);

fprintf('\nOveis Pearson sign convention test\n');

fprintf('Old Oveis:\n');
fprintf('  component Pearson = %.3f\n', pOld.r_component_mean_bias);

fprintf('New Oveis:\n');
fprintf('  component Pearson = %.3f\n', pNew.r_component_mean_bias);

%%
wrap = @(x) angle(exp(1i .* x));

% Pick representative scalar values
theta0 = 40;
f0 = f;
epsSnow0 = 1.6 + 0.01i;
epsSoil0 = 8 + 1i;
Ig0 = 0;  % or representative value, but Esg-only test can ignore volume

H1 = 0.80;
H2 = 0.90;

p = pars1;
p.snow_surface.A0 = 0;
p.soil_surface.A0 = pars1.soil_surface.A0;
p.soil_volume.enable = false;
p.As = 0;

o1 = insar_forward_profile_lut_pass_vec(theta0, f0, H1, epsSnow0, epsSoil0, Ig0, p);
o2 = insar_forward_profile_lut_pass_vec(theta0, f0, H2, epsSnow0, epsSoil0, Ig0, p);

phi_Esg_depth_test = angle(o2.E .* conj(o1.E));

fprintf('Complex Esg phase for +dH = %.4f rad\n', phi_Esg_depth_test);
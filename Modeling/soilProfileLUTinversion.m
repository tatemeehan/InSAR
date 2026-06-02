% ------------------------------------------------------------
% Build baseline soil profile for constraint
% ------------------------------------------------------------
soilprof_base = soilprof;
soilprof_base.enable = true;

% Important: use the same frozen/wet parameterization you believe gives
% near-surface eps' around 9 for wet partially frozen soil.
soilprof_base.freeze.enable = true;
soilprof_base.freeze.mode = 'physical';

% ------------------------------------------------------------
% Build VWC0 inversion LUT
% ------------------------------------------------------------
soilVWC0_LUT_surface = build_soil_VWC0_constraint_LUT( ...
    soilprof_base, ...
    f, ...
    theta_i, ...
    pars1.Lg, ...
    'Metric', 'surface', ...
    'VWC0Grid', linspace(0.05, soilprof_base.SWC, 1000), ...
    'Verbose', true);

% Optional: inspect the forward mapping
figure('Color','w');
plot(soilVWC0_LUT_surface.VWC0Grid, ...
     soilVWC0_LUT_surface.epsMetric, ...
     'k', 'LineWidth', 2);
xlabel('VWC0');
ylabel("Modeled near-surface soil \epsilon'");
title('Soil VWC0 constraint LUT');
grid on; grid minor;

% ------------------------------------------------------------
% Apply LUT to GPR/ML real soil permittivity raster
% ------------------------------------------------------------
epsRealSoil26 = ML.GPR.SoilPermittivity.mean;
validSoil = isfinite(epsRealSoil26) & epsRealSoil26 > 1;

[VWC0_eff_26, vwcDiag] = apply_soil_VWC0_constraint_LUT( ...
    epsRealSoil26, ...
    soilVWC0_LUT_surface, ...
    validSoil, ...
    'OutputClass', 'single', ...
    'ClipInput', true);

fprintf('\nVWC0 effective raster diagnostics\n');
fprintf('  nValid       = %d\n', vwcDiag.nValid);
fprintf('  nBelow LUT   = %d\n', vwcDiag.nBelowLUT);
fprintf('  nAbove LUT   = %d\n', vwcDiag.nAboveLUT);
fprintf('  VWC0 mean    = %.3f\n', vwcDiag.VWC0Mean);
fprintf('  VWC0 median  = %.3f\n', vwcDiag.VWC0Median);
fprintf('  VWC0 min/max = %.3f / %.3f\n', ...
    vwcDiag.VWC0Min, vwcDiag.VWC0Max);

% ------------------------------------------------------------
% Baseline assumption: soil fixed between passes
% ------------------------------------------------------------
VWC0_eff_23 = VWC0_eff_26;

% ------------------------------------------------------------
% Validate inversion consistency
% ------------------------------------------------------------
val = validate_soil_VWC0_constraint( ...
    epsRealSoil26, ...
    VWC0_eff_26, ...
    soilVWC0_LUT_surface, ...
    validSoil);
%% Write Raster VWC0
outDir = 'D:\IDALS\2026\soil_complex_states\';

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

geotiffwrite( ...
    fullfile(outDir, '20260226_MCS_VWC0_eff_GPRconstrained_surface.tif'), ...
    VWC0_eff_26, ...
    Coords.R, ...
    'CoordRefSysCode', Coords.EPSG);

geotiffwrite( ...
    fullfile(outDir, '20260223_MCS_VWC0_eff_GPRconstrained_surface_FIXEDSOIL.tif'), ...
    VWC0_eff_23, ...
    Coords.R, ...
    'CoordRefSysCode', Coords.EPSG);

%% Test Depth-Weighted Alternatives
soilVWC0_LUT_shallow05 = build_soil_VWC0_constraint_LUT( ...
    soilprof_base, ...
    f, ...
    theta_i, ...
    pars1.Lg, ...
    'Metric', 'shallow', ...
    'AnchorDepth', 0.0125, ...
    'VWC0Grid', linspace(0.05, soilprof_base.SWC, 1000), ...
    'Verbose', true);

[VWC0_eff_26_shallow05, diag_shallow05] = apply_soil_VWC0_constraint_LUT( ...
    epsRealSoil26, ...
    soilVWC0_LUT_shallow05, ...
    validSoil, ...
    'OutputClass', 'single', ...
    'ClipInput', true);

% Compare
figure('Color','w');

histogram(VWC0_eff_26(validSoil), 100);
hold on;
histogram(VWC0_eff_26_shallow05(validSoil), 100);

xlabel('Effective VWC0');
ylabel('Count');
legend('Surface metric', 'Shallow 5 cm metric');
title('GPR-constrained VWC0 sensitivity');
grid on;

% This reveals that the surface interpretation of the VWC0 is correct.
% The GPR/ML soil permittivity map was interpreted as a 
% near-surface/interface real-permittivity constraint. 
% We therefore inverted for an effective surface water-content parameter, 
% VWC0, while holding the deeper profile structure fixed. 
% A 5 cm shallow-average interpretation drove the inversion to saturation 
% and was not used for the baseline.

%% Forward Soil Profile Integral LUT
VWC0Grid = linspace(0.05, 0.55, 500);
thetaGrid_deg = 15:0.25:85;

soilFwdLUT = build_soil_profile_forward_LUT( ...
    soilprof_base, ...
    f, ...
    thetaGrid_deg, ...
    VWC0Grid, ...
    pars1.Lg);

soilFwdLUT_fn = 'D:\IDALS\2026\soil_complex_states\soilFwdLUT_surface_fixedSoil_20260223_20260226.mat';

save(soilFwdLUT_fn, 'soilFwdLUT', '-v7.3');

fprintf('Saved final soil forward LUT to:\n%s\n', soilFwdLUT_fn);

soilFwdLUT_tables_fn = fullfile(outDir, ...
    'soilFwdLUT_tablesOnly_VWC0_theta_fixedSoil_surfaceMetric.mat');

save(soilFwdLUT_tables_fn, ...
    '-struct', 'soilFwdLUT', ...
    'VWC0Grid', ...
    'thetaGrid_deg', ...
    'f_Hz', ...
    'Lg_default', ...
    'IgTable', ...
    'epsSurfTable', ...
    'VWC0Min', ...
    'VWC0Max', ...
    'thetaMin', ...
    'thetaMax', ...
    '-v7.3');

fprintf('Saved soil LUT tables:\n%s\n', soilFwdLUT_tables_fn);

figure('Color','w');
imagesc(soilFwdLUT.thetaGrid_deg, soilFwdLUT.VWC0Grid, real(soilFwdLUT.IgTable));
set(gca,'YDir','normal');
xlabel('Incidence angle (deg)');
ylabel('VWC0');
title('Real(Ig) soil profile LUT');
colorbar;
grid on;

figure('Color','w');
imagesc(soilFwdLUT.thetaGrid_deg, soilFwdLUT.VWC0Grid, imag(soilFwdLUT.IgTable));
set(gca,'YDir','normal');
xlabel('Incidence angle (deg)');
ylabel('VWC0');
title('Imag(Ig) soil profile LUT');
colorbar;
grid on;

% Example test point
VWC0_test  = 0.2908;
theta_test = 40;   % deg, choose representative local incidence

% LUT value
Ig_lut = soilFwdLUT.FIgReal(VWC0_test, theta_test) + ...
    1i * soilFwdLUT.FIgImag(VWC0_test, theta_test);

epsSurf_lut = soilFwdLUT.FepsReal(VWC0_test) + ...
    1i * soilFwdLUT.FepsImag(VWC0_test);

% Direct value
prof_test = soilprof_base;
prof_test.enable = true;
prof_test.VWC0 = VWC0_test;

c = 299792458;
k0 = 2*pi*soilFwdLUT.f_Hz/c;
kx = k0 * sind(theta_test);

[Ig_direct, diag_direct] = soil_profile_integral( ...
    soilFwdLUT.f_Hz, ...
    kx, ...
    prof_test, ...
    soilFwdLUT.Lg_default);

epsSurf_direct = diag_direct.eps(1);

fprintf('\nSoil profile LUT validation\n');
fprintf('Ig LUT:      %.6g %+ .6gi\n', real(Ig_lut), imag(Ig_lut));
fprintf('Ig direct:   %.6g %+ .6gi\n', real(Ig_direct), imag(Ig_direct));
fprintf('Ig abs diff: %.6g\n', abs(Ig_lut - Ig_direct));

fprintf("epsSurf LUT:      %.6f %+ .6fi\n", real(epsSurf_lut), imag(epsSurf_lut));
fprintf("epsSurf direct:   %.6f %+ .6fi\n", real(epsSurf_direct), imag(epsSurf_direct));
fprintf("epsSurf abs diff: %.6g\n", abs(epsSurf_lut - epsSurf_direct));
%% Compute IG from LUT
% Load InSAR Pair for MCS 02/26 - 02/23/2026
matDir = 'D:\InSAR\MCS\GLSAR\Export\MCS__lband-mcs_2026_TMA_2ms_200MHz_20260226T184834_TMA_TX2_RX2_V_V_V__lband-mcs_2026_TMA_2ms_200MHz_20260223T224040_TMA_TX2_RX2_V_V_V__standard__20260414T152819\mats';
geomData = io.load_export_struct(fullfile(matDir,'geomData_standard.mat'), 'geomSave');
insarData = io.load_export_struct(fullfile(matDir,'insarData_standard.mat'), 'insarSave');

incidenceAngle = (geomData.sarGeometry{1}.incidence+geomData.sarGeometry{1}.incidence2)./2;

% Oops we Have a size Mismatch on the DEM used in InSAR processing and all
% of our Model Parameters. We will Interpolte the InSAR Geometry Here to
% the Size of the LiDARGPR etc. predictors
incidenceAngle0 = mapinterp(incidenceAngle,demData.R,Coords.utmX,Coords.utmY);
phz = mapinterp(insarData.phzWrappedReferenced,demData.R,Coords.utmX,Coords.utmY);
phzUnW =  mapinterp(insarData.phzReferencedFilled,demData.R,Coords.utmX,Coords.utmY);
coh = mapinterp(insarData.coherence,demData.R,Coords.utmX,Coords.utmY);
screenWrapped = mapinterp(insarData.refScreenWrapped,demData.R,Coords.utmX,Coords.utmY);
validSoilForward = isfinite(VWC0_eff_26) & isfinite(incidenceAngle0);

[Ig_soil, epsSoilSurf] = eval_soil_profile_forward_LUT( ...
    VWC0_eff_26, ...
    incidenceAngle0, ...
    soilFwdLUT, ...
    validSoilForward, ...
    1e6);

outDir = 'D:\IDALS\2026\soil_complex_states\';

geotiffwrite(fullfile(outDir, '20260223_20260226_fixedSoil_Ig_real.tif'), ...
    real(Ig_soil), Coords.R, 'CoordRefSysCode', Coords.EPSG);

geotiffwrite(fullfile(outDir, '20260223_20260226_fixedSoil_Ig_imag.tif'), ...
    imag(Ig_soil), Coords.R, 'CoordRefSysCode', Coords.EPSG);

geotiffwrite(fullfile(outDir, '20260223_20260226_fixedSoil_epsSurf_real.tif'), ...
    real(epsSoilSurf), Coords.R, 'CoordRefSysCode', Coords.EPSG);

geotiffwrite(fullfile(outDir, '20260223_20260226_fixedSoil_epsSurf_imag.tif'), ...
    imag(epsSoilSurf), Coords.R, 'CoordRefSysCode', Coords.EPSG);

%% Recostruct Complex Soil Raster
soilDir = 'D:\IDALS\2026\soil_complex_states\';

[IgReal, R] = readgeoraster(fullfile(soilDir, '20260223_20260226_fixedSoil_Ig_real.tif'));
[IgImag, ~] = readgeoraster(fullfile(soilDir, '20260223_20260226_fixedSoil_Ig_imag.tif'));

[epsSoilReal, ~] = readgeoraster(fullfile(soilDir, '20260223_20260226_fixedSoil_epsSurf_real.tif'));
[epsSoilImag, ~] = readgeoraster(fullfile(soilDir, '20260223_20260226_fixedSoil_epsSurf_imag.tif'));

Ig_soil_fixed = complex(single(IgReal), single(IgImag));
epsSoilSurf_fixed = complex(single(epsSoilReal), single(epsSoilImag));

clear IgReal IgImag epsSoilReal epsSoilImag

validSoil = isfinite(real(Ig_soil_fixed)) & ...
            isfinite(imag(Ig_soil_fixed)) & ...
            isfinite(real(epsSoilSurf_fixed)) & ...
            isfinite(imag(epsSoilSurf_fixed));

fprintf('Valid soil pixels: %d\n', nnz(validSoil));

fprintf("epsSoilSurf real mean/min/max = %.3f / %.3f / %.3f\n", ...
    mean(real(epsSoilSurf_fixed(validSoil)), 'omitnan'), ...
    min(real(epsSoilSurf_fixed(validSoil)), [], 'omitnan'), ...
    max(real(epsSoilSurf_fixed(validSoil)), [], 'omitnan'));

fprintf("epsSoilSurf imag mean/min/max = %.3f / %.3f / %.3f\n", ...
    mean(imag(epsSoilSurf_fixed(validSoil)), 'omitnan'), ...
    min(imag(epsSoilSurf_fixed(validSoil)), [], 'omitnan'), ...
    max(imag(epsSoilSurf_fixed(validSoil)), [], 'omitnan'));

fprintf("Ig abs mean/min/max = %.5g / %.5g / %.5g\n", ...
    mean(abs(Ig_soil_fixed(validSoil)), 'omitnan'), ...
    min(abs(Ig_soil_fixed(validSoil)), [], 'omitnan'), ...
    max(abs(Ig_soil_fixed(validSoil)), [], 'omitnan'));

%% Let's Run an Insar Test over a Tile
% rr = 1500:2200;
% cc = 1200:1900;

rr = 1:size(snowDepth26,1);
cc = 1:size(snowDepth26,2);

theta_tile = incidenceAngle0(rr,cc);

HS23_tile = snowDepth23(rr,cc);
HS26_tile = snowDepth26(rr,cc);

epsSnow23_tile = epsSnow23_cases{7}(rr,cc);
[epsSnow26Real,~,~,~,~,~,~,~,~] = io.readLidarTif('D:\IDALS\2026\0226\lidarGPR\20260226_MCS_GPR_RealPermittivity.tif');
[epsSnow26Imag,~,~,~,~,~,~,~,~] = io.readLidarTif('D:\IDALS\2026\0226\lidarGPR\20260226_MCS_GPR_ImagPermittivity.tif');
epsSnow26 = complex(epsSnow26Real,epsSnow26Imag);
epsSnow26_tile = epsSnow26(rr,cc);

Ig_tile = Ig_soil_fixed(rr,cc);
epsSoilSurf_tile = epsSoilSurf_fixed(rr,cc);

% Run Passes
out23_tile = insar_forward_profile_lut_pass_vec( ...
    theta_tile, ...
    f, ...
    HS23_tile, ...
    epsSnow23_tile, ...
    epsSoilSurf_tile, ...
    Ig_tile, ...
    pars1);

out26_tile = insar_forward_profile_lut_pass_vec( ...
    theta_tile, ...
    f, ...
    HS26_tile, ...
    epsSnow26_tile, ...
    epsSoilSurf_tile, ...
    Ig_tile, ...
    pars2);

%%  InSAR Model Output
phi_model_tile = angle(out26_tile.E .* conj(out23_tile.E));
amp_model_tile = abs(out26_tile.E .* conj(out23_tile.E));

amp1 = abs(out23_tile.E);
amp2 = abs(out26_tile.E);
ifgAmp = abs(out26_tile.E .* conj(out23_tile.E));

validPhase = isfinite(phi_model_tile) & ...
             isfinite(ifgAmp) & ...
             ifgAmp > prctile(ifgAmp(:), 10, 'all');

phi_masked = phi_model_tile;
phi_masked(~validPhase) = NaN;

figure('Color','w');
imagesc(phi_masked);
axis image;
colorbar;
clim([-pi pi]);
title('Modeled interferometric phase, amplitude-masked');

figure('Color','w');
imagesc(abs(out23_tile.E));
axis image;
colorbar;
title('|E_{23}|');

figure('Color','w');
imagesc(abs(out26_tile.E));
axis image;
colorbar;
title('|E_{26}|');

figure('Color','w');
imagesc(phi_model_tile);
axis image;
colorbar;
title('Modeled interferometric phase, tile');

figure('Color','w');
imagesc(amp_model_tile);
axis image;
colorbar;
title('Modeled interferometric amplitude, tile');

components = {'Eas','Es','Esg','Egv'};

for ic = 1:numel(components)
    c = components{ic};

    figure('Color','w');
    imagesc(abs(out26_tile.(c)));
    axis image;
    colorbar;
    title(['Feb 26 |', c, '|']);
end

EtotAbs = abs(out26_tile.E);

frac_Es  = abs(out26_tile.Es)  ./ EtotAbs;
frac_Esg = abs(out26_tile.Esg) ./ EtotAbs;
frac_Egv = abs(out26_tile.Egv) ./ EtotAbs;
frac_Eas = abs(out26_tile.Eas) ./ EtotAbs;

frac_Es(~isfinite(frac_Es)) = NaN;
frac_Esg(~isfinite(frac_Esg)) = NaN;
frac_Egv(~isfinite(frac_Egv)) = NaN;
frac_Eas(~isfinite(frac_Eas)) = NaN;

figure('Color','w');
imagesc(frac_Egv);
axis image;
colorbar;
title('Fractional contribution: |Egv| / |E|');

%% Test 0.25% - 0.5% LWC Condition
% phi025 = phi_model_tile;
% phi050 = phi_model_tile;

dphi_050_minus_025 = angle(exp(1i * (phi050 - phi025)));

figure('Color','w');
imagesc(dphi_050_minus_025);
axis image;
colorbar;
clim([-pi pi]);
title('Phase sensitivity: Feb23 LWC 0.5% - 0.25%');

good = validPhase & isfinite(dphi_050_minus_025);

fprintf('LWC sensitivity median dphi = %.3f rad\n', ...
    median(dphi_050_minus_025(good), 'omitnan'));

fprintf('LWC sensitivity P05/P95 dphi = %.3f / %.3f rad\n', ...
    prctile(dphi_050_minus_025(good), 5), ...
    prctile(dphi_050_minus_025(good), 95));

%% Keep Trsting each component
cases = struct([]);

cases(1).name = 'all';
cases(1).snow_surface_A0 = pars1.snow_surface.A0;
cases(1).soil_surface_A0 = pars1.soil_surface.A0;
cases(1).soil_volume_enable = true;
cases(1).As = pars1.As;

cases(2).name = 'no_Eas';
cases(2).snow_surface_A0 = 0;
cases(2).soil_surface_A0 = pars1.soil_surface.A0;
cases(2).soil_volume_enable = true;
cases(2).As = pars1.As;

cases(3).name = 'Esg_only';
cases(3).snow_surface_A0 = 0;
cases(3).soil_surface_A0 = pars1.soil_surface.A0;
cases(3).soil_volume_enable = false;
cases(3).As = 0;

cases(4).name = 'Egv_only';
cases(4).snow_surface_A0 = 0;
cases(4).soil_surface_A0 = 0;
cases(4).soil_volume_enable = true;
cases(4).As = 0;

cases(5).name = 'Eas_only';
cases(5).snow_surface_A0 = pars1.snow_surface.A0;
cases(5).soil_surface_A0 = 0;
cases(5).soil_volume_enable = false;
cases(5).As = 0;

for ic = 1:numel(cases)

    p1 = pars1;
    p2 = pars2;

    p1.snow_surface.A0 = cases(ic).snow_surface_A0;
    p2.snow_surface.A0 = cases(ic).snow_surface_A0;

    p1.soil_surface.A0 = cases(ic).soil_surface_A0;
    p2.soil_surface.A0 = cases(ic).soil_surface_A0;

    p1.soil_volume.enable = cases(ic).soil_volume_enable;
    p2.soil_volume.enable = cases(ic).soil_volume_enable;

    p1.As = cases(ic).As;
    p2.As = cases(ic).As;

    out23_ab = insar_forward_profile_lut_pass_vec( ...
        theta_tile, f, HS23_tile, epsSnow23_tile, ...
        epsSoilSurf_tile, Ig_tile, p1);

    out26_ab = insar_forward_profile_lut_pass_vec( ...
        theta_tile, f, HS26_tile, epsSnow26_tile, ...
        epsSoilSurf_tile, Ig_tile, p2);

    phi_ab = angle(out26_ab.E .* conj(out23_ab.E));
    amp_ab = abs(out26_ab.E .* conj(out23_ab.E));

    good = isfinite(phi_ab) & amp_ab > prctile(amp_ab(:), 20, 'all');

    fprintf('\n%s\n', cases(ic).name);
    fprintf('  median phi = %.3f rad\n', median(phi_ab(good), 'omitnan'));
    fprintf('  P05/P95    = %.3f / %.3f rad\n', ...
        prctile(phi_ab(good), 5), prctile(phi_ab(good), 95));

    figure('Color','w');
    imagesc(phi_ab);
    axis image;
    colorbar;
    clim([-pi pi]);
    title(['Modeled phase: ', cases(ic).name]);

end
%%
% Storage containers
out23Cases = struct();
out26Cases = struct();
phiCases   = struct();
ampCases   = struct();
goodCases  = struct();

for ic = 1:numel(cases)

    p1 = pars1;
    p2 = pars2;

    p1.snow_surface.A0 = cases(ic).snow_surface_A0;
    p2.snow_surface.A0 = cases(ic).snow_surface_A0;

    p1.soil_surface.A0 = cases(ic).soil_surface_A0;
    p2.soil_surface.A0 = cases(ic).soil_surface_A0;

    p1.soil_volume.enable = cases(ic).soil_volume_enable;
    p2.soil_volume.enable = cases(ic).soil_volume_enable;

    p1.As = cases(ic).As;
    p2.As = cases(ic).As;

    out23_ab = insar_forward_profile_lut_pass_vec( ...
        theta_tile, f, HS23_tile, epsSnow23_tile, ...
        epsSoilSurf_tile, Ig_tile, p1);

    out26_ab = insar_forward_profile_lut_pass_vec( ...
        theta_tile, f, HS26_tile, epsSnow26_tile, ...
        epsSoilSurf_tile, Ig_tile, p2);

    phi_ab = angle(out26_ab.E .* conj(out23_ab.E));
    amp_ab = abs(out26_ab.E .* conj(out23_ab.E));

    % Robust amplitude threshold
    ampVals = amp_ab(isfinite(amp_ab));
    ampThresh = prctile(ampVals, 20);

    good = isfinite(phi_ab) & isfinite(amp_ab) & amp_ab > ampThresh;

    % Make case name safe as a struct field
    caseName = matlab.lang.makeValidName(cases(ic).name);

    % Store outputs
    out23Cases.(caseName) = out23_ab;
    out26Cases.(caseName) = out26_ab;
    phiCases.(caseName)   = phi_ab;
    ampCases.(caseName)   = amp_ab;
    goodCases.(caseName)  = good;

    fprintf('\n%s\n', cases(ic).name);
    fprintf('  median phi = %.3f rad\n', median(phi_ab(good), 'omitnan'));
    fprintf('  P05/P95    = %.3f / %.3f rad\n', ...
        prctile(phi_ab(good), 5), prctile(phi_ab(good), 95));

    figure('Color','w');
    imagesc(phi_ab);
    axis image;
    colorbar;
    clim([-pi pi]);
    title(['Modeled phase: ', cases(ic).name]);

end

% Phase Separation Diagnostics
phi_Eas = phiCases.Eas_only;
phi_Esg = phiCases.Esg_only;
phi_Egv = phiCases.Egv_only;
phi_all = phiCases.all;

good = goodCases.all & ...
       isfinite(phi_Eas) & isfinite(phi_Esg) & isfinite(phi_Egv);

dphi_Eas_Esg = angle(exp(1i .* (phi_Eas - phi_Esg)));
dphi_Eas_Egv = angle(exp(1i .* (phi_Eas - phi_Egv)));
dphi_Esg_Egv = angle(exp(1i .* (phi_Esg - phi_Egv)));

fprintf('\nPhase separation diagnostics\n');

fprintf('Eas - Esg median = %.3f rad\n', ...
    median(dphi_Eas_Esg(good), 'omitnan'));

fprintf('Eas - Esg P05/P95 = %.3f / %.3f rad\n', ...
    prctile(dphi_Eas_Esg(good), 5), ...
    prctile(dphi_Eas_Esg(good), 95));

fprintf('Esg - Egv median = %.3f rad\n', ...
    median(dphi_Esg_Egv(good), 'omitnan'));

fprintf('Esg - Egv P05/P95 = %.3f / %.3f rad\n', ...
    prctile(dphi_Esg_Egv(good), 5), ...
    prctile(dphi_Esg_Egv(good), 95));

figure('Color','w');
imagesc(dphi_Eas_Esg);
axis image;
colorbar;
clim([-pi pi]);
title('Phase separation: Eas - Esg');

figure('Color','w');
imagesc(dphi_Esg_Egv);
axis image;
colorbar;
clim([-pi pi]);
title('Phase separation: Esg - Egv');

% Coherent Cancelation Index
outAll26 = out26Cases.all;

E_sum_abs = abs(outAll26.E);

sum_abs_components = ...
    abs(outAll26.Eas) + ...
    abs(outAll26.Es)  + ...
    abs(outAll26.Esg) + ...
    abs(outAll26.Egv);

cancelIndex = E_sum_abs ./ sum_abs_components;
cancelIndex(~isfinite(cancelIndex)) = NaN;

figure('Color','w');
imagesc(cancelIndex);
axis image;
colorbar;
clim([0 1]);
title('|sum(E_i)| / sum(|E_i|): coherent cancellation index');

% Bounded Copnent Fractions
den = sum_abs_components;

f_Eas_norm = abs(outAll26.Eas) ./ den;
f_Es_norm  = abs(outAll26.Es)  ./ den;
f_Esg_norm = abs(outAll26.Esg) ./ den;
f_Egv_norm = abs(outAll26.Egv) ./ den;

f_Eas_norm(~isfinite(f_Eas_norm)) = NaN;
f_Es_norm(~isfinite(f_Es_norm))   = NaN;
f_Esg_norm(~isfinite(f_Esg_norm)) = NaN;
f_Egv_norm(~isfinite(f_Egv_norm)) = NaN;

figure('Color','w'); imagesc(f_Eas_norm); axis image; colorbar; clim([0 1]); title('|Eas| / sum(|E_i|)');
figure('Color','w'); imagesc(f_Esg_norm); axis image; colorbar; clim([0 1]); title('|Esg| / sum(|E_i|)');
figure('Color','w'); imagesc(f_Egv_norm); axis image; colorbar; clim([0 1]); title('|Egv| / sum(|E_i|)');
clim([quantile(f_Egv_norm(:),[0.01,0.99])])

Esub23 = out23Cases.all.Esg + out23Cases.all.Egv;
Esub26 = out26Cases.all.Esg + out26Cases.all.Egv;

phi_sub = angle(Esub26 .* conj(Esub23));
amp_sub = abs(Esub26 .* conj(Esub23));

phi_surface = phiCases.Eas_only;
phi_all     = phiCases.all;

dphi_surface_sub = angle(exp(1i .* (phi_surface - phi_sub)));

figure('Color','w');
imagesc(dphi_surface_sub);
axis image;
colorbar;
clim([-pi pi]);
title('Phase separation: surface - subsurface');

figure('Color','w');
imagesc(phi_sub);
axis image;
colorbar;
clim([-pi pi]);
title('Modeled phase: Esg + Egv');

A_surface = abs(out26Cases.all.Eas);
A_sub     = abs(out26Cases.all.Esg + out26Cases.all.Egv);

ratio_surface_to_sub = A_surface ./ A_sub;
ratio_surface_to_sub(~isfinite(ratio_surface_to_sub)) = NaN;

figure('Color','w');
imagesc(log10(ratio_surface_to_sub));
axis image;
colorbar;
title('log10(|Eas| / |Esg + Egv|)');

% In much of the region with local incidence less than roughly 70 degrees, 
% the model predicts that the snow/soil subsurface return is stronger than 
% the air/snow free-surface return, but the dominance is spatially 
% modulated by terrain, snow state, and interface amplitude.

%% Compare to UAS InSAR Data
phi_uas = phz;
phi_model_all = phiCases.all;
phi_model_Eas = phiCases.Eas_only;
phi_model_noEas = phiCases.no_Eas;
res_all    = angle(exp(1i .* (phi_uas - phi_model_all)));
res_Eas    = angle(exp(1i .* (phi_uas - phi_model_Eas)));
res_noEas  = angle(exp(1i .* (phi_uas - phi_model_noEas)));

valid = isfinite(phi_uas) & isfinite(phi_model_all) & ...
        isfinite(phi_model_Eas) & isfinite(phi_model_noEas);

models = struct();
models.all   = res_all;
models.Eas   = res_Eas;
models.noEas = res_noEas;

names = fieldnames(models);

for ii = 1:numel(names)

    name = names{ii};
    r = models.(name);

    rv = r(valid);

    circMean = angle(mean(exp(1i .* rv), 'omitnan'));
    circR    = abs(mean(exp(1i .* rv), 'omitnan'));
    circStd  = sqrt(-2 .* log(circR));

    maeWrapped = median(abs(rv), 'omitnan');

    fprintf('\n%s residual\n', name);
    fprintf('  circular mean residual = %.3f rad\n', circMean);
    fprintf('  circular concentration R = %.3f\n', circR);
    fprintf('  circular std = %.3f rad\n', circStd);
    fprintf('  median |wrapped residual| = %.3f rad\n', maeWrapped);

end

figure();
subplot(1,3,1)
imagesc(res_all)
colorbar; colormap(cmap)
daspect([1,1,1])
title('Phase Residual All Components')
subplot(1,3,2)
imagesc(res_Eas)
colorbar; colormap(cmap)
daspect([1,1,1])
title('Phase Residual Eas Component Only')

subplot(1,3,3)
imagesc(res_noEas)
colorbar; colormap(cmap)
daspect([1,1,1])
title('Phase Residual No Eas Component')

%% TryEas weighting
Esub23 = out23Cases.all.Esg + out23Cases.all.Egv;
Esub26 = out26Cases.all.Esg + out26Cases.all.Egv;

Eas23 = out23Cases.all.Eas;
Eas26 = out26Cases.all.Eas;

aVals = linspace(0, 1, 51);

valid = isfinite(phi_uas) & ...
    isfinite(coh) & coh >= 0.75 & ...
        isfinite(Esub23) & isfinite(Esub26) & ...
        isfinite(Eas23) & isfinite(Eas26);
valid = referenceMask;

for ia = 1:numel(aVals)

    a = aVals(ia);

    E23 = Esub23 + a .* Eas23;
    E26 = Esub26 + a .* Eas26;

    phi_model = angle(E26 .* conj(E23));
    phi_model = insar.vector_median_filter(phi_model,11);
    res = angle(exp(1i .* (insar.vector_median_filter(phi_uas,11) - phi_model)));

    rv = res(valid);

    Rr(ia) = abs(mean(exp(1i .* rv), 'omitnan'));
    circStd(ia) = sqrt(-2 .* log(Rr(ia)));
    medAbs(ia) = median(abs(rv), 'omitnan');

end

[~, iBestR] = max(Rr);
[~, iBestMed] = min(medAbs);

fprintf('Best Eas scale by R: %.3f, R = %.3f, circStd = %.3f\n', ...
    aVals(iBestR), Rr(iBestR), circStd(iBestR));

fprintf('Best Eas scale by median abs residual: %.3f, medAbs = %.3f\n', ...
    aVals(iBestMed), medAbs(iBestMed));

figure('Color','w');
plot(aVals, Rr, 'k', 'LineWidth', 2);
xlabel('Eas amplitude scale');
ylabel('Residual circular concentration R');
grid on;

figure('Color','w');
plot(aVals, medAbs, 'k', 'LineWidth', 2);
xlabel('Eas amplitude scale');
ylabel('Median |wrapped residual| (rad)');
grid on;

%Notes:
% The physically modeled air/snow free-surface phase is real and 
% can strongly affect the coherent sum, but comparison with the high-coherence 
% UAS phase suggests that the effective coherent air/snow contribution 
% is much weaker than the raw full-strength model. The observed phase is 
% primarily subsurface/interface controlled, with possible weak surface-return mixing.
% After applying the vector median filter the subsurface term controls
% the broad phase structure, while the air/snow surface term contributes non-negligibly once speckle is suppressed.
%% Median filtering
win = 11;
aVals = linspace(0, 1, 101);

[~, z_uas_filt, R_uas_filt] = insar.vector_median_filter_phase(phi_uas, win);

for ia = 1:numel(aVals)

    a = aVals(ia);

    E23 = Esub23 + a .* Eas23;
    E26 = Esub26 + a .* Eas26;

    phi_model = angle(E26 .* conj(E23));The result after phase cons

    [~, z_model_filt, R_model_filt] = insar.vector_median_filter_phase(phi_model, win);

    res = angle(z_uas_filt .* conj(z_model_filt));

    valid = isfinite(res) & isfinite(coh) & coh > 0.75 & ...
            R_uas_filt > 0.2 & R_model_filt > 0.2;

    rv = res(valid);

    Rmetric(ia) = abs(mean(exp(1i .* rv), 'omitnan'));
    medAbs(ia) = median(abs(rv), 'omitnan');

end

[~, iBestR] = max(Rmetric);
[~, iBestMed] = min(medAbs);

fprintf('Best Eas scale by R: %.3f, R = %.3f, circStd = %.3f\n', ...
    aVals(iBestR), Rmetric(iBestR), circStd(iBestR));

fprintf('Best Eas scale by median abs residual: %.3f, medAbs = %.3f\n', ...
    aVals(iBestMed), medAbs(iBestMed));

figure('Color','w');
plot(aVals, Rmetric, 'k', 'LineWidth', 2);
xlabel('Eas amplitude scale');
ylabel('Residual circular concentration R');
grid on;

figure('Color','w');
plot(aVals, medAbs, 'k', 'LineWidth', 2);
xlabel('Eas amplitude scale');
ylabel('Median |wrapped residual| (rad)');
grid on;

%% True vector median filtering
win = 11;

[phi_uas_vmf, z_uas_vmf, score_uas, R_uas] = ...
    insar.true_vector_median_filter_phase(phi_uas, win, ...
    'Normalize', true, ...
    'UseParallel', true);

[phi_model_vmf, z_model_vmf, score_model, R_model] = ...
    insar.true_vector_median_filter_phase(phiCases.all, win, ...
    'Normalize', true, ...
    'UseParallel', true);

[phi_modelNoEas_vmf, z_modelNoEas_vmf, score_model_noEas, R_model_noEas] = ...
    insar.true_vector_median_filter_phase(phiCases.no_Eas, win, ...
    'Normalize', true, ...
    'UseParallel', true);

%% Phase referencing Bias Correction
referenceMask = mapinterp(insarData.phaseFillTrustedMask,demData.R,Coords.utmX,Coords.utmY);
wrap = @(x) angle(exp(1i .* x));
res0 = wrap(phi_uas_vmf - phi_model_vmf);
bestBias = angle(mean(exp(1i .* res0(referenceMask)), 'omitnan'));
phiModel_bestBias = wrap(phi_model_vmf + bestBias);

res0noEas = wrap(phi_uas_vmf - phi_modelNoEas_vmf);
bestBiasNoEas = angle(mean(exp(1i .* res0noEas(referenceMask)), 'omitnan'));
phiModel_bestBiasNoEas = wrap(phi_modelNoEas_vmf + bestBiasNoEas);

% All-components model
resAll_bias = wrap(phi_uas_vmf - phiModel_bestBias);

% No-Eas model
resNoEas_bias = wrap(phi_uas_vmf - phiModel_bestBiasNoEas);

validAll = referenceMask & isfinite(resAll_bias);
validNoEas = referenceMask & isfinite(resNoEas_bias);

R_all = abs(mean(exp(1i .* resAll_bias(validAll)), 'omitnan'));
R_noEas = abs(mean(exp(1i .* resNoEas_bias(validNoEas)), 'omitnan'));

med_all = median(abs(resAll_bias(validAll)), 'omitnan');
med_noEas = median(abs(resNoEas_bias(validNoEas)), 'omitnan');

fprintf('\nBias-corrected residuals over referenceMask\n');

fprintf('All components:\n');
fprintf('  bias = %.4f rad\n', bestBias);
fprintf('  R = %.3f\n', R_all);
fprintf('  median |res| = %.3f rad\n', med_all);

fprintf('No Eas:\n');
fprintf('  bias = %.4f rad\n', bestBiasNoEas);
fprintf('  R = %.3f\n', R_noEas);
fprintf('  median |res| = %.3f rad\n', med_noEas);
%% Export InSAR forward-model diagnostic outputs
% Assumes these exist:
%   out23Cases, out26Cases, phiCases, ampCases, goodCases
%   Coords.R, Coords.EPSG
%   pars1, pars2
%   targetMeanLWC23, ii   % optional, for naming

% -------------------------
% Output directory
% -------------------------
baseOutDir = 'D:\IDALS\2026\forward_model_runs\';

if exist('targetMeanLWC23','var') && exist('ii','var')
    lwcLabel = strrep(sprintf('LWC23mean_%0.4f', targetMeanLWC23(ii)), '.', 'p');
else
    lwcLabel = 'LWC23_currentCase';
end

runTag = ['20260223_20260226_fixedSoil_surfaceHeight_', lwcLabel];

outDir = fullfile(baseOutDir, runTag);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fprintf('Exporting forward-model outputs to:\n%s\n', outDir);

% -------------------------
% Reconstruct key phasors
% -------------------------
Eas23 = out23Cases.all.Eas;
Eas26 = out26Cases.all.Eas;

Esub23 = out23Cases.all.Esg + out23Cases.all.Egv;
Esub26 = out26Cases.all.Esg + out26Cases.all.Egv;

Eall23 = out23Cases.all.E;
Eall26 = out26Cases.all.E;

% -------------------------
% Key modeled phases
% -------------------------
phi_all = angle(Eall26 .* conj(Eall23));
phi_surface = angle(Eas26 .* conj(Eas23));
phi_subsurface = angle(Esub26 .* conj(Esub23));

% If already stored, these should match:
% phi_all       = phiCases.all;
% phi_surface   = phiCases.Eas_only;

% -------------------------
% Amplitude diagnostics
% -------------------------
amp_all = abs(Eall26 .* conj(Eall23));
amp_surface = abs(Eas26 .* conj(Eas23));
amp_subsurface = abs(Esub26 .* conj(Esub23));

% Two-pass geometric mean amplitude ratio:
% surface/subsurface return dominance
A_surface_ifg = sqrt(abs(Eas23) .* abs(Eas26));
A_sub_ifg     = sqrt(abs(Esub23) .* abs(Esub26));

ratioLog_ifg = log10(A_surface_ifg ./ A_sub_ifg);
ratioLog_ifg(~isfinite(ratioLog_ifg)) = NaN;

% Component phase separation
dphi_surface_subsurface = angle(exp(1i .* (phi_surface - phi_subsurface)));

% Coherent cancellation index, pass 26
sumAbs26 = ...
    abs(out26Cases.all.Eas) + ...
    abs(out26Cases.all.Es)  + ...
    abs(out26Cases.all.Esg) + ...
    abs(out26Cases.all.Egv);

cancelIndex26 = abs(out26Cases.all.E) ./ sumAbs26;
cancelIndex26(~isfinite(cancelIndex26)) = NaN;

% Bounded component fractions, pass 26
f_Eas_26 = abs(out26Cases.all.Eas) ./ sumAbs26;
f_Es_26  = abs(out26Cases.all.Es)  ./ sumAbs26;
f_Esg_26 = abs(out26Cases.all.Esg) ./ sumAbs26;
f_Egv_26 = abs(out26Cases.all.Egv) ./ sumAbs26;

f_Eas_26(~isfinite(f_Eas_26)) = NaN;
f_Es_26(~isfinite(f_Es_26))   = NaN;
f_Esg_26(~isfinite(f_Esg_26)) = NaN;
f_Egv_26(~isfinite(f_Egv_26)) = NaN;

% -------------------------
% Return-regime classification
% -------------------------
% ratioLog_ifg < -0.5 : subsurface > ~3.16x surface
% ratioLog_ifg >  0.5 : surface > ~3.16x subsurface
validRatio = isfinite(ratioLog_ifg);

returnRegime = zeros(size(ratioLog_ifg), 'single');
returnRegime(validRatio & ratioLog_ifg < -0.5) = 1;  % subsurface dominated
returnRegime(validRatio & abs(ratioLog_ifg) <= 0.5) = 2;  % mixed
returnRegime(validRatio & ratioLog_ifg > 0.5) = 3;   % surface dominated
returnRegime(~validRatio) = NaN;

fprintf('\nReturn-regime summary:\n');
fprintf('  Subsurface dominated: %.2f %%\n', ...
    100 * nnz(returnRegime == 1) / nnz(validRatio));
fprintf('  Mixed return:         %.2f %%\n', ...
    100 * nnz(returnRegime == 2) / nnz(validRatio));
fprintf('  Surface dominated:    %.2f %%\n', ...
    100 * nnz(returnRegime == 3) / nnz(validRatio));

% -------------------------
% Export GeoTIFFs
% -------------------------
write_tif(single(phi_all), ...
    fullfile(outDir, [runTag, '_phi_all.tif']), Coords);

write_tif(single(phi_surface), ...
    fullfile(outDir, [runTag, '_phi_surface_Eas.tif']), Coords);

write_tif(single(phi_subsurface), ...
    fullfile(outDir, [runTag, '_phi_subsurface_EsgEgv.tif']), Coords);

write_tif(single(dphi_surface_subsurface), ...
    fullfile(outDir, [runTag, '_dphi_surface_minus_subsurface.tif']), Coords);

write_tif(single(amp_all), ...
    fullfile(outDir, [runTag, '_amp_all.tif']), Coords);

write_tif(single(amp_surface), ...
    fullfile(outDir, [runTag, '_amp_surface_Eas.tif']), Coords);

write_tif(single(amp_subsurface), ...
    fullfile(outDir, [runTag, '_amp_subsurface_EsgEgv.tif']), Coords);

write_tif(single(ratioLog_ifg), ...
    fullfile(outDir, [runTag, '_log10_Eas_over_EsgEgv.tif']), Coords);

write_tif(single(cancelIndex26), ...
    fullfile(outDir, [runTag, '_cancelIndex_pass26.tif']), Coords);

write_tif(single(f_Eas_26), ...
    fullfile(outDir, [runTag, '_frac_Eas_pass26.tif']), Coords);

write_tif(single(f_Esg_26), ...
    fullfile(outDir, [runTag, '_frac_Esg_pass26.tif']), Coords);

write_tif(single(f_Egv_26), ...
    fullfile(outDir, [runTag, '_frac_Egv_pass26.tif']), Coords);

write_tif(single(returnRegime), ...
    fullfile(outDir, [runTag, '_returnRegime_1sub_2mixed_3surf.tif']), Coords);

% Optional: export individual ablation phases if they exist
caseNames = fieldnames(phiCases);

for ic = 1:numel(caseNames)

    cName = caseNames{ic};

    write_tif(single(phiCases.(cName)), ...
        fullfile(outDir, [runTag, '_phi_', cName, '.tif']), Coords);

    if isfield(ampCases, cName)
        write_tif(single(ampCases.(cName)), ...
            fullfile(outDir, [runTag, '_amp_', cName, '.tif']), Coords);
    end
end

% -------------------------
% Save compact MAT package
% -------------------------
runInfo = struct();

runInfo.runTag = runTag;
runInfo.outDir = outDir;
runInfo.created = datetime('now');

if exist('targetMeanLWC23','var') && exist('ii','var')
    runInfo.LWC23_targetMean = targetMeanLWC23(ii);
    runInfo.LWC23_caseIndex = ii;
end

runInfo.description = ...
    ['Fixed-soil baseline; GPR-constrained VWC0 soil state; ', ...
     'snow surface height phase enabled; products include all, surface, ', ...
     'subsurface, amplitude ratio, cancellation, and return regime maps.'];

runInfo.returnRegimeLegend = struct( ...
    'subsurfaceDominated', 1, ...
    'mixedReturn', 2, ...
    'surfaceDominated', 3);

runInfo.ratioLogDefinition = ...
    'log10(sqrt(|Eas23|*|Eas26|) / sqrt(|Esub23|*|Esub26|)), Esub = Esg + Egv';

runInfo.phaseDefinition = ...
    'phi = angle(E26 .* conj(E23))';

runInfo.pars1 = pars1;
runInfo.pars2 = pars2;

if exist('soilFwdLUT','var')
    runInfo.soilFwdLUT_file_note = ...
        'soilFwdLUT exists in workspace; save separately if not already saved.';
end

quickProducts = struct();
quickProducts.phi_all = single(phi_all);
quickProducts.phi_surface = single(phi_surface);
quickProducts.phi_subsurface = single(phi_subsurface);
quickProducts.dphi_surface_subsurface = single(dphi_surface_subsurface);
quickProducts.ratioLog_ifg = single(ratioLog_ifg);
quickProducts.cancelIndex26 = single(cancelIndex26);
quickProducts.returnRegime = single(returnRegime);

save(fullfile(outDir, [runTag, '_quickReload.mat']), ...
    'runInfo', ...
    'quickProducts', ...
    '-v7.3');

fprintf('\nExport complete.\nSaved quick reload file:\n%s\n', ...
    fullfile(outDir, [runTag, '_quickReload.mat']));
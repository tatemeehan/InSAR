%% Coherence Hole Analysis
% Inputs from your pair (i,b1)-(j,b2) via geomData.sarGeometry{idx}
B = geomData.sarGeometry{geomData.slcGeomIndex.idx}.Bperp;      % m
r = geomData.sarGeometry{geomData.slcGeomIndex.idx}.slant;      % m
th = deg2rad(geomData.sarGeometry{geomData.slcGeomIndex.idx}.incidence); % rad
kz = (4*pi/params.lambda) .* (B ./ (r .* max(sin(th),1e-3)));
C = insarData.coherence;
mask = isfinite(C);

rangeIndex1 = geomData.sarGeometry{geomData.slcGeomIndex(1,1).idx}.closestIndex_master;
rangeIndex2 = geomData.sarGeometry{geomData.slcGeomIndex(1,1).idx}.closestIndex_slave;
slant1 = geomData.sarGeometry{geomData.slcGeomIndex(1,1).idx}.slant;
slant2 = geomData.sarGeometry{geomData.slcGeomIndex(1,1).idx}.slant2;



% Inputs you already have:
% C  = insarData(pairIdx).coherence;       % measured coherence
% kz = kz_map;                              % from your geometry (rad/m)

% Optional SNR correction (recommended)
% gammaSNR = [];%your_gamma_snr_map;            % or [] to skip
% [g1SNR, NES1] = insar.snr_from_gamma0_empirical(sarData(1).pow{1}, mask);
% [g2SNR, NES2] = insar.snr_from_gamma0_empirical(sarData(1).pow{2}, mask);


% [g1SNRs, NES1s] = insar.snr_from_gamma0_by_slant(sarData(1).pow{1}, slant1, mask);
% [g2SNRs, NES2s] = insar.snr_from_gamma0_by_slant(sarData(1).pow{2}, slant2, mask);
% 
% 
% gammaSNR = [];%your_gamma_snr_map;            % or [] to skip

% [g1SNR, NES1, prof] = insar.snr_from_gamma0_by_slant2D(sarData(1).pow{1}, slant1, rangeIndex1, mask);
% [g1SNR, NES1,diagOut] = insar.snr_from_gamma0_knn(sarData(1).pow{1}, slant1, rangeIndex1, mask);
% KNN Search Fast via Decimation
tic
[g1SNR, NES1,diagOut1] = insar.snr_from_gamma0_knn_fast(sarData(1).pow{1}, slant1, rangeIndex1, mask);
toc
tic
[g2SNR, NES2,diagOut2] = insar.snr_from_gamma0_knn_fast(sarData(1).pow{2}, slant2, rangeIndex2, mask);
toc

gammaSNR = insar.coherence_ceiling_from_snr(g1SNR, g2SNR);


% --- Single component ---
% res1 = insar.fit_sigma_z_single(C, kz, struct('mask', mask,'gammaSNR',gammaSNR,'robust',true,'doplot',true));
% fprintf('Single-comp sigma_z = %.3f m (R^2=%.2f, n=%d)\n', res1.sigma_z, res1.r2, res1.n);
opts.gamma_snr = gammaSNR;
out = insar.fit_sigma_gridsearch_single(C, kz, mask, opts);
nBoot = 100;
[sigma_hat, sigma_ci] = insar.sigma1_bootstrap(kz, C, gammaSNR, nBoot);

% --- Two-component RVoG (magnitude) ---
% res2 = insar.fit_sigma_z_rvog(C, kz, struct( ...
%     'gammaSNR',gammaSNR, 'robust',true, 'plot',true, ...
%     'bounds',struct('sigma_z',[0 2],'alpha',[0 1],'z0',[-1 1]) ));
% fprintf('RVoG: sigma_z=%.3f m, alpha=%.2f, z0=%.2f m\n', res2.sigma_z, res2.alpha, res2.z0);

% Single Parameter Model is Sufficient!!
% Two Component Fit Struggles even with Tikhonov Regularization and Clamped
% Parameter Space
% out2 = insar.fit_sigma_gridsearch_two(C, kz, mask, opts);

% Veritcal Roughness Model
coherenceMod1 = exp(-0.5.*(out.sigma_z.*kz).^2);
residualC = (coherenceMod1(:)-C(:));
nanmedian(residualC)
figure();
subplot(1,3,1)
binscatter(kz(:),residualC(:))
xlim([0 10])
xlabel('Vertical Scattering');ylabel('Coherence Residual')
subplot(1,3,2)
binscatter(kz(:),nanmean([slant1(:),slant2(:)],2))
xlabel('Vertical Scattering');ylabel('Slant Range')
xlim([0,10])
ylim([0 800])
subplot(1,3,3)
binscatter(kz(:),rad2deg(th(:)))
xlabel('Vertical Scattering');ylabel('Slant Range')
xlim([0,10])

% coherenceMod2 = abs(out2.frac.*exp(-0.5.*(out2.sigma1.*kz).^2)+(1-out2.frac).*exp(-0.5.*(out2.sigma2.*kz).^2));

sigmaZ = out.sigma_z; % m  <-- try a few values
gamma_vol = exp(-0.5*(kz.*sigmaZ).^2);


figure(); binscatter(gamma_vol(mask), C(mask)); grid on
xlabel(['\gamma_{vol} model \sigma = ',num2str(sigmaZ)]); ylabel('Measured coherence');
title('Does vertical-wavenumber decorrelation explain the drop?');

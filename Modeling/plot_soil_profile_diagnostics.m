function plot_soil_profile_diagnostics(diag, f_Hz)
%PLOT_SOIL_PROFILE_DIAGNOSTICS Visualize soil profile effects on Ig.
% Adds frozen-model diagnostics if present:
%   diag.theta_w_liq, diag.theta_w_ice, diag.fL, (optional) diag.sigma

z = diag.z;

% ---------------- Figure 1: Profiles ----------------
figure;
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(diag.VWC, z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
xlabel('VWC (m^3/m^3)'); ylabel('Depth z (m)'); title('VWC(z)');
ylim([0,max(diag.z)])


nexttile;
plot(diag.Tc, z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
xlabel('T (^{\circ}C)'); ylabel('Depth z (m)'); title('T(z)');
ylim([0,max(diag.z)])

nexttile;
plot(real(diag.eps), z, 'LineWidth', 1.5); hold on;
plot(imag(diag.eps), z, 'LineWidth', 1.5);
grid on; set(gca,'YDir','reverse');
xlabel('\epsilon'); ylabel('Depth z (m)');
legend('\epsilon''','\epsilon''''','Location','best'); title('\epsilon(z)');
ylim([0,max(diag.z)])

nexttile;
plot(diag.alpha, z, 'LineWidth', 1.5); hold on;
plot(diag.beta,  z, 'LineWidth', 1.5);
grid on; set(gca,'YDir','reverse');
xlabel('\alpha (1/m), \beta (rad/m)'); ylabel('Depth z (m)');
legend('\alpha','\beta','Location','best'); title('k_z(z) components');
ylim([0,max(diag.z)])


% ---------------- Figure 2: Frozen diagnostics (if present) ----------------
hasFrozen = isfield(diag,'theta_w_liq') && isfield(diag,'theta_w_ice') && isfield(diag,'fL');

if hasFrozen
    figure;
    tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

    % (A) Liquid vs Ice volumetric contents
    nexttile;
    plot(diag.theta_w_liq, z, 'LineWidth', 1.5); hold on;
    plot(diag.theta_w_ice, z, 'LineWidth', 1.5);
    grid on; set(gca,'YDir','reverse');
    xlabel('\theta_w (m^3/m^3)'); ylabel('Depth z (m)');
    legend('\theta_{w,liq}','\theta_{w,ice}','Location','best');
    title('Liquid vs Ice water content');

    % (B) Unfrozen fraction
    nexttile;
    plot(diag.fL, z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
    xlabel('f_L'); ylabel('Depth z (m)');
    title('Unfrozen water fraction f_L(z)');
    xlim([0 1]);

    % (C) Conductivity vs depth (or sigma_eff if you later add it)
    nexttile;
    if isfield(diag,'sigma')
        plot(diag.sigma, z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
        xlabel('\sigma (S/m)'); ylabel('Depth z (m)');
        title('\sigma(z) used in model');
        xlim([0 0.1])
    else
        text(0.1,0.5,'diag.sigma not found','Units','normalized');
        axis off;
    end
end

% ---------------- Figure 3: cumulative attenuation/phase ----------------
figure;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(diag.A_cum, z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
xlabel('Two-way attenuation exponent 2\int_0^z \alpha dz'); ylabel('Depth z (m)');
title('Cumulative attenuation');

nexttile;
plot(diag.B_cum, z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
xlabel('Two-way phase 2\int_0^z \beta dz (rad)'); ylabel('Depth z (m)');
title('Cumulative phase');

% ---------------- Figure 4: integrand magnitude/phase ----------------
figure;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile;
plot(abs(diag.integrand), z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
xlabel('|integrand|'); ylabel('Depth z (m)');
title('Contribution magnitude vs depth');

nexttile;
plot(unwrap(angle(diag.integrand)), z, 'LineWidth', 1.5); grid on; set(gca,'YDir','reverse');
xlabel('arg(integrand) unwrapped (rad)'); ylabel('Depth z (m)');
title('Contribution phase vs depth');

% ---------------- Print summary ----------------
fprintf('f = %.3f GHz\n', f_Hz/1e9);
fprintf('Ig = %.6g%+.6gi  (|Ig|=%.4g)\n', real(diag.Ig), imag(diag.Ig), abs(diag.Ig));
fprintf('Dg_eff = %.6g%+.6gi\n', real(diag.Dg_eff), imag(diag.Dg_eff));
if isfield(diag,'sigma_eff')
    fprintf('sigma_eff (energy-weighted) = %.6g S/m\n', diag.sigma_eff);
end
end

% function plot_soil_profile_diagnostics(diag, f_Hz)
% %PLOT_SOIL_PROFILE_DIAGNOSTICS Visualize soil profile effects on Ig.
% 
% z = diag.z;
% 
% figure;
% subplot(2,2,1);
% plot(diag.VWC, z); grid on; set(gca,'YDir','reverse');
% xlabel('VWC (m^3/m^3)'); ylabel('Depth z (m)'); title('VWC(z)');
% 
% subplot(2,2,2);
% plot(diag.Tc, z); grid on; set(gca,'YDir','reverse');
% xlabel('T (^{\circ}C)'); ylabel('Depth z (m)'); title('T(z)');
% 
% subplot(2,2,3);
% plot(real(diag.eps), z); hold on; plot(imag(diag.eps), z);
% grid on; set(gca,'YDir','reverse');
% xlabel('\epsilon'); ylabel('Depth z (m)');
% legend('\epsilon''','\epsilon'''''); title('\epsilon(z)');
% 
% subplot(2,2,4);
% plot(diag.alpha, z); hold on; plot(diag.beta, z);
% grid on; set(gca,'YDir','reverse');
% xlabel('\alpha (1/m), \beta (rad/m)'); ylabel('Depth z (m)');
% legend('\alpha','\beta'); title('k_z(z) components');
% 
% figure;
% subplot(1,2,1);
% plot(diag.A_cum, z); grid on; set(gca,'YDir','reverse');
% xlabel('Two-way attenuation exponent 2\int_0^z \alpha dz'); ylabel('Depth z (m)');
% title('Cumulative attenuation');
% 
% subplot(1,2,2);
% plot(diag.B_cum, z); grid on; set(gca,'YDir','reverse');
% xlabel('Two-way phase 2\int_0^z \beta dz (rad)'); ylabel('Depth z (m)');
% title('Cumulative phase');
% 
% figure;
% subplot(1,2,1);
% plot(abs(diag.integrand), z); grid on; set(gca,'YDir','reverse');
% xlabel('|integrand|'); ylabel('Depth z (m)');
% title('Contribution magnitude vs depth');
% 
% subplot(1,2,2);
% plot(unwrap(angle(diag.integrand)), z); grid on; set(gca,'YDir','reverse');
% xlabel('arg(integrand) unwrapped (rad)'); ylabel('Depth z (m)');
% title('Contribution phase vs depth');
% 
% fprintf('f = %.3f GHz\n', f_Hz/1e9);
% fprintf('Ig = %.6g%+.6gi  (|Ig|=%.4g)\n', real(diag.Ig), imag(diag.Ig), abs(diag.Ig));
% fprintf('Dg_eff = %.6g%+.6gi\n', real(diag.Dg_eff), imag(diag.Dg_eff));
% end

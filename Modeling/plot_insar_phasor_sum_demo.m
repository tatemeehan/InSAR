function plot_insar_phasor_sum_demo()
%PLOT_INSAR_PHASOR_SUM_DEMO
% Conceptual phasor-addition figure for the 4-component snow-soil model:
%   Ebar = Eas + Es + Esg + Egv
%
% Poster-style, grayscale, head-to-tail vector addition.

close all;

%% --- Define component phasors (conceptual example) ---
% You can replace these with real modeled values later.
Eas = 0.75 * exp(1i*deg2rad(82));   % air-snow surface
Es  = 0.42 * exp(1i*deg2rad(28));   % snow volume
Esg = 0.58 * exp(1i*deg2rad(3));    % snow-soil interface
Egv = 0.30 * exp(1i*deg2rad(-32));  % soil volume

% Partial sums
S1 = Eas;
S2 = Eas + Es;
S3 = Eas + Es + Esg;
Ebar = Eas + Es + Esg + Egv;

%% --- Styling ---
colEas = [0.15 0.15 0.15];
colEs  = [0.45 0.45 0.45];
colEsg = [0.55 0.55 0.55];
colEgv = [0.35 0.35 0.35];
colRes = [0.00 0.00 0.00];
colAux = [0.70 0.70 0.70];

lwComp = 3.0;
lwRes  = 4.5;
lwAux  = 1.5;

%% --- Figure setup ---
figure('Color','w','Position',[100 100 950 720]);
ax = axes;
hold(ax,'on');
axis(ax,'equal');
box(ax,'on');

% Limits chosen to leave room for labels
xlim([-0.2 2.2]);
ylim([-0.2 1.8]);

% Light grid
xticks(0:0.4:2.0);
yticks(0:0.4:1.6);
grid(ax,'on');
ax.GridColor = [0.82 0.82 0.82];
ax.GridAlpha = 0.9;
ax.LineWidth = 1.2;
ax.FontName = 'Times New Roman';
ax.FontSize = 16;

xlabel('Real', 'FontWeight','bold', 'FontSize',20);
ylabel('Imaginary', 'FontWeight','bold', 'FontSize',20);
title('Coherent Phasor Summation', 'FontWeight','bold', 'FontSize',28);

%% --- Origin and partial-sum points ---
O  = 0 + 0i;
P1 = S1;
P2 = S2;
P3 = S3;
P4 = Ebar;

%% --- Optional faint partial-sum guides from origin ---
plot([0 real(P1)], [0 imag(P1)], '--', 'Color', colAux, 'LineWidth', lwAux);
plot([0 real(P2)], [0 imag(P2)], '--', 'Color', colAux, 'LineWidth', lwAux);
plot([0 real(P3)], [0 imag(P3)], '--', 'Color', colAux, 'LineWidth', lwAux);

%% --- Head-to-tail component phasors ---
draw_arrow_complex(O,  Eas, colEas, lwComp);
draw_arrow_complex(P1, Es,  colEs,  lwComp);
draw_arrow_complex(P2, Esg, colEsg, lwComp);
draw_arrow_complex(P3, Egv, colEgv, lwComp);

%% --- Resultant phasor from origin ---
draw_arrow_complex(O, Ebar, colRes, lwRes);

%% --- Mark points lightly ---
plot(real([O P1 P2 P3 P4]), imag([O P1 P2 P3 P4]), 'o', ...
    'MarkerSize',5, 'MarkerFaceColor',[0.2 0.2 0.2], ...
    'MarkerEdgeColor','none');

%% --- Labels for component phasors ---
label_phasor(O,  Eas, '$E_{as}$', [-0.10, 0.08], 28);
label_phasor(P1, Es,  '$E_s$',    [-0.03, 0.08], 28);
label_phasor(P2, Esg, '$E_{sg}$', [ 0.02,-0.02], 28);
label_phasor(P3, Egv, '$E_{gv}$', [ 0.02,-0.03], 28);

% Resultant label
midRes = 0.58 * Ebar;
text(real(midRes)+0.10, imag(midRes)+0.02, '$\bar{E}$', ...
    'Interpreter','latex', 'FontSize',30, 'FontWeight','bold', 'Color',colRes);

text(real(midRes)+0.12, imag(midRes)-0.20, 'Measured phase', ...
    'FontName','Times New Roman', 'FontSize',22, 'FontWeight','bold', 'Color',colRes);

%% --- Measured phase angle arc ---
phi = angle(Ebar);
rArc = 0.32;
th = linspace(0, phi, 200);
plot(rArc*cos(th), rArc*sin(th), '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 2.2);

text(0.34*cos(phi/2)+0.03, 0.34*sin(phi/2)+0.02, '$\angle \bar{E}$', ...
    'Interpreter','latex', 'FontSize',22, 'Color',[0.20 0.20 0.20]);

%% --- Optional equation on plot ---
text(0.06, 1.63, '$\bar{E} = E_{as} + E_s + E_{sg} + E_{gv}$', ...
    'Interpreter','latex', 'FontSize',24, 'FontWeight','bold', 'Color',[0.1 0.1 0.1]);

%% --- Clean tick labels if you want minimalist style ---
% Uncomment next two lines for cleaner poster look:
% xticklabels({});
% yticklabels({});

end

function draw_arrow_complex(z0, dz, colorIn, lw)
%DRAW_ARROW_COMPLEX Draw arrow from z0 to z0+dz using quiver.
quiver(real(z0), imag(z0), real(dz), imag(dz), 0, ...
    'Color', colorIn, ...
    'LineWidth', lw, ...
    'MaxHeadSize', 0.25, ...
    'AutoScale', 'off');
end

function label_phasor(z0, dz, labelStr, offset, fs)
%LABEL_PHASOR Place label near the midpoint of a phasor.
zm = z0 + 0.55*dz;
text(real(zm)+offset(1), imag(zm)+offset(2), labelStr, ...
    'Interpreter','latex', ...
    'FontSize', fs, ...
    'FontWeight','bold', ...
    'Color', [0.15 0.15 0.15]);
end
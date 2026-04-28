function coh = coherence_mc_4component(theta_i, f, Hs, eps_snow, eps_g0, pars1, pars2, Nmc)
%COHERENCE_MC_4COMPONENT Monte Carlo coherence for 4-component snow/soil model.
%
% Coherent sum model:
%   E = Eas + Es + Esg + Egv
%
% This MC applies partially correlated random phases to the components that
% represent temporal variability (surface diffuse leakage, cross-pol leakage,
% and optional co-pol phase knobs), then estimates complex coherence:
%   gamma = <E1 * conj(E2)> / sqrt(<|E1|^2><|E2|^2>)
%
% Output:
%   coh.HH, coh.VV, coh.HV, coh.VH, coh.mat (2x2)

% Force matrix mode for consistent outputs
pars1.pol_mode = "matrix";
pars2.pol_mode = "matrix";

% Temporal decorrelation knobs (0..1)
[rho_snow, rho_soil] = read_temporal_rhos(pars1);

S12 = complex(zeros(2,2));
S11 = zeros(2,2);
S22 = zeros(2,2);

for n = 1:Nmc
    % Draw phases for pass 1 and independent phases
    ph = draw_component_phases();

    % Build pass-2 phases with desired temporal correlation
    ph2 = mix_component_phases(ph, rho_snow, rho_soil);

    % Apply phases into pars for each pass (component-aware)
    p1 = apply_phases_to_pars_4comp(pars1, ph,  "pass1");
    p2 = apply_phases_to_pars_4comp(pars2, ph2, "pass2");

    % Run forward model for each pass
    o1 = insar_forward_snow_soil_4component(theta_i, f, Hs, eps_snow, eps_g0, p1);
    o2 = insar_forward_snow_soil_4component(theta_i, f, Hs, eps_snow, eps_g0, p2);

    E1 = ensure_2x2(o1.E);
    E2 = ensure_2x2(o2.E);

    S12 = S12 + E1 .* conj(E2);
    S11 = S11 + abs(E1).^2;
    S22 = S22 + abs(E2).^2;
end

S12 = S12 / Nmc; S11 = S11 / Nmc; S22 = S22 / Nmc;
cohMat = S12 ./ sqrt(S11 .* S22);

coh = struct();
coh.HH = cohMat(1,1);
coh.HV = cohMat(1,2);
coh.VH = cohMat(2,1);
coh.VV = cohMat(2,2);
coh.mat = cohMat;
coh.rho_snow = rho_snow;
coh.rho_soil = rho_soil;
end

% ========================= helpers =========================

function [rho_snow, rho_soil] = read_temporal_rhos(pars)
rho_snow = 1.0; rho_soil = 1.0;

if isfield(pars,'temporal') && ~isempty(pars.temporal)
    t = pars.temporal;
    if isfield(t,'rho_snow') && ~isempty(t.rho_snow), rho_snow = t.rho_snow; end
    if isfield(t,'rho_soil') && ~isempty(t.rho_soil), rho_soil = t.rho_soil; end
    if isfield(t,'rho') && ~isempty(t.rho)
        % fallback
        rho_snow = t.rho;
        rho_soil = t.rho;
    end
end

rho_snow = max(0,min(1,rho_snow));
rho_soil = max(0,min(1,rho_soil));
end

function ph = draw_component_phases()
% Pass-1 phases (reference) + independent phases for mixing
% Surface kernels have: psi (diffuse), psi_x (xpol diffuse), psi0 (coh knob)
% Snow volume has: snow.psi0 and optionally snow.psi_x
% Soil volume: allow optional phase knob (soil_volume.psi0) if you want later

ph = struct();

% --- Snow surface ---
ph.snowSurf.psi1   = 2*pi*rand;  ph.snowSurf.psiI   = 2*pi*rand;
ph.snowSurf.psi_x1 = 2*pi*rand;  ph.snowSurf.psi_xI = 2*pi*rand;
ph.snowSurf.psi01  = 2*pi*rand;  ph.snowSurf.psi0I  = 2*pi*rand;

% --- Snow volume ---
ph.snowVol.psi01   = 2*pi*rand;  ph.snowVol.psi0I   = 2*pi*rand;
ph.snowVol.psi_x1  = 2*pi*rand;  ph.snowVol.psi_xI  = 2*pi*rand;

% --- Soil surface ---
ph.soilSurf.psi1   = 2*pi*rand;  ph.soilSurf.psiI   = 2*pi*rand;
ph.soilSurf.psi_x1 = 2*pi*rand;  ph.soilSurf.psi_xI = 2*pi*rand;
ph.soilSurf.psi01  = 2*pi*rand;  ph.soilSurf.psi0I  = 2*pi*rand;

% --- Soil volume (optional) ---
ph.soilVol.psi01   = 2*pi*rand;  ph.soilVol.psi0I   = 2*pi*rand;
end

function ph2 = mix_component_phases(ph, rho_snow, rho_soil)
% Create pass-2 phases with partial correlation
ph2 = ph;

% Snow surface
ph2.snowSurf.psi2   = phase_mix(ph.snowSurf.psi1,   ph.snowSurf.psiI,   rho_snow);
ph2.snowSurf.psi_x2 = phase_mix(ph.snowSurf.psi_x1, ph.snowSurf.psi_xI, rho_snow);
ph2.snowSurf.psi02  = phase_mix(ph.snowSurf.psi01,  ph.snowSurf.psi0I,  rho_snow);

% Snow volume
ph2.snowVol.psi02   = phase_mix(ph.snowVol.psi01,   ph.snowVol.psi0I,   rho_snow);
ph2.snowVol.psi_x2  = phase_mix(ph.snowVol.psi_x1,  ph.snowVol.psi_xI,  rho_snow);

% Soil surface
ph2.soilSurf.psi2   = phase_mix(ph.soilSurf.psi1,   ph.soilSurf.psiI,   rho_soil);
ph2.soilSurf.psi_x2 = phase_mix(ph.soilSurf.psi_x1, ph.soilSurf.psi_xI, rho_soil);
ph2.soilSurf.psi02  = phase_mix(ph.soilSurf.psi01,  ph.soilSurf.psi0I,  rho_soil);

% Soil volume (optional)
ph2.soilVol.psi02   = phase_mix(ph.soilVol.psi01,   ph.soilVol.psi0I,   rho_soil);
end

function p = apply_phases_to_pars_4comp(pars, ph, whichPass)
% Apply component phases into pars for the 4-component engine.
% Also performs backward-compatible aliasing if user only provided pars.rough.

p = pars;

% Ensure the engine sees snow_surface / soil_surface
p = ensure_surface_structs(p);

% Determine phase fields for pass1 vs pass2
if strcmp(whichPass,"pass1")
    snowSurf_psi   = ph.snowSurf.psi1;
    snowSurf_psi_x = ph.snowSurf.psi_x1;
    snowSurf_psi0  = ph.snowSurf.psi01;

    snowVol_psi0   = ph.snowVol.psi01;
    snowVol_psi_x  = ph.snowVol.psi_x1;

    soilSurf_psi   = ph.soilSurf.psi1;
    soilSurf_psi_x = ph.soilSurf.psi_x1;
    soilSurf_psi0  = ph.soilSurf.psi01;

    soilVol_psi0   = ph.soilVol.psi01;
else
    snowSurf_psi   = ph.snowSurf.psi2;
    snowSurf_psi_x = ph.snowSurf.psi_x2;
    snowSurf_psi0  = ph.snowSurf.psi02;

    snowVol_psi0   = ph.snowVol.psi02;
    snowVol_psi_x  = ph.snowVol.psi_x2;

    soilSurf_psi   = ph.soilSurf.psi2;
    soilSurf_psi_x = ph.soilSurf.psi_x2;
    soilSurf_psi0  = ph.soilSurf.psi02;

    soilVol_psi0   = ph.soilVol.psi02;
end

% --- Snow surface phases ---
if isfield(p,'snow_surface') && ~isempty(p.snow_surface)
    if isfield(p.snow_surface,'Cdiff') && p.snow_surface.Cdiff > 0
        p.snow_surface.psi = snowSurf_psi;
    end
    if isfield(p.snow_surface,'Xpol') && p.snow_surface.Xpol > 0
        p.snow_surface.psi_x = snowSurf_psi_x;
    end
    % coherent knob (optional, can be strong)
    if isfield(p.snow_surface,'psi0')
        p.snow_surface.psi0 = snowSurf_psi0;
    else
        % allow it anyway
        p.snow_surface.psi0 = snowSurf_psi0;
    end
end

% --- Snow volume phases ---
if ~isfield(p,'snow') || isempty(p.snow), p.snow = struct(); end
p.snow.psi0 = snowVol_psi0;
% only matters if you enable snow.Xpol > 0 in the engine
p.snow.psi_x = snowVol_psi_x;

% --- Soil surface phases ---
if isfield(p,'soil_surface') && ~isempty(p.soil_surface)
    if isfield(p.soil_surface,'Cdiff') && p.soil_surface.Cdiff > 0
        p.soil_surface.psi = soilSurf_psi;
    end
    if isfield(p.soil_surface,'Xpol') && p.soil_surface.Xpol > 0
        p.soil_surface.psi_x = soilSurf_psi_x;
    end
    if isfield(p.soil_surface,'psi0')
        p.soil_surface.psi0 = soilSurf_psi0;
    else
        p.soil_surface.psi0 = soilSurf_psi0;
    end
end

% --- Soil volume (optional phase knob) ---
% If later you want soil volume temporal randomness, define soil_volume.psi0 usage in the engine.
if isfield(p,'soil_volume') && ~isempty(p.soil_volume)
    p.soil_volume.psi0 = soilVol_psi0;
end

end

function p = ensure_surface_structs(p)
% Backward compatibility:
% - If user only set pars.rough, mirror it into soil_surface.
% - If snow_surface not defined, you can optionally mirror soil_surface.
if ~isfield(p,'soil_surface') || isempty(p.soil_surface)
    if isfield(p,'rough') && ~isempty(p.rough)
        p.soil_surface = p.rough;
        p.soil_surface.enable = true;
    end
end
if ~isfield(p,'snow_surface') || isempty(p.snow_surface)
    % Default: mirror soil_surface params but leave disabled unless explicitly enabled
    if isfield(p,'soil_surface') && ~isempty(p.soil_surface)
        p.snow_surface = p.soil_surface;
        if ~isfield(p.snow_surface,'enable'), p.snow_surface.enable = false; end
    end
end
end

function X2 = ensure_2x2(X)
if isempty(X), X2 = complex(zeros(2,2));
elseif isscalar(X), X2 = [X 0; 0 X];
else, X2 = X;
end
end

function psi2 = phase_mix(psi1, psi_ind, rho)
a = rho;
b = sqrt(max(1 - rho^2, 0));
z = a*exp(1i*psi1) + b*exp(1i*psi_ind);
psi2 = angle(z);
end

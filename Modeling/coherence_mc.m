function coh = coherence_mc(theta_i,f,Hs,eps_snow,eps_g0,pars1,pars2,Nmc)
% Returns coherence for HH,VV,HV,VH from Monte Carlo ensemble.

% Force matrix mode so we always get 2x2
pars1.pol_mode = "matrix";
pars2.pol_mode = "matrix";

S12 = complex(zeros(2,2));
S11 = zeros(2,2);
S22 = zeros(2,2);

% --- Temporal decorrelation knobs (0..1) ---
rho_soil = 1.0;
rho_snow = 1.0;

if isfield(pars1,'temporal') && ~isempty(pars1.temporal)
    if isfield(pars1.temporal,'rho_soil') && ~isempty(pars1.temporal.rho_soil)
        rho_soil = pars1.temporal.rho_soil;
    elseif isfield(pars1.temporal,'rho') && ~isempty(pars1.temporal.rho)
        rho_soil = pars1.temporal.rho;
    end
    if isfield(pars1.temporal,'rho_snow') && ~isempty(pars1.temporal.rho_snow)
        rho_snow = pars1.temporal.rho_snow;
    elseif isfield(pars1.temporal,'rho') && ~isempty(pars1.temporal.rho)
        rho_snow = pars1.temporal.rho;
    end
end

rho_soil = max(0,min(1,rho_soil));
rho_snow = max(0,min(1,rho_snow));

% for n = 1:Nmc
%     o1 = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars1);
%     o2 = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars2);
% 
%     E1 = ensure_2x2(o1.E);
%     E2 = ensure_2x2(o2.E);
% 
%     S12 = S12 + E1 .* conj(E2);
%     S11 = S11 + abs(E1).^2;
%     S22 = S22 + abs(E2).^2;
% end

for n = 1:Nmc

% ---- Pass 1 phases (reference) ----
psi1_soil   = 2*pi*rand;      % soil diffuse
psi1_xsoil  = 2*pi*rand;      % soil cross-pol diffuse
psi1_xsnow  = 2*pi*rand;      % snow cross-pol
psi1_snow0  = 2*pi*rand;      % NEW: snow co-pol (HH/VV) phase

% ---- Independent phases (for decorrelation) ----
psiI_soil   = 2*pi*rand;
psiI_xsoil  = 2*pi*rand;
psiI_xsnow  = 2*pi*rand;
psiI_snow0  = 2*pi*rand;      % NEW: independent co-pol snow phase

% ---- Pass 2 phases (partially correlated) ----
psi2_soil   = phase_mix(psi1_soil,  psiI_soil,  rho_soil);
psi2_xsoil  = phase_mix(psi1_xsoil, psiI_xsoil, rho_soil);
psi2_xsnow  = phase_mix(psi1_xsnow, psiI_xsnow, rho_snow);
psi2_snow0  = phase_mix(psi1_snow0, psiI_snow0, rho_snow);   % NEW

psi1_soil0 = 2*pi*rand;
psiI_soil0 = 2*pi*rand;
psi2_soil0 = phase_mix(psi1_soil0, psiI_soil0, rho_soil);

pars1n = pars1;
pars2n = pars2;

pars1n.rough.psi0 = psi1_soil0;
pars2n.rough.psi0 = psi2_soil0;

% --- SOIL diffuse phase ---
if isfield(pars1n,'rough')
    if isfield(pars1n.rough,'Cdiff') && pars1n.rough.Cdiff > 0
        pars1n.rough.psi = psi1_soil;
        pars2n.rough.psi = psi2_soil;
    end
    if isfield(pars1n.rough,'Xpol') && pars1n.rough.Xpol > 0
        pars1n.rough.psi_x = psi1_xsoil;
        pars2n.rough.psi_x = psi2_xsoil;
    end
end

% --- SNOW phases (cross-pol + co-pol) ---
if ~isfield(pars1n,'snow') || isempty(pars1n.snow)
    pars1n.snow = struct();
end
if ~isfield(pars2n,'snow') || isempty(pars2n.snow)
    pars2n.snow = struct();
end

% NEW: co-pol snow phase used by Es0
pars1n.snow.psi0 = psi1_snow0;
pars2n.snow.psi0 = psi2_snow0;

% Existing: cross-pol snow phase (only used if Xpol>0)
if isfield(pars1n.snow,'Xpol') && pars1n.snow.Xpol > 0
    pars1n.snow.psi_x = psi1_xsnow;
    pars2n.snow.psi_x = psi2_xsnow;
end

    o1 = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars1n);
    o2 = forward_realization(theta_i,f,Hs,eps_snow,eps_g0,pars2n);

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
end

function X2 = ensure_2x2(X)
if isempty(X), X2 = complex(zeros(2,2));
elseif isscalar(X), X2 = [X 0; 0 X];
else, X2 = X;
end
end

function psi2 = phase_mix(psi1, psi_ind, rho)
%PHASE_MIX Create partially correlated phase for pass 2.
% rho=1 -> psi2=psi1, rho=0 -> psi2~psi_ind
a = rho;
b = sqrt(max(1 - rho^2, 0));
z = a*exp(1i*psi1) + b*exp(1i*psi_ind);
psi2 = angle(z);
end


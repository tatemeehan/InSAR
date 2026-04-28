function out = insar_forward_snow_soil_4component(theta_i_deg, f_Hz, Hs, eps_s, eps_g, pars)
%INSAR_FORWARD_SNOW_SOIL_4COMPONENT
% Coherent 4-component forward model:
%   E = Eas + Es + Esg + Egv
%
% Eas: air/snow rough surface return (top interface)
% Es : snow volume return (depth-integrated benchmark)
% Esg: snow/soil rough surface return (bottom interface)
% Egv: soil volume return (distributed within soil)
%
% Notes:
% - Specular-pixel coherent model: single conserved kx, medium-specific kz.
% - Roughness kernels are "effective returns into the pixel" (coh + diffuse leakage proxy).

% --- Snow kz + geometry ---
[kzs, beta_s, alpha_s, k0, kx, lambda] = kz_from_eps(eps_s, theta_i_deg, f_Hz);

% Air kz (eps=1)
kza = sqrt(k0^2 - kx^2);
if imag(kza) < 0, kza = -kza; end

% --- Soil kz (uniform default; overridden by profile mode) ---
use_profile = isfield(pars,'soil_profile') && isfield(pars.soil_profile,'enable') && pars.soil_profile.enable;

Ig = [];
profDiag = [];
eps_g_surf = eps_g;

if use_profile
    % Ig models soil volume integral; profDiag provides surface kz, eps
    [Ig, profDiag] = soil_profile_integral(f_Hz, kx, pars.soil_profile, pars.Lg);
    kzg = profDiag.kz(1);
    eps_g_surf = profDiag.eps(1);
else
    [kzg, ~, ~] = kz_from_eps(eps_g, theta_i_deg, f_Hz);
end

% --- Transmission coefficients (air<->snow) ---
[T_as_HH, T_sa_HH] = fresnel_T(kza, kzs, 1, eps_s, "HH");  % TE
[T_as_VV, T_sa_VV] = fresnel_T(kza, kzs, 1, eps_s, "VV");  % TM

TasTsa_mat = [T_as_HH*T_sa_HH, T_as_HH*T_sa_VV;
              T_as_VV*T_sa_HH, T_as_VV*T_sa_VV];

% Resolve polarization mode
pol_mode = "single";
if isfield(pars,'pol_mode') && ~isempty(pars.pol_mode)
    pol_mode = string(pars.pol_mode);
end
pol_mode = lower(pol_mode);

pol = "HH";
if isfield(pars,'pol') && ~isempty(pars.pol), pol = string(pars.pol); end
pol = upper(pol);
if pol=="TE", pol="HH"; end
if pol=="TM", pol="VV"; end

if pol_mode=="matrix"
    TasTsa = TasTsa_mat;
else
    switch pol
        case "HH", TasTsa = TasTsa_mat(1,1);
        case "VV", TasTsa = TasTsa_mat(2,2);
        case "HV", TasTsa = TasTsa_mat(1,2);
        case "VH", TasTsa = TasTsa_mat(2,1);
        otherwise, error("Unknown pars.pol = %s", pol);
    end
end

% --- Two-way snow propagation to bottom interface ---
% Ts = exp(-(2*alpha_s + 1i*2*beta_s)*Hs);
% Two-way *excess* phase/attenuation replacing air with snow
% Ts = exp( -(2*alpha_s + 1i*2*(beta_s - kza))*Hs );
Ts_abs    = exp(-(2*alpha_s + 1i*2*beta_s)*Hs);          % absolute round-trip in snow
Ts_excess = exp(-(2*alpha_s + 1i*2*(beta_s - kza))*Hs);  % excess round-trip vs air
Ts_amp = abs(Ts_abs);                 % attenuation (two-way)
Ts_ph  = Ts_excess ./ abs(Ts_excess); % unit-magnitude phasor carrying excess phase
Ts_eff = Ts_amp .* Ts_ph;             % combine: realistic amplitude + correct InSAR phase

% % --- Two-way propagation through snow ---
% Ts_abs = exp(-(2*alpha_s + 1i*2*beta_s)*Hs);   % absolute round-trip in snow
% 
% % --- Data-matching excess phase convention ---
% % Positive modeled phase for increased snow-induced delay in the later pass,
% % when interferogram is formed as angle(E2 .* conj(E1))
% Ts_excess = exp(-(2*alpha_s - 1i*2*(beta_s - kza))*Hs);
% 
% % Split amplitude and phase
% Ts_amp = abs(Ts_abs);                   % physical two-way attenuation
% Ts_ph  = Ts_excess ./ abs(Ts_excess);   % unit-magnitude excess-phase phasor
% 
% % Final effective snow factor for subsurface terms
% Ts_eff = Ts_amp .* Ts_ph;

% =========================
% 1) Air/Snow rough surface
% =========================
Eas = 0;
snowSurfDbg = struct();
if isfield(pars,'snow_surface') && isfield(pars.snow_surface,'enable') && pars.snow_surface.enable
    % Use roughness kernel to get effective surface return into pixel
    if pol_mode=="matrix"
        % build 2x2 Eas matrix
        Eas_HH = snow_roughness_kernel(kza,kzs,kx,pars.snow_surface,1,eps_s,"HH");
        Eas_VV = snow_roughness_kernel(kza,kzs,kx,pars.snow_surface,1,eps_s,"VV");
        Eas_HV = snow_roughness_kernel(kza,kzs,kx,pars.snow_surface,1,eps_s,"HV");
        Eas_VH = snow_roughness_kernel(kza,kzs,kx,pars.snow_surface,1,eps_s,"VH");
        Eas = [Eas_HH, Eas_HV; Eas_VH, Eas_VV];
    else
        [Eas, snowSurfDbg] = snow_roughness_kernel(kza,kzs,kx,pars.snow_surface,1,eps_s,pol);
    end
else
    % Optional: smooth Fresnel top reflection if you want it on by default
    if isfield(pars,'snow_surface') && isfield(pars.snow_surface,'A0') && pars.snow_surface.A0 ~= 0
        % If user provided scale but not enabled, treat as smooth specular
        Ras = fresnel_R(kza,kzs,1,eps_s,pol);
        Eas = pars.snow_surface.A0 * Ras;
    end
end

% =================
% 2) Snow volume Es
% =================
Ls = pars.Ls; invLs = 0; if isfinite(Ls) && Ls > 0, invLs = 1/Ls; end
As = pars.As;

Ds = (invLs + 2*alpha_s) + 1i*(2*beta_s);
if Hs <= 0 || As == 0
    Es0 = 0;
else
    Es0 = As * (1 - exp(-Ds*Hs)) / Ds;
end

% optional snow phase knob (you already use psi0)
psi_s0 = 0.0;
if isfield(pars,'snow') && isfield(pars.snow,'psi0') && ~isempty(pars.snow.psi0)
    psi_s0 = pars.snow.psi0;
end
Es0 = Es0 * exp(1i*psi_s0);

% Build Es in scalar or matrix mode, using your existing depol knobs if desired
Es = Es0;
if pol_mode=="matrix"
    Xsnow = 0; psi_sx = 0; VV_scale = 1;
    if isfield(pars,'snow')
        if isfield(pars.snow,'Xpol'), Xsnow = pars.snow.Xpol; end
        if isfield(pars.snow,'psi_x'), psi_sx = pars.snow.psi_x; end
        if isfield(pars.snow,'VV_scale'), VV_scale = pars.snow.VV_scale; end
    end
    Es = [Es0, Xsnow*Es0*exp(1i*psi_sx);
          Xsnow*Es0*exp(1i*psi_sx), VV_scale*Es0];
else
    % single: map channel similarly
    if pol=="HV" || pol=="VH"
        Xsnow = 0; psi_sx = 0;
        if isfield(pars,'snow')
            if isfield(pars.snow,'Xpol'), Xsnow = pars.snow.Xpol; end
            if isfield(pars.snow,'psi_x'), psi_sx = pars.snow.psi_x; end
        end
        Es = Xsnow * Es0 * exp(1i*psi_sx);
    elseif pol=="VV"
        VV_scale = 1;
        if isfield(pars,'snow') && isfield(pars.snow,'VV_scale'), VV_scale = pars.snow.VV_scale; end
        Es = VV_scale * Es0;
    else
        Es = Es0;
    end
end

% Apply air-snow coupling (your request): snow volume should include TasTsa
if isfield(pars,'apply_TasTsa_to_snowvol') && ~isempty(pars.apply_TasTsa_to_snowvol)
    if pars.apply_TasTsa_to_snowvol
        Es = TasTsa * Es;
    end
else
    % default ON (per your request)
    Es = TasTsa * Es;
end

% =========================
% 3) Snow/Soil rough surface
% =========================
Esg = 0;
soilSurfDbg = struct();
if isfield(pars,'soil_surface') && isfield(pars.soil_surface,'enable') && pars.soil_surface.enable
    if pol_mode=="matrix"
        Asg_HH = soil_roughness_kernel(kzs,kzg,kx,pars.soil_surface,eps_s,eps_g_surf,"HH");
        Asg_VV = soil_roughness_kernel(kzs,kzg,kx,pars.soil_surface,eps_s,eps_g_surf,"VV");
        Asg_HV = soil_roughness_kernel(kzs,kzg,kx,pars.soil_surface,eps_s,eps_g_surf,"HV");
        Asg_VH = soil_roughness_kernel(kzs,kzg,kx,pars.soil_surface,eps_s,eps_g_surf,"VH");
        Asg = [Asg_HH, Asg_HV; Asg_VH, Asg_VV];
    else
        [Asg, soilSurfDbg] = soil_roughness_kernel(kzs,kzg,kx,pars.soil_surface,eps_s,eps_g_surf,pol);
    end

    % Two-way through snow + coupling through air-snow
    % Esg = Ts * TasTsa * Asg;
    Esg = Ts_eff * TasTsa * Asg;
end

% ===============
% 4) Soil volume
% ===============
Egv = 0;
if isfield(pars,'soil_volume') && isfield(pars.soil_volume,'enable') && pars.soil_volume.enable
    Agv = pars.soil_volume.Ag;   % interpreted as soil-volume strength (not interface Fresnel)
    Lg = pars.Lg; invLg = 0; if isfinite(Lg) && Lg > 0, invLg = 1/Lg; end

    if use_profile
        % include snow-soil transmission (both ways)
        [T_sg_HH, T_gs_HH] = fresnel_T(kzs, kzg, eps_s, eps_g_surf, "HH");
        [T_sg_VV, T_gs_VV] = fresnel_T(kzs, kzg, eps_s, eps_g_surf, "VV");
        TsgTgs_mat = [T_sg_HH*T_gs_HH, T_sg_HH*T_gs_VV;
                      T_sg_VV*T_gs_HH, T_sg_VV*T_gs_VV];

        if pol_mode=="matrix"
            TsgTgs = TsgTgs_mat;
        else
            switch pol
                case "HH", TsgTgs = TsgTgs_mat(1,1);
                case "VV", TsgTgs = TsgTgs_mat(2,2);
                case "HV", TsgTgs = TsgTgs_mat(1,2);
                case "VH", TsgTgs = TsgTgs_mat(2,1);
            end
        end

        % Egv = Ts * TasTsa * TsgTgs * Agv * Ig;
        Egv = Ts_eff * TasTsa * TsgTgs * Agv * Ig;   % profile
    else
        % uniform soil volume benchmark (your Dg form)
        [~, beta_g, alpha_g] = kz_from_eps(eps_g, theta_i_deg, f_Hz);
        Dg = (invLg + 2*alpha_g) + 1i*(2*beta_g);
        % Egv = Ts * TasTsa * Agv / Dg;
        Egv = Ts_eff * TasTsa * Agv / Dg;            % uniform
    end
end

% --- Total ---
E = Eas + Es + Esg + Egv;

% --- Pack outputs ---
out = struct();
out.E = E;
out.Eas = Eas;
out.Es  = Es;
out.Esg = Esg;
out.Egv = Egv;

out.phi = angle(E);

out.meta.lambda = lambda;
out.meta.k0 = k0;
out.meta.kx = kx;
out.kz.kza = kza;
out.kz.kzs = kzs;
out.kz.kzg = kzg;
out.eps.eps_s = eps_s;
out.eps.eps_g = eps_g;
out.eps.eps_g_surf = eps_g_surf;

out.Ts_abs = Ts_abs;
out.Ts_excess = Ts_excess;
out.Ts_eff = Ts_eff;
out.TasTsa = TasTsa;

out.profile.enable = use_profile;
out.profile.Ig = Ig;
out.profile.diag = profDiag;

out.dbg.snow_surface = snowSurfDbg;
out.dbg.soil_surface = soilSurfDbg;

% --- Optional output phase-convention mapping ---
out.output_phase_convention = "native";
if isfield(pars,'output_phase_convention') && ~isempty(pars.output_phase_convention)
    out = apply_output_phase_convention(out, string(pars.output_phase_convention));
end
end

% ---- Fresnel helpers ----
function R12 = fresnel_R(kz1, kz2, eps1, eps2, pol)
pol = upper(string(pol));
if pol=="HH", pol="TE"; end
if pol=="VV", pol="TM"; end

if pol=="TE"
    R12 = (kz1 - kz2) ./ (kz1 + kz2);
elseif pol=="TM"
    R12 = (eps2.*kz1 - eps1.*kz2) ./ (eps2.*kz1 + eps1.*kz2);
else
    error("Unknown pol for fresnel_R: %s", pol);
end
end

function [T12, T21] = fresnel_T(kz1, kz2, eps1, eps2, pol)
pol = upper(string(pol));
if pol=="HH", pol="TE"; end
if pol=="VV", pol="TM"; end

if pol=="TE"
    T12 = 2*kz1 ./ (kz1 + kz2);
    T21 = 2*kz2 ./ (kz2 + kz1);
elseif pol=="TM"
    T12 = 2*eps2.*kz1 ./ (eps2.*kz1 + eps1.*kz2);
    T21 = 2*eps1.*kz2 ./ (eps2.*kz1 + eps1.*kz2);
else
    error("Unknown pol for fresnel_T: %s", pol);
end
end

function [kz, beta, alpha, k0, kx, lambda, theta_t] = kz_from_eps(eps_r, theta_i_deg, f_Hz)
c = 299792458;
lambda = c / f_Hz;
k0 = 2*pi / lambda;

theta = deg2rad(theta_i_deg);

% conserved tangential component (set in air)
kx = k0 * sin(theta);

kz = sqrt(k0^2 * eps_r - kx^2);
if imag(kz) < 0, kz = -kz; end

beta  = real(kz);
alpha = imag(kz);

n = sqrt(eps_r);
theta_t = asin(kx ./ (k0 .* n));
end

function out = apply_output_phase_convention(out, mode)
%APPLY_OUTPUT_PHASE_CONVENTION
% Maps output phasors into the requested convention without changing the
% native internal physics.
%
% mode:
%   "native"       -> unchanged
%   "data_matched" -> conjugate phase-carrying outputs so angle() flips sign

mode = lower(string(mode));

switch mode
    case "native"
        out.output_phase_convention = "native";

    case "data_matched"
        % --- Conjugate all phase-carrying complex outputs ---
        if isfield(out,'E'),          out.E = conj(out.E); end
        if isfield(out,'Eas'),        out.Eas = conj(out.Eas); end
        if isfield(out,'Es'),         out.Es = conj(out.Es); end
        if isfield(out,'Esg'),        out.Esg = conj(out.Esg); end
        if isfield(out,'Egv'),        out.Egv = conj(out.Egv); end

        if isfield(out,'Ts_abs'),     out.Ts_abs = conj(out.Ts_abs); end
        if isfield(out,'Ts_excess'),  out.Ts_excess = conj(out.Ts_excess); end
        if isfield(out,'Ts_eff'),     out.Ts_eff = conj(out.Ts_eff); end
        if isfield(out,'TasTsa'),     out.TasTsa = conj(out.TasTsa); end

        if isfield(out,'profile') && isfield(out.profile,'Ig') && ~isempty(out.profile.Ig)
            out.profile.Ig = conj(out.profile.Ig);
        end

        % Optional: conjugate complex debug phasors if present
        if isfield(out,'dbg') && isfield(out.dbg,'snow_surface')
            out.dbg.snow_surface = conjugate_complex_fields_recursive(out.dbg.snow_surface);
        end
        if isfield(out,'dbg') && isfield(out.dbg,'soil_surface')
            out.dbg.soil_surface = conjugate_complex_fields_recursive(out.dbg.soil_surface);
        end

        % Recompute total phase from transformed output
        out.phi = angle(out.E);
        out.output_phase_convention = "data_matched";

    otherwise
        error('Unknown output_phase_convention: %s', mode);
end
end

function s = conjugate_complex_fields_recursive(s)
% Conjugate all numeric complex fields inside a struct recursively.
if ~isstruct(s), return; end

fns = fieldnames(s);
for i = 1:numel(fns)
    fn = fns{i};
    val = s.(fn);

    if isstruct(val)
        s.(fn) = conjugate_complex_fields_recursive(val);
    elseif isnumeric(val) && ~isreal(val)
        s.(fn) = conj(val);
    else
        s.(fn) = val;
    end
end
end
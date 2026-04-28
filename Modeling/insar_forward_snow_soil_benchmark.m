function out = insar_forward_snow_soil_benchmark(theta_i_deg, f_Hz, Hs, eps_s, eps_g, pars)
%INSAR_FORWARD_SNOW_SOIL_BENCHMARK
% Coherent benchmark model:
%   E = Es + Eg
%
% Snow volume scattering:
%   Es = As * (1 - exp(-Ds*Hs)) / Ds
%   Ds = (1/Ls + 2*alpha_s) + i(2*beta_s)
%
% Soil return:
%   Uniform: Eg = Ts * Ag / Dg
%            Dg = (1/Lg + 2*alpha_g) + i(2*beta_g)
%
%   Profile (optional): Eg = Ts * Ag * Ig
%     Ig is computed numerically from eps_g(z) via soil_profile_integral.m
%     Interface Ag uses surface-layer kzg(z=0) (physically natural).

% --- Snow kz + geometry ---
[kzs, beta_s, alpha_s, k0, kx, lambda] = kz_from_eps(eps_s, theta_i_deg, f_Hz);
% Air vertical wavenumber (eps=1)
kza = sqrt(k0^2 - kx^2);
if imag(kza) < 0, kza = -kza; end

% --- Two-way air-snow transmission matrix (2x2, compatible with HH/VV/HV/VH) ---
% Fresnel transmission does NOT create cross-pol; this just applies correct
% in/out factors to each matrix element.

% Compute TE/TM transmission coefficients once
[T_as_HH, T_sa_HH] = fresnel_T(kza, kzs, 1, eps_s, "HH");  % TE
[T_as_VV, T_sa_VV] = fresnel_T(kza, kzs, 1, eps_s, "VV");  % TM

% Build 2x2 mapping for [HH HV; VH VV]
TasTsa_mat = [T_as_HH*T_sa_HH, T_as_HH*T_sa_VV;
              T_as_VV*T_sa_HH, T_as_VV*T_sa_VV];

% In single-pol mode, pick the requested channel’s scalar factor
TasTsa = 1;
if isfield(pars,'pol_mode') && string(pars.pol_mode)=="matrix"
    TasTsa = TasTsa_mat;
else
    pol = "HH";
    if isfield(pars,'pol') && ~isempty(pars.pol), pol = string(pars.pol); end
    pol = upper(pol);
    if pol=="TE", pol="HH"; end
    if pol=="TM", pol="VV"; end

    switch pol
        case {"HH"}
            TasTsa = TasTsa_mat(1,1);
        case {"VV"}
            TasTsa = TasTsa_mat(2,2);
        case {"HV"}
            TasTsa = TasTsa_mat(1,2);
        case {"VH"}
            TasTsa = TasTsa_mat(2,1);
        otherwise
            error("Unknown pars.pol = %s", pol);
    end
end



% --- Soil kz (uniform default; overridden by profile if enabled) ---
kzg = []; beta_g = []; alpha_g = [];


% --- Length scales ---
Ls = pars.Ls; invLs = 0; if isfinite(Ls) && Ls > 0, invLs = 1/Ls; end
Lg = pars.Lg; invLg = 0; if isfinite(Lg) && Lg > 0, invLg = 1/Lg; end

As = pars.As;        % snow scattering amplitude
Ag = pars.Ag;        % soil amplitude (fallback if roughness off)

% --- Optional soil profile mode ---
use_profile = isfield(pars,'soil_profile') && isfield(pars.soil_profile,'enable') && pars.soil_profile.enable;
% Soil scattering interpretation in profile mode:
%   "interface" : reflect at snow-soil boundary (use R12-based roughness kernel, ignore Ig)
%   "volume"    : distributed scattering within soil (use Ig, DO NOT use R12 inside Ag)
% soil_scatter_mode = "volume";
soil_scatter_mode = "interface";

if isfield(pars,'soil_profile') && isfield(pars.soil_profile,'scatter_mode') && ~isempty(pars.soil_profile.scatter_mode)
    soil_scatter_mode = string(pars.soil_profile.scatter_mode);
end
soil_scatter_mode = lower(soil_scatter_mode);

Ig = [];
profDiag = [];
eps_g_surf = eps_g;   % default (uniform)

if use_profile
    [Ig, profDiag] = soil_profile_integral(f_Hz, kx, pars.soil_profile, Lg);

    % Use surface-layer properties for interface
    kzg = profDiag.kz(1);
    beta_g  = real(kzg);
    alpha_g = imag(kzg);
    eps_g_surf = profDiag.eps(1);     % IMPORTANT for TM
else
    % Uniform soil
    [kzg, beta_g, alpha_g] = kz_from_eps(eps_g, theta_i_deg, f_Hz);
end

% roughDbg = struct();   % default
% % --- Roughness kernel (replaces Ag) ---
% if isfield(pars,'use_roughness') && pars.use_roughness
%     % Ag = soil_roughness_kernel(kzs, kzg, kx, pars.rough);
%     pol = "TE";
%     if isfield(pars,'pol') && ~isempty(pars.pol), pol = string(pars.pol); end
%     [Ag, roughDbg] = soil_roughness_kernel(kzs, kzg, kx, pars.rough, eps_s, eps_g_surf, pol);
% end

% --- Roughness kernel (replaces Ag) ---
roughDbg = struct(); % always defined

if isfield(pars,'use_roughness') && pars.use_roughness
    if pars.pol_mode=="single"
        pol = "HH";
        if isfield(pars,'pol') && ~isempty(pars.pol), pol = string(pars.pol); end
        [Ag, roughDbg] = soil_roughness_kernel(kzs, kzg, kx, pars.rough, eps_s, eps_g_surf, pol);
    else
        % matrix mode: compute Ag for HH,VV,HV,VH
        [Ag_HH, dbgHH] = soil_roughness_kernel(kzs, kzg, kx, pars.rough, eps_s, eps_g_surf, "HH");
        [Ag_VV, dbgVV] = soil_roughness_kernel(kzs, kzg, kx, pars.rough, eps_s, eps_g_surf, "VV");
        [Ag_HV, dbgHV] = soil_roughness_kernel(kzs, kzg, kx, pars.rough, eps_s, eps_g_surf, "HV");
        [Ag_VH, dbgVH] = soil_roughness_kernel(kzs, kzg, kx, pars.rough, eps_s, eps_g_surf, "VH");

        roughDbg = struct('HH',dbgHH,'VV',dbgVV,'HV',dbgHV,'VH',dbgVH);

        % Pack into 2x2 Ag matrix
        Ag = [Ag_HH, Ag_HV; Ag_VH, Ag_VV];
    end
end


% % --- Snow volume scattering ---
% Ds = (invLs + 2*alpha_s) + 1i*(2*beta_s);
% if Hs <= 0 || As == 0
%     Es = 0;
% else
%     Es = As * (1 - exp(-Ds*Hs)) / Ds;
% end

% --- Snow volume scattering (now polarization-aware) ---
Ds = (invLs + 2*alpha_s) + 1i*(2*beta_s);

if Hs <= 0 || As == 0
    Es0 = 0;
else
    Es0 = As * (1 - exp(-Ds*Hs)) / Ds;
end

% NEW: apply co-pol snow phase (temporal decorrelation)
psi_s0 = 0.0;
if isfield(pars,'snow') && ~isempty(pars.snow) && isfield(pars.snow,'psi0') && ~isempty(pars.snow.psi0)
    psi_s0 = pars.snow.psi0;
end
Es0 = Es0 * exp(1i*psi_s0);

% Snow depolarization controls (defaults)
Xsnow  = 0.0;
psi_sx = 0.0;
VV_scale = 1.0;



if isfield(pars,'snow') && ~isempty(pars.snow)
    if isfield(pars.snow,'Xpol')  && ~isempty(pars.snow.Xpol),  Xsnow  = pars.snow.Xpol;  end
    if isfield(pars.snow,'psi_x') && ~isempty(pars.snow.psi_x), psi_sx = pars.snow.psi_x; end
    if isfield(pars.snow,'VV_scale') && ~isempty(pars.snow.VV_scale), VV_scale = pars.snow.VV_scale; end
end

% Build Es to match Eg shape:
% - single mode: Es is scalar for HH/VV, or reduced for HV/VH
% - matrix mode: Es is 2x2 with optional cross-pol leakage
if isfield(pars,'pol_mode') && string(pars.pol_mode)=="matrix"
    % 2x2 snow scattering matrix
    Es_HH = Es0;
    Es_VV = VV_scale * Es0;
    Es_HV = Xsnow * Es0 * exp(1i*psi_sx);
    Es_VH = Es_HV;  % symmetric for now (can split later)
    Es = [Es_HH, Es_HV; Es_VH, Es_VV];
else
    % single-pol channel
    pol_s = "HH";
    if isfield(pars,'pol') && ~isempty(pars.pol), pol_s = string(pars.pol); end
    pol_s = upper(pol_s);

    if pol_s=="HV" || pol_s=="VH"
        Es = Xsnow * Es0 * exp(1i*psi_sx);   % cross-pol snow only if enabled
    elseif pol_s=="VV" || pol_s=="TM"
        Es = VV_scale * Es0;
    else
        Es = Es0; % HH/TE default
    end
end


% --- Two-way snow transmission to soil interface ---
Ts = exp(-(2*alpha_s + 1i*2*beta_s)*Hs);

% --- Soil return ---
% if Ag == 0
%     Eg = 0;
%     Dg = (invLg + 2*alpha_g) + 1i*(2*beta_g);
% else
%     if use_profile
%         Eg = Ts * Ag * Ig;
%         if abs(Ig) > 0
%             Dg = 1 / Ig;     % effective “denominator” for phase decomposition diagnostics
%         else
%             Dg = complex(inf, inf);
%         end
%     else
%         Dg = (invLg + 2*alpha_g) + 1i*(2*beta_g);
%         Eg = Ts * Ag / Dg;
%     end
% end

% --- Soil return ---
Dg = []; % scalar diagnostic "denominator" for the *propagation integral* (same for all pol)
if use_profile
    if abs(Ig) > 0, Dg = 1/Ig; else, Dg = complex(inf,inf); end
else
    Dg = (invLg + 2*alpha_g) + 1i*(2*beta_g);
end

% Eg computation: Ag can be scalar (single) or 2x2 (matrix mode)
if isempty(Ag) || (isscalar(Ag) && Ag==0) || (ismatrix(Ag) && all(Ag(:)==0))
    if isfield(pars,'pol_mode') && string(pars.pol_mode)=="matrix"
        Eg = complex(zeros(2,2));
    else
        Eg = 0;
    end
else
    % % if use_profile
    % %     Eg = Ts * Ag * Ig;
    % % else
    % %     Eg = Ts * Ag / Dg;
    % % end
    % if use_profile
    %     Eg = Ts * TasTsa * Ag * Ig;
    % else
    %     Eg = Ts * TasTsa * Ag / Dg;
    % end
    if use_profile
        if soil_scatter_mode=="interface"
            % Interface reflection at snow-soil boundary:
            % Use Ag as-is (R12-based roughness kernel), and ignore Ig.
            Eg = Ts * TasTsa * Ag;   % (Ts already two-way in snow)
        elseif soil_scatter_mode=="volume"
            % Volume scattering within soil:
            % MUST include snow<->soil transmission, and Ag must be a soil-volume strength
            % (not an interface reflection coefficient).
            [T_sg_HH, T_gs_HH] = fresnel_T(kzs, kzg, eps_s, eps_g_surf, "HH");
            [T_sg_VV, T_gs_VV] = fresnel_T(kzs, kzg, eps_s, eps_g_surf, "VV");
            TsgTgs_mat = [T_sg_HH*T_gs_HH, T_sg_HH*T_gs_VV;
                T_sg_VV*T_gs_HH, T_sg_VV*T_gs_VV];

            % Apply correct scalar/matrix form
            if isfield(pars,'pol_mode') && string(pars.pol_mode)=="matrix"
                TsgTgs = TsgTgs_mat;
            else
                pol = "HH";
                if isfield(pars,'pol') && ~isempty(pars.pol), pol = string(pars.pol); end
                pol = upper(pol);
                if pol=="TE", pol="HH"; end
                if pol=="TM", pol="VV"; end
                switch pol
                    case "HH", TsgTgs = TsgTgs_mat(1,1);
                    case "VV", TsgTgs = TsgTgs_mat(2,2);
                    case "HV", TsgTgs = TsgTgs_mat(1,2);
                    case "VH", TsgTgs = TsgTgs_mat(2,1);
                end
            end

            Eg = Ts * TasTsa * TsgTgs * Ag * Ig;
        else
            error("Unknown soil_scatter_mode: %s", soil_scatter_mode);
        end
    else
        Eg = Ts * TasTsa * Ag / Dg;
    end
end

E = Es + Eg;

% --- Pack outputs ---
out = struct();
out.E = E; out.Es = Es; out.Eg = Eg;
out.phi = angle(E);

out.kzs = kzs; out.kzg = kzg;
out.beta_s = beta_s; out.alpha_s = alpha_s;
out.beta_g = beta_g; out.alpha_g = alpha_g;
out.eps_g_surf = eps_g_surf;
out.kzg_surf = kzg;

out.meta.lambda = lambda;
out.meta.k0 = k0;
out.meta.kx = kx;

out.Ag = Ag;
out.Ts = Ts;
out.Dg = Dg;

out.profile = struct();
out.profile.enable = use_profile;
out.profile.Ig = Ig;
out.profile.diag = profDiag;
out.rough = roughDbg;
end

function [T12, T21] = fresnel_T(kz1, kz2, eps1, eps2, pol)
% Returns field transmission coefficients:
% T12: medium 1 -> medium 2
% T21: medium 2 -> medium 1
pol = upper(string(pol));
if pol=="HH", pol="TE"; end
if pol=="VV", pol="TM"; end

if pol=="TE"
    T12 = 2*kz1 ./ (kz1 + kz2);
    T21 = 2*kz2 ./ (kz2 + kz1);
elseif pol=="TM"
    % r_p = (eps2*kz1 - eps1*kz2)/(eps2*kz1 + eps1*kz2)
    % t_p (1->2) = 2*eps2*kz1/(eps2*kz1 + eps1*kz2)
    T12 = 2*eps2.*kz1 ./ (eps2.*kz1 + eps1.*kz2);
    T21 = 2*eps1.*kz2 ./ (eps2.*kz1 + eps1.*kz2);
else
    error("Unknown pol for fresnel_T: %s", pol);
end
end


% function out = insar_forward_snow_soil_benchmark(theta_i_deg, f_Hz, Hs, eps_s, eps_g, pars)
% 
% [kzs, beta_s, alpha_s, k0, kx, lambda] = kz_from_eps(eps_s, theta_i_deg, f_Hz);
% [kzg, beta_g, alpha_g]                = kz_from_eps(eps_g, theta_i_deg, f_Hz);
% 
% Ls = pars.Ls;  invLs = 0; if isfinite(Ls) && Ls > 0, invLs = 1/Ls; end
% Lg = pars.Lg;  invLg = 0; if isfinite(Lg) && Lg > 0, invLg = 1/Lg; end
% 
% As = pars.As;   % real (Level 1) or complex (Level 2)
% Ag = pars.Ag;
% 
% % Replace Ag with roughness-physics Ag_eff if enabled
% if isfield(pars,'use_roughness') && pars.use_roughness
%     Ag = soil_roughness_kernel(kzs, kzg, kx, pars.rough);
% end
% 
% % Snow volume scattering: z in [0, Hs]
% Ds = (invLs + 2*alpha_s) + 1i*(2*beta_s);
% if Hs <= 0 || As == 0
%     Es = 0;
% else
%     Es = As * (1 - exp(-Ds*Hs)) / Ds;
% end
% 
% % Two-way snow transmission to soil interface
% Ts = exp(-(2*alpha_s + 1i*2*beta_s)*Hs);
% 
% % Soil scattering: z' in [0, inf)
% Dg = (invLg + 2*alpha_g) + 1i*(2*beta_g);
% if Ag == 0
%     Eg = 0;
% else
%     Eg = Ts * Ag / Dg;
% end
% 
% E = Es + Eg;
% 
% out = struct();
% out.E = E; out.Es = Es; out.Eg = Eg;
% out.phi = angle(E);
% 
% out.kzs = kzs; out.kzg = kzg;
% out.beta_s = beta_s; out.alpha_s = alpha_s;
% out.beta_g = beta_g; out.alpha_g = alpha_g;
% out.meta.lambda = lambda;
% out.meta.k0 = k0;
% out.meta.kx = kx;
% out.Ag = Ag;
% out.Ts = Ts;
% out.Dg = Dg;
% end

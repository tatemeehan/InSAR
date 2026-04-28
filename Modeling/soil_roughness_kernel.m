function [Ag_eff, dbg] = soil_roughness_kernel(kzs, kzg, kx, pars, eps_s, eps_g, pol)
%SOIL_ROUGHNESS_KERNEL
% Effective complex soil scattering weight for a benchmark PolSAR kernel.
%
% Co-pol:
%   Ag = A0 * ( R12(pol)*Gamma_coh + Cdiff*Gamma_diff*exp(i*psi) )
%
% Cross-pol (benchmark depolarization):
%   Ag_xpol = A0 * ( Xpol*Gamma_diff*exp(i*psi_x) )
%
% Notes:
% - R12 is Fresnel (no cross-pol for isotropic interface)
% - Gamma_coh comes from Gaussian height stats:
%       <exp(i*dkz*h)> = exp(-0.5*sigma_h^2*dkz^2)
% - Gamma_diff is an energy-aware amplitude proxy:
%       sqrt(1-|Gamma_coh|^2) * sqrt(PSDshape(q))
% - PSDshape assumes Gaussian correlation function => Gaussian PSD
%   (still a benchmark; we can swap PSD models later)

% --- required pars ---
sigma_h = pars.sigma_h;   % RMS height [m]
ell     = pars.ell;       % corr length [m]
A0      = pars.A0;        % overall soil scattering scale
Cdiff   = pars.Cdiff;     % diffuse leakage strength (0..1)

% --- optional pars ---
psi   = 0;
if isfield(pars,'psi') && ~isempty(pars.psi), psi = pars.psi; end

Xpol  = 0;
if isfield(pars,'Xpol') && ~isempty(pars.Xpol), Xpol = pars.Xpol; end

psi_x = 0;
if isfield(pars,'psi_x') && ~isempty(pars.psi_x), psi_x = pars.psi_x; end

% --- normalize pol strings ---
if nargin < 7 || isempty(pol), pol = "HH"; end
pol = upper(string(pol));
if pol=="TE", pol="HH"; end
if pol=="TM", pol="VV"; end

% --- Fresnel reflection (co-pol only) ---
switch pol
    case "HH"  % TE
        R12 = (kzs - kzg) ./ (kzs + kzg);
    case "VV"  % TM
        R12 = (eps_g.*kzs - eps_s.*kzg) ./ (eps_g.*kzs + eps_s.*kzg);
    case {"HV","VH"}
        % No Fresnel x-pol for isotropic interface:
        R12 = 0;
    otherwise
        error('Unknown pol: %s (use HH,VV,HV,VH or TE,TM)', pol);
end

% --- Coherent roughness factor ---
% dkz = (kzs - kzg);
% Gamma_coh = exp(-0.5 * (sigma_h * dkz).^2);

qz = 2*real(kzs);                           % two-way vertical wavenumber in snow
Gamma_coh = exp(-0.5 * (sigma_h*qz).^2);    % real, positive coherent attenuation

% qz = real(kzs) + real(kzg);                 % “round trip” vertical component idea
% Gamma_coh = exp(-0.5 * (sigma_h*qz).^2);


R_coh = R12 .* Gamma_coh;

psi0 = 0;
if isfield(pars,'psi0') && ~isempty(pars.psi0)
    psi0 = pars.psi0;
end
R_coh = R_coh .* exp(1i*psi0);


% --- Diffuse amplitude proxy from PSDshape(q) ---
% For monostatic backscatter, a simple benchmark is q ≈ |2*kx|
q = abs(2*kx);
PSDshape = exp(-(q*ell).^2 / 4);             % Gaussian corr -> Gaussian PSD shape

Gamma_diff = sqrt(max(1 - abs(Gamma_coh).^2, 0)) .* sqrt(PSDshape);

% --- Cross-pol benchmark ---
if pol=="HV" || pol=="VH"
    Ag_eff = A0 .* (Xpol .* Gamma_diff .* exp(1i*psi_x));
else
    R_diff = Cdiff .* Gamma_diff .* exp(1i*psi);
    Ag_eff = A0 .* (R_coh + R_diff);
end

% --- Debug outputs ---
if nargout > 1
    dbg = struct();
    dbg.pol = pol;
    dbg.R12 = R12;
    dbg.qz = qz;
    dbg.Gamma_coh  = Gamma_coh;
    dbg.Gamma_diff = Gamma_diff;
    dbg.PSDshape   = PSDshape;
    dbg.R_coh      = R_coh;
    if ~(pol=="HV" || pol=="VH")
        dbg.R_diff = Cdiff .* Gamma_diff .* exp(1i*psi);
    else
        dbg.R_xpol = Xpol .* Gamma_diff .* exp(1i*psi_x);
    end
end
end

% function [Ag_eff,dbg] = soil_roughness_kernel(kzs, kzg, kx, pars, eps_s, eps_g, pol)
% %SOIL_ROUGHNESS_KERNEL
% % Effective complex soil scattering weight:
% %   Ag = A0 * ( R12*Gamma_coh + Cdiff*Gamma_diff*exp(i*psi) )
% %
% % Physics:
% %   - R12 is the smooth Fresnel-like reflection coefficient using kz's
% %   - Gamma_coh is the coherent attenuation from Gaussian height statistics:
% %       Gamma_coh = <exp(i Δkz h)> = exp(-0.5*sigma_h^2*Δkz^2),  Δkz = kzs-kzg
% %   - Gamma_diff is an explicit diffuse amplitude (benchmark, but tied to physics):
% %       Gamma_diff ~ sqrt(1-|Gamma_coh|^2) * sqrt(PSDshape(q))
% %       q is the scattering vector magnitude; for monostatic backscatter use q=|2*kx|
% %       PSDshape for Gaussian correlation length ell: exp(-(q*ell)^2/4)
% %
% % Tunables:
% %   - Cdiff (0..1): how much diffuse leaks back into the coherent pixel
% %   - psi: random phase for diffuse term (fixed for a realization)
% 
% % Required parameters
% sigma_h = pars.sigma_h;   % RMS height [m]
% ell     = pars.ell;       % corr length [m]
% A0      = pars.A0;        % overall soil scattering scale
% Cdiff   = pars.Cdiff;     % diffuse leakage strength (0..1)
% 
% % % Smooth reflection using vertical wavenumbers
% % R12 = (kzs - kzg) / (kzs + kzg);
% 
% % ---- Polarization handling ----
% if nargin < 7 || isempty(pol), pol = "TE"; end
% pol = upper(string(pol));
% 
% % Allow canonical SAR labels
% % HH -> TE (s-pol), VV -> TM (p-pol)
% % HV/VH -> cross-pol (benchmark diffuse-only here)
% isXpol = (pol=="HV") || (pol=="VH");
% 
% pol_eff = pol;
% if pol_eff=="HH", pol_eff="TE"; end
% if pol_eff=="VV", pol_eff="TM"; end
% 
% % ---- Smooth Fresnel reflection (complex eps, complex kz) ----
% % Only needed for co-pol (TE/TM). Cross-pol returns later.
% if ~isXpol
%     if pol_eff=="TE"
%         R12 = (kzs - kzg) ./ (kzs + kzg);
%     elseif pol_eff=="TM"
%         R12 = (eps_g.*kzs - eps_s.*kzg) ./ (eps_g.*kzs + eps_s.*kzg);
%     else
%         error('Unknown pol: %s (use TE/TM or HH/VV or HV/VH)', pol);
%     end
% else
%     R12 = complex(0); % not used for cross-pol (benchmark choice)
% end
% 
% 
% % fprintf('TE: |R12|=%.3f, ang=%.2f deg\n', abs(R12_TE), rad2deg(angle(R12_TE)));
% % fprintf('TM: |R12|=%.3f, ang=%.2f deg\n', abs(R12_TM), rad2deg(angle(R12_TM)));
% 
% 
% % --- Coherent roughness factor (derived for Gaussian heights) ---
% dkz = (kzs - kzg);
% Gamma_coh = exp(-0.5 * (sigma_h * dkz).^2);
% 
% R_coh = R12 .* Gamma_coh;
% 
% % --- Diffuse kernel shape via PSD(q) ---
% q = abs(2*kx);                               % monostatic reduction of |ks-ki|
% PSDshape = exp(-(q*ell).^2 / 4);             % Gaussian correlation -> Gaussian PSD shape
% 
% % Explicit diffuse amplitude (energy-aware benchmark)
% Gamma_diff = sqrt(max(1 - abs(Gamma_coh).^2, 0)) .* sqrt(PSDshape);
% 
% Xpol = 0;
% if isfield(pars,'Xpol') && ~isempty(pars.Xpol), Xpol = pars.Xpol; end
% psi_x = 0;
% if isfield(pars,'psi_x') && ~isempty(pars.psi_x), psi_x = pars.psi_x; end
% 
% if isXpol
%     R_xpol = Xpol .* Gamma_diff .* exp(1i*psi_x);
%     Ag_eff = A0 .* R_xpol;
% 
%     if nargout > 1
%         dbg = struct();
%         dbg.pol = pol;
%         dbg.isXpol = true;
%         dbg.R12 = R12;                 % = 0 (unused)
%         dbg.dkz = dkz;
%         dbg.Gamma_coh  = Gamma_coh;
%         dbg.Gamma_diff = Gamma_diff;
%         dbg.PSDshape   = PSDshape;
%         dbg.R_coh  = complex(0);
%         dbg.R_diff = complex(0);
%         dbg.R_xpol = R_xpol;
%     end
%     return
% end
% 
% psi = 0;
% if isfield(pars,'psi') && ~isempty(pars.psi)
%     psi = pars.psi;
% end
% 
% R_diff = Cdiff .* Gamma_diff .* exp(1i*psi);
% 
% Ag_eff = A0 .* (R_coh + R_diff);
% 
% if nargout > 1
%     dbg = struct();
%     dbg.R12 = R12;
%     dbg.dkz = dkz;
%     dbg.Gamma_coh  = Gamma_coh;
%     dbg.Gamma_diff = Gamma_diff;
%     dbg.PSDshape   = PSDshape;
%     dbg.R_coh  = R_coh;
%     dbg.R_diff = R_diff;
% end
% end
% 
% % function Ag_eff = soil_roughness_kernel(kzs, kzg, kx, pars)
% % %SOIL_ROUGHNESS_KERNEL Effective complex soil scattering weight from
% % % smooth Fresnel + coherent roughness attenuation + diffuse kernel.
% % 
% % % Required parameters
% % sigma_h = pars.sigma_h;   % RMS height [m]
% % ell     = pars.ell;       % corr length [m]
% % A0      = pars.A0;        % overall soil scattering scale
% % Cdiff   = pars.Cdiff;     % diffuse scale (0..1)
% % 
% % % Smooth reflection using vertical wavenumbers
% % R12 = (kzs - kzg) / (kzs + kzg);
% % 
% % % Coherent roughness factor (Gaussian heights)
% % % Gamma_coh = exp(-2 * (kzs * sigma_h).^2);   % complex kzs OK
% % dkz = kzs - kzg;
% % Gamma_coh = exp(-0.5 * (dkz * sigma_h).^2);
% % 
% % R_coh = R12 .* Gamma_coh;
% % 
% % % Diffuse kernel shaped by PSD at q ~ 2*kx (Gaussian correlation)
% % q = 2*kx;
% % PSDshape = exp(-(q*ell/2).^2);              % normalized shape in [0,1]
% % 
% % Gamma_diff = sqrt(max(1 - abs(Gamma_coh).^2, 0)) .* sqrt(PSDshape);
% % 
% % psi = 0;
% % if isfield(pars,'psi') && ~isempty(pars.psi)
% %     psi = pars.psi;  % user can set or randomize outside
% % end
% % 
% % R_diff = Cdiff .* Gamma_diff .* exp(1i*psi);
% % 
% % Ag_eff = A0 .* (R_coh + R_diff);
% % end

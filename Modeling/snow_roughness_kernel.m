function [Aas_eff, dbg] = snow_roughness_kernel(kza, kzs, kx, pars, eps_a, eps_s, pol)
%SNOW_ROUGHNESS_KERNEL
% Effective complex air-snow surface scattering weight, benchmark PolSAR kernel.
%
% Co-pol:
%   Aas_eff = A0 * ( R12(pol)*Gamma_coh + Cdiff*Gamma_diff*exp(i*psi) )
%
% Cross-pol:
%   Aas_xpol = A0 * ( Xpol*Gamma_diff*exp(i*psi_x) )
%
% Notes:
% - R12 is Fresnel reflection at air->snow
% - Gamma_coh uses Gaussian height stats, applied to coherent component
% - Gamma_diff is an energy-aware amplitude proxy shaped by PSD(q)

% --- required pars ---
sigma_h = pars.sigma_h;   % RMS height [m]
ell     = pars.ell;       % corr length [m]
A0      = pars.A0;        % overall surface scattering scale
Cdiff   = pars.Cdiff;     % diffuse leakage strength (0..1)

% --- optional pars ---
psi   = 0;
if isfield(pars,'psi') && ~isempty(pars.psi), psi = pars.psi; end

Xpol  = 0;
if isfield(pars,'Xpol') && ~isempty(pars.Xpol), Xpol = pars.Xpol; end

psi_x = 0;
if isfield(pars,'psi_x') && ~isempty(pars.psi_x), psi_x = pars.psi_x; end

psi0 = 0;
if isfield(pars,'psi0') && ~isempty(pars.psi0), psi0 = pars.psi0; end

% --- normalize pol strings ---
if nargin < 7 || isempty(pol), pol = "HH"; end
pol = upper(string(pol));
if pol=="TE", pol="HH"; end
if pol=="TM", pol="VV"; end

% --- Fresnel reflection (co-pol only) ---
switch pol
    case "HH"  % TE
        R12 = (kza - kzs) ./ (kza + kzs);
    case "VV"  % TM
        R12 = (eps_s.*kza - eps_a.*kzs) ./ (eps_s.*kza + eps_a.*kzs);
    case {"HV","VH"}
        R12 = 0;
    otherwise
        error('Unknown pol: %s (use HH,VV,HV,VH or TE,TM)', pol);
end

% --- Coherent roughness factor ---
% For top interface, use incident-side (air) two-way vertical wavenumber
qz = 2*real(kza);
Gamma_coh = exp(-0.5 * (sigma_h*qz).^2);   % real, positive

R_coh = R12 .* Gamma_coh;
R_coh = R_coh .* exp(1i*psi0);

% --- Diffuse amplitude proxy from PSDshape(q) ---
q = abs(2*kx);                              % monostatic proxy
PSDshape   = exp(-(q*ell).^2 / 4);          % Gaussian correlation -> Gaussian PSD shape
Gamma_diff = sqrt(max(1 - abs(Gamma_coh).^2, 0)) .* sqrt(PSDshape);

% --- Cross-pol benchmark ---
if pol=="HV" || pol=="VH"
    Aas_eff = A0 .* (Xpol .* Gamma_diff .* exp(1i*psi_x));
else
    R_diff  = Cdiff .* Gamma_diff .* exp(1i*psi);
    Aas_eff = A0 .* (R_coh + R_diff);
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

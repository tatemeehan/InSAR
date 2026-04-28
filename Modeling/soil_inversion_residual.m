function r = soil_inversion_residual(p, data, wparams)
% SOIL_INVERSION_RESIDUAL
% Residual vector for lsqnonlin, using soilPermittivity_advanced.
%
% p: parameter vector
% data: struct with fields
%   .VWC, .SWC, .EC_dSm, .Tc
%   .freq_vna      (Nv×1)
%   .eps_vna_meas  (Nv×1 complex)
%   .freq_hydra    (optional scalar)
%   .eps_hydra     (optional complex)
%   .w_vna_real, .w_vna_imag (optional weights)
%   .w_hydra       (scalar weight)
%
% wparams: water parameter struct (fixed parts)

    % --- Map p -> soilParams, tauBW, clayParams ---

    soilParams.eps_real_soil = p(1);
    soilParams.eps_real_clay = p(2);
    soilParams.frac_clay     = p(3);

    wparams.frac_bound       = p(4);

    tauBW.free_scale  = 1.0;        % keep fixed for now
    tauBW.bound_scale = p(5);
    tauBW.B_bound     = p(6);
    tauBW.max_tau     = data.tau_max;   % e.g. 5e-10

    clayParams.eps_inf_clay = data.clay_eps_inf;  % e.g. 5
    clayParams.tau_clay     = p(7);
    clayParams.alpha_clay   = p(8);

    % --- Unpack data ---

    VWC       = data.VWC;
    SWC       = data.SWC;
    EC_dSm    = data.EC_dSm;
    Tc        = data.Tc;
    sigma_Sm  = EC_dSm(:) * 0.1;        % dS/m -> S/m

    freq_vna  = data.freq_vna(:);
    eps_meas  = data.eps_vna_meas(:);
    Nv        = numel(freq_vna);

    % Weights
    if isfield(data,'w_vna_real'), w_vna_real = data.w_vna_real;
    else,                          w_vna_real = 1.0;     end
    if isfield(data,'w_vna_imag'), w_vna_imag = data.w_vna_imag;
    else,                          w_vna_imag = 1.0;     end

    % --- VNA residuals ---

    eps_model_vna = zeros(Nv,1);
    for k = 1:Nv
        f = freq_vna(k);
        eps_model_vna(k) = soilPermittivity_advanced( ...
            VWC, SWC, sigma_Sm, soilParams, Tc, f, ...
            wparams, tauBW, clayParams);
    end

    re_err = real(eps_model_vna) - real(eps_meas);
    im_err = imag(eps_model_vna) - imag(eps_meas);

    r_vna = [ w_vna_real * re_err ;
              w_vna_imag * im_err ];

    % --- HydraGO anchor at 50 MHz (optional) ---

    r_h = [];
    if isfield(data,'freq_hydra') && ~isempty(data.freq_hydra) ...
       && isfield(data,'eps_hydra') && ~isempty(data.eps_hydra)

        fH = data.freq_hydra;
        eps_model_h = soilPermittivity_advanced( ...
            VWC, SWC, sigma_Sm, soilParams, Tc, fH, ...
            wparams, tauBW, clayParams);

        eH = eps_model_h - data.eps_hydra;

        w_h = data.w_hydra;     % scalar
        r_h = w_h * [real(eH); imag(eH)];
    end

    % --- Total residual vector for lsqnonlin ---
    r = [r_vna; r_h];
end

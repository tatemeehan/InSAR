function J = soil_forward_cost(p, data, wparams)

    % Unpack parameters you want to tune
    soilParams.eps_real_soil = p(1);
    soilParams.eps_real_clay = p(2);
    soilParams.frac_clay     = p(3);
    tauBW.A         = p(4);
    tauBW.B         = p(5);
    tauBW.tau_scale = p(6);
    tauBW.max_tau   = data.tau_max;   % can keep fixed or include in p

    VWC    = data.VWC;
    SWC    = data.SWC;
    EC_dSm = data.EC_dSm;
    Tc     = data.Tc;

    % --- VNA frequencies ---
    freqs_vna  = data.freq_vna;
    eps_vna_meas = data.eps_vna_meas;   % complex vector
    
    % Forward model at all VNA freqs
    eps_model_vna = zeros(size(eps_vna_meas));
    for k = 1:numel(freqs_vna)
        eps_model_vna(k) = soilPermittivity_colecole_CRIM_bound( ...
            VWC, SWC, EC_dSm.*0.1, ...
            soilParams.eps_real_soil, soilParams.eps_real_clay, soilParams.frac_clay, ...
            Tc, freqs_vna(k), wparams, tauBW);
    end

    % Error at VNA freqs
    w_vna = 1.0;
    err_vna = eps_model_vna - eps_vna_meas;
    cost_vna = w_vna * mean(abs(err_vna).^2);

    % --- HydraGO anchor at 50 MHz (optional) ---
    if isfield(data,'freq_hydra') && ~isempty(data.freq_hydra)
        fH = data.freq_hydra;
        eps_model_h = soilPermittivity_colecole_CRIM_bound( ...
            VWC, SWC, EC_dSm*0.1, ...
            soilParams.eps_real_soil, soilParams.eps_real_clay, soilParams.frac_clay, ...
            Tc, fH, wparams, tauBW);

        err_h = eps_model_h - data.eps_hydra;    % complex error at 50 MHz
        w_h   = data.weight_hydra;
        cost_h = w_h * (real(err_h)^2 + imag(err_h)^2);
    else
        cost_h = 0;
    end

    % Optional small regularization to keep parameters in sane range
    reg = 0.0;
    if isfield(data,'lambda_reg')
        reg = data.lambda_reg * sum(p.^2);
    end

    J = cost_vna + cost_h + reg;
end

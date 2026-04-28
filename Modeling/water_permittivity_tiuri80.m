function Ew = water_permittivity_tiuri80(frequency, temperature)
    % Calculates the complex water dielectric constant based on Tiuri 1980.
    %
    % Parameters:
    % frequency: Frequency in Hz
    % temperature: Temperature in Kelvin
    %
    % Returns:
    % Ew: Complex dielectric constant of water

    % Constants
    GHz = 1e9; % GHz to Hz conversion factor
    FREEZING_POINT = 373.15; % Freezing point of water in Kelvin
    freqGHz = frequency / GHz; % Convert frequency to GHz
    tempC = temperature - FREEZING_POINT; % Convert temperature to Celsius

    % Validate temperature
    if any(tempC > 0)
        error('The water temperature must not be higher or equal to %f K', FREEZING_POINT);
    end

    % Parameters from Tiuri (1980)
    e2 = 4.903e-2; % Static loss factor
    e1 = 87.74 - 0.4008 * tempC + 9.398e-4 * tempC.^2 + 1.410e-6 * tempC.^3;

    % Relaxation frequency based on Liebe 1991
    theta = 1 - 300.0 ./ temperature;
    f1 = 20.2 + 146.4 .* theta + 316 .* theta.^2;

    % Calculate complex dielectric constant
    Ew = e2 + (e1 - e2) ./ complex(1, -freqGHz ./ f1);
end

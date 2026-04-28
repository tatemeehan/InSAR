
function Ew = water_permittivity_maetzler87(frequency, temperature)
% Calculates the complex water dielectric constant depending on the frequency and temperature
% Based on Mätzler, C., & Wegmuller, U. (1987). Dielectric properties of freshwater
% ice at microwave frequencies. Journal of Physics D: Applied Physics, 20(12), 1623-1630.
% 
% Inputs:
%   - frequency: Frequency in Hz
%   - temperature: Temperature in K
% 
% Outputs:
%   - Ew: Complex permittivity of pure ice
% 
% Raises:
%   - Error if temperature > 273.15 K (liquid water not suitable for this model)

    % Check if temperature is within the valid range
    if temperature > 273.15
        error('The model is not suitable for temperatures above freezing (liquid water).');
    end

    % Convert frequency to GHz
    freqGHz = frequency / 1e9;

    % Compute theta (dimensionless temperature parameter)
    theta = 1 - 300.0 / temperature;

    % Temperature-dependent parameters
    e0 = 77.66 - 103.3 * theta; % Static permittivity
    e1 = 0.0671 * e0;           % High-frequency permittivity

    % Relaxation frequencies (GHz)
    f1 = 20.2 + 146.4 * theta + 316 * theta^2;
    f2 = 39.8 * f1;

    % Asymptotic high-frequency value
    e2 = 3.52 + 7.52 * theta;

    % Compute the complex permittivity
    Ew = e2 + (e1 - e2) ./ (1 - 1i * freqGHz / f2) + (e0 - e1) ./ (1 - 1i * freqGHz / f1);
end
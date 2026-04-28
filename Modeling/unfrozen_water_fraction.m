function [fL, info] = unfrozen_water_fraction(Tc, VWC, SWC, soilParams, uwf)
%UNFROZEN_WATER_FRACTION  Liquid fraction of pore water under freezing.
%
% Returns fL in [fmin, 1], where:
%   theta_w_liq = fL .* VWC
%   theta_w_ice = (1-fL) .* VWC
%
% Inputs:
%   Tc        : temperature [C], scalar or vector
%   VWC       : total volumetric water content [m^3/m^3] at that depth
%   SWC       : porosity [m^3/m^3]
%   soilParams: struct, may include frac_clay, frac_silt, etc.
%   uwf       : struct with fields:
%       .mode = 'logistic' | 'temp_only' | 'texture'
%
% Notes:
%   - These are *phenomenological* UWF curves. They can be upgraded later
%     to include matric potential / Clapeyron-based soil freezing curves.

Tc  = Tc(:);
VWC = VWC(:);

% Defaults
if nargin < 5 || isempty(uwf), uwf = struct(); end
if ~isfield(uwf,'mode') || isempty(uwf.mode), uwf.mode = 'logistic'; end

% clamp helper
clamp01 = @(x) min(max(x,0),1);

switch lower(uwf.mode)

    case 'logistic'
        % --- Simple, tunable ---
        % fL = fmin + (1-fmin) / (1 + exp(-(T-T0)/dT))
        if ~isfield(uwf,'T0'),   uwf.T0 = 0.0; end     % C
        if ~isfield(uwf,'dT'),   uwf.dT = 0.5; end     % C
        if ~isfield(uwf,'fmin'), uwf.fmin = 0.02; end  % residual liquid fraction
        dT = max(uwf.dT, realmin);

        fL = uwf.fmin + (1-uwf.fmin) .* (1 ./ (1 + exp(-(Tc - uwf.T0)./dT)));

        info = uwf;

    case 'temp_only'
        % --- Fewer knobs: temperature-only power curve ---
        % Idea: residual unfrozen fraction decays with subzero temperature:
        % fL = max(fmin, exp(-gamma * max(0, Tm - T)))
        %
        % Parameters:
        %   Tm    : "melt point" (C), default 0
        %   gamma : 1/C, sets decay rate, default 1.2  (moderate)
        %   fmin  : residual
        if ~isfield(uwf,'Tm'),    uwf.Tm = 0.0; end
        if ~isfield(uwf,'gamma'), uwf.gamma = 1.2; end
        if ~isfield(uwf,'fmin'),  uwf.fmin = 0.02; end

        dT = max(0, uwf.Tm - Tc);           % only below freezing
        fL = exp(-uwf.gamma .* dT);
        fL = max(fL, uwf.fmin);
        fL = clamp01(fL);

        info = uwf;

    case 'texture'
        % --- Texture-based: Anderson/Koopmans-style power law ---
        % Common empirical form for unfrozen water content (volumetric):
        %   theta_u(T) = A * |T|^(-B)   for T<0
        % then fL = min(1, max(fmin, theta_u / VWC))
        %
        % We choose A,B from soil texture proxies (clay fraction).
        %
        % Parameters:
        %   frac_clay : from soilParams.frac_clay if present (0..1)
        %   A0,B0,A1,B1 : mapping coefficients (defaults below)
        %   fmin
        if ~isfield(uwf,'fmin'), uwf.fmin = 0.02; end

        frac_clay = 0.15;
        if nargin >= 4 && isstruct(soilParams) && isfield(soilParams,'frac_clay') && ~isempty(soilParams.frac_clay)
            frac_clay = soilParams.frac_clay;
        end
        frac_clay = min(max(frac_clay,0),1);

        % Default mapping: more clay => larger A and B => more unfrozen water at subzero
        if ~isfield(uwf,'A0'), uwf.A0 = 0.02; end
        if ~isfield(uwf,'A1'), uwf.A1 = 0.10; end
        if ~isfield(uwf,'B0'), uwf.B0 = 0.6;  end
        if ~isfield(uwf,'B1'), uwf.B1 = 1.0;  end

        A = uwf.A0 + uwf.A1 * frac_clay;     % volumetric
        B = uwf.B0 + uwf.B1 * frac_clay;     % exponent

        theta_u = zeros(size(Tc));
        subzero = (Tc < 0);

        % Avoid T=0 singularity with a small epsilon
        Tabs = abs(Tc(subzero)) + 1e-6;
        theta_u(subzero) = A .* (Tabs).^(-B);

        % Above 0C: all liquid
        theta_u(~subzero) = VWC(~subzero);

        % Convert to fraction of the *existing* pore water
        fL = ones(size(Tc));
        denom = max(VWC, realmin);
        fL(subzero) = theta_u(subzero) ./ denom(subzero);

        % Clamp
        fL = clamp01(fL);
        fL = max(fL, uwf.fmin);

        info = uwf;
        info.frac_clay_used = frac_clay;
        info.A_used = A;
        info.B_used = B;

    otherwise
        error('Unknown uwf.mode: %s', uwf.mode);
end
end

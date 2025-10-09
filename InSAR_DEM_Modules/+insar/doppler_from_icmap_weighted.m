function fdc = doppler_from_icmap_weighted(P, V, lambda, X, Y, Z, icMap, lookMask, varargin)
% Doppler centroid map from fractional closest-approach indices.
% Inputs:
%   P,V       [N x 3] platform pos/vel
%   lambda    wavelength [m]
%   X,Y,Z     DEM grids (same size)
%   icMap     fractional closest index per pixel (same size)
%   lookMask  logical mask (same size)
%   opts.windowSamples / windowMeters / windowTime (choose one)
%   opts.t (if using windowTime)
%   opts.useR2Weight (default true)
%
% Output:
%   fdc       (Hz), NaN where outside mask

    % arguments
    %     P double
    %     V double
    %     lambda (1,1) double
    %     X double
    %     Y double
    %     Z double
    %     icMap double
    %     lookMask logical
    %     opts.windowSamples (1,1) double = 51
    %     opts.windowMeters  = []
    %     opts.windowTime    = []
    %     opts.t = []
    %     opts.useR2Weight logical = true
    % end
% Accept NV pairs or a single struct as 9th arg
if ~isempty(varargin) && isstruct(varargin{1}) && numel(varargin)==1
    S = varargin{1};            % struct style
    varargin = reshape([fieldnames(S)'; struct2cell(S)'], 1, []); % -> NV pairs
end

% Defaults
p = inputParser;  p.FunctionName = mfilename;
addParameter(p,'windowSamples',51,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'windowMeters',[],@(x)isempty(x)||isscalar(x));
addParameter(p,'windowTime',[],@(x)isempty(x)||isscalar(x));
addParameter(p,'t',[],@(x)isempty(x)||isvector(x));
addParameter(p,'useR2Weight',true,@(x)islogical(x)||ismember(x,[0 1]));
parse(p, varargin{:});
opts = p.Results;

    [ny,nx] = size(X);
    fdc = nan(ny,nx);

    % Decide half-window in samples
    if ~isempty(opts.windowMeters)
        ds = vecnorm(diff(P),2,2); ds_med = median(ds(isfinite(ds))); if isempty(ds_med)||ds_med<=0, ds_med = 1; end
        halfW = max(3, round(opts.windowMeters / ds_med));
    elseif ~isempty(opts.windowTime) && ~isempty(opts.t)
        dt = median(diff(opts.t)); if ~isfinite(dt) || dt<=0, dt = 1; end
        halfW = max(3, round(opts.windowTime / dt));
    else
        halfW = max(3, round(opts.windowSamples));
    end
    sigma = max(1, halfW/2.5);   % Gaussian width

    % Process masked pixels row-by-row
    for iy = 1:ny
        cols = find(lookMask(iy,:));
        if isempty(cols), continue; end

        % fractional closest indices and linear weights
        ic = icMap(iy, cols);           % Kx1
        i0 = floor(ic);                 % integer base
        a  = ic - i0;                   % fractional part
        a = a(:);
        i1 = max(1, i0);
        i2 = min(size(P,1), i0+1);
        P1 = P(i1,:); P2 = P(i2,:);       % K×3
        % if numel(a) > 1
        %     keyboard
        % end
        P_ic = P1 + a.*(P2 - P1);
        % platform at ic: linear interpolation
        % P_ic = (1-a') .* P(i1,:) + a' .* P(i2,:);     % [K x 3]
        % LOS unit vectors at ic
        Gk = [X(iy,cols)', Y(iy,cols)', Z(iy,cols)'];
        R  = Gk-P_ic; Rn = sqrt(sum(R.^2,2)); Rhat = R ./ Rn;

        % windowed average of V around ic, then project onto LOS
        fLine = nan(numel(cols),1);

        for k = 1:numel(cols)
            c  = ic(k);
            c  = max(1, min(size(P,1), c));
            iL = max(1, floor(c) - halfW);
            iU = min(size(P,1), floor(c) + halfW);
            idx = (iL:iU).';
            u   = idx - c;

            w = exp(-0.5*(u./sigma).^2);
            if opts.useR2Weight
                w = w ./ max(Rn(k).^2, 1);
            end
            ws = sum(w);
            if ~isfinite(ws) || ws<=0, continue; end
            w = w / ws;

            Vbar = w' * V(idx,:);
            if any(~isfinite(Vbar)) || any(~isfinite(Rhat(k,:))), continue; end

            fLine(k) = (2/lambda) * dot(Vbar, Rhat(k,:));
        end

        fdc(iy, cols) = fLine;
    end
end

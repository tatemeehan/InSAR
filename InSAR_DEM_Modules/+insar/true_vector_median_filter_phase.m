function [phaseOut, zOut, scoreOut, Rlocal] = true_vector_median_filter_phase(phase, window_size, varargin)
%TRUE_VECTOR_MEDIAN_FILTER_PHASE Exact vector median filter for wrapped phase.
%
% This is a true vector median filter:
%
%   For each local window, choose the actual complex vector in the window
%   that minimizes the sum of Euclidean distances to all other vectors.
%
% For wrapped phase input, the phase is first mapped to the unit circle:
%
%   z = exp(1i * phase)
%
% Inputs
%   phase       : wrapped phase image [rad], or complex phasor/image
%   window_size : odd integer window size, e.g. 3, 5, 7, 11
%
% Name-value options
%   'Normalize'   : true/false, default true
%                   If true, complex input is normalized to unit phasors.
%                   Use true for phase filtering.
%
%   'UseParallel' : true/false, default false
%                   Use parfor over rows if Parallel Toolbox is available.
%
% Outputs
%   phaseOut : filtered wrapped phase [rad]
%   zOut     : selected vector median complex phasor/image
%   scoreOut : mean distance of selected vector to all vectors in window
%   Rlocal   : local circular/phasor concentration, abs(mean(z_window))
%
% Notes
%   This is exact but expensive. For an 11x11 window, each pixel solves
%   a 121 x 121 distance problem. Use on tiles first.

% -------------------------
% Parse inputs
% -------------------------
p = inputParser;
p.addParameter('Normalize', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('UseParallel', false, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});

normalizeInput = logical(p.Results.Normalize);
useParallel    = logical(p.Results.UseParallel);

if mod(window_size, 2) == 0
    error('window_size must be odd.');
end

halfWin = floor(window_size / 2);

% -------------------------
% Convert to complex vectors
% -------------------------
if isreal(phase)
    validInput = isfinite(phase);
    z = complex(nan(size(phase)), nan(size(phase)));
    z(validInput) = exp(1i .* double(phase(validInput)));
else
    z = complex(double(real(phase)), double(imag(phase)));
    validInput = isfinite(real(z)) & isfinite(imag(z));

    if normalizeInput
        mag = abs(z);
        good = validInput & mag > 0;
        z(good) = z(good) ./ mag(good);
        z(~good) = NaN;
    end
end

[nRows, nCols] = size(z);

% -------------------------
% Symmetric padding
% -------------------------
zPad = padarray(z, [halfWin halfWin], 'symmetric', 'both');

% -------------------------
% Allocate outputs
% -------------------------
zOut = complex(nan(nRows, nCols), nan(nRows, nCols));
scoreOut = nan(nRows, nCols);
Rlocal = nan(nRows, nCols);

% -------------------------
% Main filter
% -------------------------
if useParallel

    parfor rr = 1:nRows

        zRow = complex(nan(1, nCols), nan(1, nCols));
        scoreRow = nan(1, nCols);
        Rrow = nan(1, nCols);

        for cc = 1:nCols

            block = zPad(rr:rr+2*halfWin, cc:cc+2*halfWin);
            v = block(:);

            good = isfinite(real(v)) & isfinite(imag(v));
            v = v(good);

            if isempty(v)
                continue
            end

            if normalizeInput
                mag = abs(v);
                goodMag = mag > 0;
                v = v(goodMag) ./ mag(goodMag);
            end

            if isempty(v)
                continue
            end

            % True vector median: choose sample minimizing total distance
            D = abs(v - v.');
            sumD = sum(D, 2);

            [bestScore, idxBest] = min(sumD);

            zRow(cc) = v(idxBest);
            scoreRow(cc) = bestScore ./ numel(v);
            Rrow(cc) = abs(mean(v, 'omitnan'));
        end

        zOut(rr, :) = zRow;
        scoreOut(rr, :) = scoreRow;
        Rlocal(rr, :) = Rrow;
    end

else

    for rr = 1:nRows
        for cc = 1:nCols

            block = zPad(rr:rr+2*halfWin, cc:cc+2*halfWin);
            v = block(:);

            good = isfinite(real(v)) & isfinite(imag(v));
            v = v(good);

            if isempty(v)
                continue
            end

            if normalizeInput
                mag = abs(v);
                goodMag = mag > 0;
                v = v(goodMag) ./ mag(goodMag);
            end

            if isempty(v)
                continue
            end

            % True vector median: choose sample minimizing total distance
            D = abs(v - v.');
            sumD = sum(D, 2);

            [bestScore, idxBest] = min(sumD);

            zOut(rr, cc) = v(idxBest);
            scoreOut(rr, cc) = bestScore ./ numel(v);
            Rlocal(rr, cc) = abs(mean(v, 'omitnan'));
        end
    end
end

phaseOut = angle(zOut);

end
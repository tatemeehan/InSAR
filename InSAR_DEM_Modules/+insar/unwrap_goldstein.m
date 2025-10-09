function unwrapped = unwrap_goldstein(wrappedPhase, alpha, windowSize)
%UNWRAP_GOLDSTEIN Applies Goldstein's branch-cut phase unwrapping algorithm.
%
%   unwrapped = unwrap_goldstein(wrappedPhase)
%   unwrapped = unwrap_goldstein(wrappedPhase, alpha, windowSize)
%
%   Inputs:
%       wrappedPhase : 2D wrapped phase image (radians)
%       alpha         : High-pass filter exponent (default = 0.5)
%       windowSize    : Filter window size (default = 5)
%
%   Output:
%       unwrapped     : 2D unwrapped phase image (radians)

if nargin < 2 || isempty(alpha), alpha = 0.5; end
if nargin < 3 || isempty(windowSize), windowSize = 5; end

% Compute phase quality map using Goldstein's method
quality = compute_goldstein_quality(wrappedPhase, alpha, windowSize);

% Identify residues (2pi discontinuities)
residues = compute_phase_residues(wrappedPhase);

% Place branch cuts to connect residues
branchCuts = place_branch_cuts(residues, quality);

% Mask out branch cuts
wrappedMasked = wrappedPhase;
wrappedMasked(branchCuts) = NaN;

% Apply flood-fill style unwrapping avoiding branch cuts
unwrapped = flood_fill_unwrap(wrappedMasked, quality);

end

function quality = compute_goldstein_quality(phase, alpha, windowSize)
    % High-pass filter the wrapped phase gradient magnitude
    dx = wrapToPi([diff(phase,1,2), zeros(size(phase,1),1)]);
    dy = wrapToPi([diff(phase,1,1); zeros(1,size(phase,2))]);
    gradMag = sqrt(dx.^2 + dy.^2);

    % High-pass filter: subtract boxcar low-pass filtered version
    h = fspecial('average', windowSize);
    lowpass = imfilter(gradMag, h, 'replicate');
    highpass = abs(gradMag - lowpass);

    quality = highpass .^ alpha;
end

function residues = compute_phase_residues(phase)
    % Compute 2D phase residues (4-pixel loops)
    [M, N] = size(phase);
    residues = zeros(M-1, N-1);
    for i = 1:M-1
        for j = 1:N-1
            loop = [phase(i,j), phase(i,j+1), phase(i+1,j+1), phase(i+1,j)];
            dphi = diff(unwrap([loop loop(1)]));
            residues(i,j) = round(sum(wrapToPi(dphi))/(2*pi));
        end
    end
end

function cuts = place_branch_cuts(residues, quality)
    % Naive branch cut placement: connect +1 and -1 residues by quality map descent
    cuts = false(size(residues) + 1);
    [posY, posX] = find(residues == 1);
    [negY, negX] = find(residues == -1);

    for k = 1:min(numel(posX), numel(negX))
        p1 = [posY(k), posX(k)];
        p2 = [negY(k), negX(k)];
        path = bresenham(p1, p2);
        for idx = 1:size(path,1)
            y = path(idx,1); x = path(idx,2);
            cuts(y,x) = true;
        end
    end
end

function unwrapped = flood_fill_unwrap(wrapped, quality)
    % Region growing unwrapping avoiding branch cuts
    unwrapped = nan(size(wrapped));
    visited = false(size(wrapped));
    [~, seedIdx] = max(quality(:));
    [sy, sx] = ind2sub(size(wrapped), seedIdx);

    queue = [sy sx];
    unwrapped(sy,sx) = wrapped(sy,sx);
    visited(sy,sx) = true;

    while ~isempty(queue)
        [y, x] = deal(queue(1,1), queue(1,2));
        queue(1,:) = [];

        neighbors = [y-1,x; y+1,x; y,x-1; y,x+1];
        for n = 1:4
            ny = neighbors(n,1);
            nx = neighbors(n,2);
            if ny < 1 || ny > size(wrapped,1) || nx < 1 || nx > size(wrapped,2)
                continue;
            end
            if visited(ny,nx) || isnan(wrapped(ny,nx))
                continue;
            end
            dphi = wrapToPi(wrapped(ny,nx) - wrapped(y,x));
            unwrapped(ny,nx) = unwrapped(y,x) + dphi;
            visited(ny,nx) = true;
            queue(end+1,:) = [ny nx];
        end
    end
end

function path = bresenham(p1, p2)
    % Bresenham's line algorithm for integer grid paths
    x1 = p1(2); y1 = p1(1);
    x2 = p2(2); y2 = p2(1);
    dx = abs(x2 - x1); sx = sign(x2 - x1);
    dy = -abs(y2 - y1); sy = sign(y2 - y1);
    err = dx + dy;
    path = [];

    while true
        path(end+1,:) = [y1, x1];
        if x1 == x2 && y1 == y2, break; end
        e2 = 2*err;
        if e2 >= dy
            err = err + dy;
            x1 = x1 + sx;
        end
        if e2 <= dx
            err = err + dx;
            y1 = y1 + sy;
        end
    end
end

function out = optimize_antenna_angle(flatIncidenceDeg, slantRangeM, lookmask, varargin)
%OPTIMIZE_ANTENNA_TILT_GAUSSIAN
% Choose antenna elevation tilt (0=horizontal, -90=nadir) to maximize coverage
% using a physically realistic Gaussian-like beam taper.
%
% Beam model:
%   G(d) = exp(-(d/alpha)^p)  where d = |beta - beta0| in degrees
% alpha and p are fit so that:
%   G(20 deg) = 10^(-6/20)   (=-6 dB)
%   G(30 deg) = 10^(-12/20) (=-12 dB)
%
% Inputs:
%   flatIncidenceDeg  : LOS elevation above horizontal ("flat-earth incidence"), deg (same size as lookmask)
%   slantRangeM       : slant range, meters (same size)
%   lookmask          : logical mask of valid sidelooking pixels (same size)
%
% Name-Value options:
%   'CR'              : optional CR table/struct with incidence column/field (deg) named
%                       'incidence' (case-insensitive). Can be [].
%   'BeamSpec'        : [d1 dB1 d2 dB2] default [20 -6 30 -12]
%   'BetaStepDeg'     : grid step for beta center search, default 0.5
%   'BetaCentersDeg'  : explicit vector of centers (deg). If provided overrides step/range.
%   'BetaRangeDeg'    : [min max] range for center search, default [0 90]
%   'WeightMode'      : 'none' | 'R2' | 'R4' | 'R2clamp' | 'R4clamp', default 'none'
%   'ClampPrctile'    : percentile for clamp modes (e.g., 10), default 10
%   'MinCRGain'       : minimum allowed beam gain at each CR (0..1). default [] (no constraint)
%   'TiltStepDeg'     : discrete antenna step (deg), default 5
%   'PreferLessNadir' : tie-breaker for discrete tilt: true -> prefer less nadir (closer to 0),
%                       false -> prefer more nadir. default true
%
% Output struct 'out' fields:
%   .betaOptDeg_cont         continuous optimal beta center (deg)
%   .tiltOptDeg_cont         continuous optimal tilt (deg) = -betaOptDeg_cont
%   .betaOptDeg_discrete     snapped beta center (deg)
%   .tiltOptDeg_discrete     snapped tilt (deg)
%   .alpha, .p               beam parameters
%   .beamSpec                spec used
%   .scoreCenters            score for each beta center
%   .betaCentersDeg          centers evaluated
%   .crBetaDeg               CR betas used (deg) (if provided)
%   .crGain_cont, .crGain_discrete gains at CRs for chosen solutions (if CR provided)
%   .notes                   text notes
%
% Notes:
% - flatIncidenceDeg is the correct geometry variable for antenna elevation pointing.
% - This function optimizes *uniform coverage* by default using mean beam gain over ROI.
% “Antenna elevation was optimized by maximizing mean beam-weighted footprint coverage in flat-earth LOS elevation (β),
% using a generalized Gaussian beam model fit to −6 dB at ±20° and −12 dB at ±30°, 
% with corner reflectors required to receive ≥ −12 dB relative gain.”

% ---------------- Parse options ----------------
p = inputParser;
p.addParameter('CR', [], @(x) true);
p.addParameter('BeamSpec', [20 -6 30 -12], @(x) isnumeric(x) && numel(x)==4);
p.addParameter('BetaStepDeg', 0.5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('BetaCentersDeg', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('BetaRangeDeg', [0 90], @(x) isnumeric(x) && numel(x)==2);
p.addParameter('WeightMode', 'none', @(x) ischar(x) || isstring(x));
p.addParameter('ClampPrctile', 10, @(x) isnumeric(x) && isscalar(x) && x>0 && x<100);
p.addParameter('MinCRGain', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>=0 && x<=1));
p.addParameter('TiltStepDeg', 5, @(x) isnumeric(x) && isscalar(x) && x>0);
p.addParameter('PreferLessNadir', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

% ---------------- Validate sizes ----------------
assert(isequal(size(flatIncidenceDeg), size(slantRangeM), size(lookmask)), ...
    'flatIncidenceDeg, slantRangeM, and lookmask must be the same size.');

% ---------------- Vectorize ROI ----------------
roiMask = lookmask & isfinite(flatIncidenceDeg) & isfinite(slantRangeM);
beta = flatIncidenceDeg(roiMask);  beta = beta(:);
R    = slantRangeM(roiMask);       R    = R(:);

% ---------------- Weights ----------------
w = compute_weights(R, opts.WeightMode, opts.ClampPrctile);

% Normalize weights for stability (optional; doesn’t change argmax for mean gain if constant factor)
w = w ./ (mean(w) + eps);

% ---------------- Beam fit (generalized Gaussian) ----------------
beamSpec = opts.BeamSpec(:).';
d1 = beamSpec(1); dB1 = beamSpec(2);
d2 = beamSpec(3); dB2 = beamSpec(4);

G1 = 10^(dB1/20);   % amplitude gain (dB is negative)
G2 = 10^(dB2/20);

% Fit p and alpha for G(d)=exp(-(d/alpha)^p)
% lnG = -(d/alpha)^p
L1 = -log(G1);  % positive
L2 = -log(G2);  % positive
p_exp = log(L2/L1) / log(d2/d1);
alpha = d1 / (L1)^(1/p_exp);

% ---------------- Candidate centers ----------------
if ~isempty(opts.BetaCentersDeg)
    betaCenters = opts.BetaCentersDeg(:);
else
    br = opts.BetaRangeDeg;
    betaCenters = (br(1):opts.BetaStepDeg:br(2)).';
end

% ---------------- Optional CR betas ----------------
crBeta = [];
if ~isempty(opts.CR)
    crBeta = extract_cr_beta(opts.CR);  % column vector deg
end

% ---------------- Score each center ----------------
score = -inf(size(betaCenters));

for k = 1:numel(betaCenters)
    b0 = betaCenters(k);

    % CR constraint (optional)
    if ~isempty(crBeta) && ~isempty(opts.MinCRGain)
        Gcr = beam_gain(abs(crBeta - b0), alpha, p_exp);
        if any(Gcr < opts.MinCRGain)
            continue
        end
    end

    % Mean beam gain over ROI (cartographer objective)
    d = abs(beta - b0);
    G = beam_gain(d, alpha, p_exp);

    % Weighted mean gain (w defaults to uniform if WeightMode='none')
    score(k) = sum(w .* G) / sum(w);
end

% ---------------- Pick continuous optimum ----------------
[~, ix] = max(score);
betaOpt_cont = betaCenters(ix);
tiltOpt_cont = -betaOpt_cont;

% ---------------- Snap to discrete tilt increments ----------------
tiltCandidates = snap_candidates(tiltOpt_cont, opts.TiltStepDeg);
% Evaluate candidates by same score rule (and same CR constraint)
[tiltOpt_disc, betaOpt_disc] = choose_best_discrete(tiltCandidates, beta, w, crBeta, opts.MinCRGain, alpha, p_exp, opts.PreferLessNadir);

% ---------------- Report CR gains ----------------
crGain_cont = [];
crGain_disc = [];
if ~isempty(crBeta)
    crGain_cont = beam_gain(abs(crBeta - betaOpt_cont), alpha, p_exp);
    crGain_disc = beam_gain(abs(crBeta - betaOpt_disc), alpha, p_exp);
end

% ---------------- Package output ----------------
out = struct();
out.betaOptDeg_cont     = betaOpt_cont;
out.tiltOptDeg_cont     = tiltOpt_cont;
out.betaOptDeg_discrete = betaOpt_disc;
out.tiltOptDeg_discrete = tiltOpt_disc;
out.alpha               = alpha;
out.p                   = p_exp;
out.beamSpec            = beamSpec;
out.betaCentersDeg      = betaCenters;
out.scoreCenters        = score;
out.crBetaDeg           = crBeta;
out.crGain_cont         = crGain_cont;
out.crGain_discrete     = crGain_disc;

out.notes = sprintf(['Beam model: G(d)=exp(-(d/alpha)^p), fit to %gdB@%gdeg and %gdB@%gdeg.\n' ...
                     'Objective: maximize weighted mean gain over ROI (cartographer). WeightMode=%s.\n'], ...
                     dB1, d1, dB2, d2, string(opts.WeightMode));
end

% ========================= Helper functions =========================

function w = compute_weights(R, mode, clampPrct)
mode = lower(string(mode));
switch mode
    case "none"
        w = ones(size(R));
    case "r2"
        w = 1 ./ max(R,1).^2;
    case "r4"
        w = 1 ./ max(R,1).^4;
    case "r2clamp"
        Rmin = prctile(R, clampPrct);
        w = 1 ./ max(R, Rmin).^2;
    case "r4clamp"
        Rmin = prctile(R, clampPrct);
        w = 1 ./ max(R, Rmin).^4;
    otherwise
        error("Unknown WeightMode: %s", mode);
end
end

function G = beam_gain(dDeg, alpha, p_exp)
% Generalized Gaussian taper in angle domain
G = exp(-(dDeg./alpha).^p_exp);
end

function crBeta = extract_cr_beta(CR)
% Accepts:
%  - table with variable 'incidence' (case-insensitive)
%  - struct with field 'incidence' or 'Incidence'
%  - numeric vector (assumed already beta/flatInc degrees)
if istable(CR)
    vnames = CR.Properties.VariableNames;
    idx = find(strcmpi(vnames, 'losElevationAngle'), 1);
    if isempty(idx)
        error("CR table must contain a column named 'losElevationAngle' (case-insensitive).");
    end
    crBeta = CR.(vnames{idx});
elseif isstruct(CR)
    f = fieldnames(CR);
    idx = find(strcmpi(f, 'losElevationAngle'), 1);
    if isempty(idx)
        error("CR struct must contain a field named 'losElevationAngle' (case-insensitive).");
    end
    crBeta = CR.(f{idx});
elseif isnumeric(CR)
    crBeta = CR;
else
    error("Unsupported CR input type. Provide a table, struct, or numeric vector.");
end

crBeta = crBeta(:);
crBeta = crBeta(isfinite(crBeta));
end

function tilts = snap_candidates(tiltContDeg, stepDeg)
% Provide the two nearest step candidates (and clip to [-90,0] range)
t0 = stepDeg * floor(tiltContDeg/stepDeg);
t1 = t0 + stepDeg;
tilts = unique([t0; t1]);

% Ensure within feasible mechanical range
tilts = max(min(tilts, 0), -90);
tilts = unique(tilts);
end

function [tiltBest, betaBest] = choose_best_discrete(tiltCandidates, beta, w, crBeta, minCRGain, alpha, p_exp, preferLessNadir)
% Evaluate a small set of discrete tilts; pick the best score, with CR constraint if provided.
bestScore = -inf;
tiltBest = tiltCandidates(1);
betaBest = -tiltBest;

for i = 1:numel(tiltCandidates)
    tilt = tiltCandidates(i);
    b0 = -tilt; % beta center

    % CR constraint (optional)
    if ~isempty(crBeta) && ~isempty(minCRGain)
        Gcr = beam_gain(abs(crBeta - b0), alpha, p_exp);
        if any(Gcr < minCRGain)
            continue
        end
    end

    d = abs(beta - b0);
    G = beam_gain(d, alpha, p_exp);
    s = sum(w .* G) / sum(w);

    if s > bestScore + 1e-12
        bestScore = s;
        tiltBest = tilt;
        betaBest = b0;
    elseif abs(s - bestScore) <= 1e-12
        % Tie-break: prefer less nadir (tilt closer to 0) or more nadir
        if preferLessNadir
            if tilt > tiltBest
                tiltBest = tilt;
                betaBest = b0;
            end
        else
            if tilt < tiltBest
                tiltBest = tilt;
                betaBest = b0;
            end
        end
    end
end
end

% function [tiltOptDeg, betaOptDeg, score, betaCenters] = optimize_antenna_angle(flatIncDeg, rSlant, lookmask, beamHalfWidthDeg, betaStepDeg)
% % flatIncDeg: LOS elevation above horizontal (your flat-earth incidence), degrees
% % rSlant:     slant range (m)
% % lookmask:   logical mask
% % beamHalfWidthDeg: e.g., 20
% % betaStepDeg: e.g., 0.5 or 1
% 
% if nargin < 4 || isempty(beamHalfWidthDeg), beamHalfWidthDeg = 20; end
% if nargin < 5 || isempty(betaStepDeg), betaStepDeg = 0.5; end
% 
% % Vectorize valid samples
% beta = flatIncDeg(lookmask);
% R    = rSlant(lookmask);
% 
% beta = beta(:);
% R    = R(:);
% 
% % Basic two-way spreading weight
% % w = 1 ./ max(R,1).^2;
% w = ones(numel(R),1);
% % Rmin = prctile(R, 10);
% % w = 1 ./ max(R, Rmin).^2;
% 
% % Remove NaNs
% ok = isfinite(beta) & isfinite(w);
% beta = beta(ok);
% w    = w(ok);
% 
% % Search centers
% betaCenters = 0:betaStepDeg:90;
% score = zeros(size(betaCenters));
% 
% for k = 1:numel(betaCenters)
%     b0 = betaCenters(k);
%     inBeam = abs(beta - b0) <= beamHalfWidthDeg;
%     score(k) = sum(w(inBeam));
% end
% 
% % Best center
% [~,ix] = max(score);
% betaOptDeg = betaCenters(ix);
% 
% % Your antenna convention: 0 = horizontal, -90 = nadir
% tiltOptDeg = -betaOptDeg;
% end

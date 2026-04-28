function out = sweep_initial_lwc_phase(theta_i, f, snow1_base, snow2_base, eps_g0, pars1, pars2, lwc1_vec, chan)
%SWEEP_INITIAL_LWC_PHASE
% Sweep pass-1 LWC while holding pass-2 fixed, and diagnose wrapped/unwrapped
% interferometric phase for total and component terms.
%
% Inputs:
%   theta_i   : incidence angle [deg]
%   f         : frequency [Hz]
%   snow1_base: struct for pass 1 snow state (fields Hs,rho,lwc,Tk)
%   snow2_base: struct for pass 2 snow state
%   eps_g0    : soil/background input for run_pass
%   pars1     : parameter struct for pass 1
%   pars2     : parameter struct for pass 2
%   lwc1_vec  : vector of pass-1 LWC values to test (fraction, e.g. 0:0.0005:0.02)
%   chan      : 'HH','VV','HV','VH'   (default 'VV')
%
% Output:
%   out: struct with wrapped/unwrapped phase curves, component magnitudes,
%        wrap indices, and saved model outputs.

if nargin < 9 || isempty(chan)
    chan = 'VV';
end
chan = upper(string(chan));

% Force matrix mode so indexing is always consistent
pars1.pol_mode = "matrix";
pars2.pol_mode = "matrix";

% Channel index
[ii,jj] = pol_to_idx(chan);

n = numel(lwc1_vec);

% Preallocate
Gtot = complex(n,1);
Gts  = complex(n,1);
Gg   = complex(n,1);
Gas  = complex(n,1);
Gsv  = complex(n,1);
Gsg  = complex(n,1);
Ggv  = complex(n,1);

E1mag = zeros(n,1);
E2mag = zeros(n,1);
Eas1mag = zeros(n,1);
Es1mag  = zeros(n,1);
Eg1mag  = zeros(n,1);
Esg1mag = zeros(n,1);
Egv1mag = zeros(n,1);

o1save = cell(n,1);
o2save = cell(n,1);

for k = 1:n
    snow1 = snow1_base;
    snow2 = snow2_base;

    snow1.lwc = lwc1_vec(k);

    p1 = pars1; p2 = pars2;
    p1.snowState = snow1;
    p2.snowState = snow2;

    o1 = run_pass(theta_i, f, snow1, eps_g0, p1);
    o2 = run_pass(theta_i, f, snow2, eps_g0, p2);

    % Save if you want to inspect particular cases later
    o1save{k} = o1;
    o2save{k} = o2;

    % Extract channel
    E1   = get_chan(o1.E,   ii, jj);
    E2   = get_chan(o2.E,   ii, jj);
    Eas1 = get_chan(o1.Eas, ii, jj);
    Eas2 = get_chan(o2.Eas, ii, jj);
    Es1  = get_chan(o1.Es,  ii, jj);
    Es2  = get_chan(o2.Es,  ii, jj);
    Esg1 = get_chan(o1.Esg, ii, jj);
    Esg2 = get_chan(o2.Esg, ii, jj);
    Egv1 = get_chan(o1.Egv, ii, jj);
    Egv2 = get_chan(o2.Egv, ii, jj);

    Eg1 = Esg1 + Egv1;
    Eg2 = Esg2 + Egv2;

    % Interferometric phasors
    Gtot(k) = E2 * conj(E1);
    Gg(k)   = Eg2 * conj(Eg1);
    Gas(k)  = Eas2 * conj(Eas1);
    Gsv(k)  = Es2 * conj(Es1);
    Gsg(k)  = Esg2 * conj(Esg1);
    Ggv(k)  = Egv2 * conj(Egv1);

    % Ts_eff or fallback to Ts
    if isfield(o1,'Ts_eff') && isfield(o2,'Ts_eff')
        Ts1 = o1.Ts_eff;
        Ts2 = o2.Ts_eff;
    else
        Ts1 = o1.Ts;
        Ts2 = o2.Ts;
    end
    Gts(k) = Ts2 * conj(Ts1);

    % Magnitudes for cancellation diagnostics
    E1mag(k)   = abs(E1);
    E2mag(k)   = abs(E2);
    Eas1mag(k) = abs(Eas1);
    Es1mag(k)  = abs(Es1);
    Eg1mag(k)  = abs(Eg1);
    Esg1mag(k) = abs(Esg1);
    Egv1mag(k) = abs(Egv1);
end

% Wrapped phases
ph_tot = angle(Gtot);
ph_ts  = angle(Gts);
ph_g   = angle(Gg);
ph_as  = angle(Gas);
ph_sv  = angle(Gsv);
ph_sg  = angle(Gsg);
ph_gv  = angle(Ggv);

% Unwrapped phases
phu_tot = unwrap(ph_tot);
phu_ts  = unwrap(ph_ts);
phu_g   = unwrap(ph_g);
phu_as  = unwrap(ph_as);
phu_sv  = unwrap(ph_sv);
phu_sg  = unwrap(ph_sg);
phu_gv  = unwrap(ph_gv);

% Wrap detection
wrap_idx_tot = find(abs(diff(ph_tot)) > pi);
wrap_idx_ts  = find(abs(diff(ph_ts))  > pi);
wrap_idx_g   = find(abs(diff(ph_g))   > pi);

wrap_lwc_tot = 0.5*(lwc1_vec(wrap_idx_tot) + lwc1_vec(wrap_idx_tot+1));
wrap_lwc_ts  = 0.5*(lwc1_vec(wrap_idx_ts)  + lwc1_vec(wrap_idx_ts+1));
wrap_lwc_g   = 0.5*(lwc1_vec(wrap_idx_g)   + lwc1_vec(wrap_idx_g+1));

% Likely cancellation locations
[~, idx_min_E1] = min(E1mag);
[~, idx_min_E2] = min(E2mag);

out = struct();
out.chan = chan;
out.lwc1 = lwc1_vec(:);

out.G = struct('tot',Gtot,'ts',Gts,'g',Gg,'as',Gas,'sv',Gsv,'sg',Gsg,'gv',Ggv);
out.ph = struct('tot',ph_tot,'ts',ph_ts,'g',ph_g,'as',ph_as,'sv',ph_sv,'sg',ph_sg,'gv',ph_gv);
out.phu = struct('tot',phu_tot,'ts',phu_ts,'g',phu_g,'as',phu_as,'sv',phu_sv,'sg',phu_sg,'gv',phu_gv);

out.mag = struct('E1',E1mag,'E2',E2mag,'Eas1',Eas1mag,'Es1',Es1mag,'Eg1',Eg1mag,'Esg1',Esg1mag,'Egv1',Egv1mag);

out.wrap = struct( ...
    'idx_tot',wrap_idx_tot, 'lwc_tot',wrap_lwc_tot, ...
    'idx_ts', wrap_idx_ts,  'lwc_ts', wrap_lwc_ts, ...
    'idx_g',  wrap_idx_g,   'lwc_g',  wrap_lwc_g);

out.minmag = struct( ...
    'idx_E1', idx_min_E1, 'lwc_E1', lwc1_vec(idx_min_E1), 'E1min', E1mag(idx_min_E1), ...
    'idx_E2', idx_min_E2, 'lwc_E2', lwc1_vec(idx_min_E2), 'E2min', E2mag(idx_min_E2));

out.o1 = o1save;
out.o2 = o2save;
end

function z = get_chan(X, ii, jj)
if isempty(X)
    z = 0;
elseif isscalar(X)
    z = X;
else
    z = X(ii,jj);
end
end

function [ii,jj] = pol_to_idx(chan)
switch upper(string(chan))
    case "HH"
        ii = 1; jj = 1;
    case "HV"
        ii = 1; jj = 2;
    case "VH"
        ii = 2; jj = 1;
    case "VV"
        ii = 2; jj = 2;
    otherwise
        error('Unknown channel: %s', chan);
end
end
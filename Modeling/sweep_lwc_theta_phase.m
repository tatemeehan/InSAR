function out = sweep_lwc_theta_phase(theta_vec, lwc1_vec, f, snow1_base, snow2_base, ...
    eps_g0, pars1, pars2, obs_phase_ref, chan)
%SWEEP_LWC_THETA_PHASE
% 2D sweep over incidence angle and pass-1 initial LWC.
%
% Convention:
% - Use output_phase_convention = "data_matched"
% - Interferometric phase formed as angle(E_late .* conj(E_early))
%
% Inputs
%   theta_vec      : vector of incidence angles [deg]
%   lwc1_vec       : vector of initial pass-1 LWC values [fraction]
%   f              : frequency [Hz]
%   snow1_base     : pass-1 snow struct template
%   snow2_base     : pass-2 snow struct template
%   eps_g0         : soil/background arg for run_pass
%   pars1, pars2   : parameter structs
%   obs_phase_ref  : observed referenced phase [rad], scalar or []
%   chan           : 'HH','VV','HV','VH' (default 'VV')
%
% Outputs
%   out struct with wrapped/unwrapped total and components

if nargin < 9 || isempty(obs_phase_ref)
    obs_phase_ref = [];
end
if nargin < 10 || isempty(chan)
    chan = 'VV';
end
chan = upper(string(chan));

[ii,jj] = pol_to_idx(chan);

nT = numel(theta_vec);
nL = numel(lwc1_vec);

% Force data-matched output convention and matrix mode
pars1.pol_mode = "matrix";
pars2.pol_mode = "matrix";
pars1.output_phase_convention = "data_matched";
pars2.output_phase_convention = "data_matched";

% Preallocate interferometric phasors
Gtot = complex(zeros(nT,nL));
Gts  = complex(zeros(nT,nL));
Gg   = complex(zeros(nT,nL));
Gas  = complex(zeros(nT,nL));
Gsv  = complex(zeros(nT,nL));

magE1 = zeros(nT,nL);
magE2 = zeros(nT,nL);
magEg1 = zeros(nT,nL);
magEas1 = zeros(nT,nL);
magEs1 = zeros(nT,nL);

for it = 1:nT
    theta_i = theta_vec(it);

    for il = 1:nL
        snow1 = snow1_base;
        snow2 = snow2_base;

        snow1.lwc = lwc1_vec(il);

        p1 = pars1; p2 = pars2;
        p1.snowState = snow1;
        p2.snowState = snow2;

        % EARLY then LATE, always
        o1 = run_pass(theta_i, f, snow1, eps_g0, p1);   % early
        o2 = run_pass(theta_i, f, snow2, eps_g0, p2);   % late

        E1   = get_chan(o1.E,   ii, jj);
        E2   = get_chan(o2.E,   ii, jj);
        Eas1 = get_chan(o1.Eas, ii, jj);
        Eas2 = get_chan(o2.Eas, ii, jj);
        Es1  = get_chan(o1.Es,  ii, jj);
        Es2  = get_chan(o2.Es,  ii, jj);

        Eg1 = get_chan(o1.Esg + o1.Egv, ii, jj);
        Eg2 = get_chan(o2.Esg + o2.Egv, ii, jj);

        % Data-matched interferogram convention:
        % positive phase for increased later-pass snow delay
        Gtot(it,il) = E2 * conj(E1);
        Gts(it,il)  = o2.Ts_excess * conj(o1.Ts_excess);
        Gg(it,il)   = Eg2 * conj(Eg1);
        Gas(it,il)  = Eas2 * conj(Eas1);
        Gsv(it,il)  = Es2 * conj(Es1);

        magE1(it,il)   = abs(E1);
        magE2(it,il)   = abs(E2);
        magEg1(it,il)  = abs(Eg1);
        magEas1(it,il) = abs(Eas1);
        magEs1(it,il)  = abs(Es1);
    end
end

% Wrapped phases
ph_wrap.tot = angle(Gtot);
ph_wrap.ts  = angle(Gts);
ph_wrap.g   = angle(Gg);
ph_wrap.as  = angle(Gas);
ph_wrap.sv  = angle(Gsv);

% Unwrap along LWC dimension for each incidence angle
ph_unwrap.tot = unwrap(ph_wrap.tot, [], 2);
ph_unwrap.ts  = unwrap(ph_wrap.ts,  [], 2);
ph_unwrap.g   = unwrap(ph_wrap.g,   [], 2);
ph_unwrap.as  = unwrap(ph_wrap.as,  [], 2);
ph_unwrap.sv  = unwrap(ph_wrap.sv,  [], 2);

% Wrap count / branch index relative to wrapped phase
wrap_index = round((ph_unwrap.tot - ph_wrap.tot)/(2*pi));

% Wrapped misfit to observed referenced phase
misfit_wrap = [];
if ~isempty(obs_phase_ref)
    misfit_wrap = angle(exp(1i*(ph_wrap.tot - obs_phase_ref)));
end

out = struct();
out.theta = theta_vec(:);
out.lwc1  = lwc1_vec(:);

out.G = struct('tot',Gtot,'ts',Gts,'g',Gg,'as',Gas,'sv',Gsv);
out.ph_wrap = ph_wrap;
out.ph_unwrap = ph_unwrap;
out.wrap_index = wrap_index;

out.mag = struct( ...
    'E1',magE1, 'E2',magE2, ...
    'Eg1',magEg1, 'Eas1',magEas1, 'Es1',magEs1);

out.misfit_wrap = misfit_wrap;
out.obs_phase_ref = obs_phase_ref;
out.chan = chan;
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
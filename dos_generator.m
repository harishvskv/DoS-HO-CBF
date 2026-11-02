function A = dos_generator(p, t, mode, opt)
%DOS_GENERATOR  Link availability, packet arrival and age-of-data sequences for an outage pattern.

if nargin < 3 || isempty(mode), mode = 'none'; end
if nargin < 4, opt = struct(); end
opt = setdef(opt,'seed',1);
opt = setdef(opt,'starts',120);
opt = setdef(opt,'durations',3.0);
opt = setdef(opt,'n_range',[1 3]);
opt = setdef(opt,'start_range',[115 125]);
opt = setdef(opt,'dur_range',[2.0 3.4]);
opt = setdef(opt,'spike_range',[0.5 2.0]);
opt = setdef(opt,'jitter','uniform');

t = t(:)';
N = numel(t);
dt = t(2) - t(1);

age_healthy = p.cyber.T_pkt + p.cyber.tau_nom_hi;
assert(age_healthy < p.etm.tau_trigger, ...
    ['dos_generator: worst healthy age of data %.4f s is not below the ' ...
     'event trigger %.4f s. The trigger would fire on normal jitter.'], ...
     age_healthy, p.etm.tau_trigger);

draws = struct('mode',mode,'seed',opt.seed);
switch lower(mode)
    case 'none'
        W = zeros(0,2);
        draws.spec = 'no attack';

    case 'deterministic'
        s0 = opt.starts(:);  d0 = opt.durations(:);
        if isscalar(d0), d0 = repmat(d0,size(s0)); end
        W = [s0, s0 + d0];
        draws.spec = 'deterministic start times and durations';

    case 'stochastic'
        rng(opt.seed,'twister');
        nA = randi(opt.n_range);
        s0 = opt.start_range(1) + diff(opt.start_range)*rand(nA,1);
        d0 = opt.dur_range(1)   + diff(opt.dur_range)*rand(nA,1);
        W  = [s0, s0 + d0];
        draws.spec = struct( ...
            'count',    sprintf('uniform integer on [%d, %d]', opt.n_range), ...
            'start',    sprintf('uniform on [%.2f, %.2f] s', opt.start_range), ...
            'duration', sprintf('uniform on [%.2f, %.2f] s', opt.dur_range), ...
            'spike',    sprintf('uniform on [%.2f, %.2f] s', opt.spike_range));
        draws.count = nA;

    otherwise
        error('dos_generator:mode','Unknown mode "%s"', mode);
end
draws.windows = W;

alpha = true(1,N);
for k = 1:size(W,1)
    alpha(t >= W(k,1) & t <= W(k,2)) = false;
end

switch lower(mode)
    case 'stochastic'
        tau_net = p.cyber.tau_nom_lo + ...
                  (p.cyber.tau_nom_hi - p.cyber.tau_nom_lo)*rand(1,N);
        spike = opt.spike_range(1) + diff(opt.spike_range)*rand(1,N);
        tau_net(~alpha) = tau_net(~alpha) + spike(~alpha);
    otherwise
        mid = 0.5*(p.cyber.tau_nom_lo + p.cyber.tau_nom_hi);
        amp = 0.5*(p.cyber.tau_nom_hi - p.cyber.tau_nom_lo);
        tau_net = mid + amp*sin(0.5*t);
end
tau_net = min(tau_net, opt.spike_range(2) + p.cyber.tau_nom_hi);

pkt     = false(1,N);
t_stamp = nan(1,N);
kGen    = round(p.cyber.T_pkt/dt);
kGen    = max(kGen,1);

last_stamp = t(1);
for n = 1:N
    if mod(n-1, kGen) == 0
        tg = t(n);
        na = n + round(tau_net(n)/dt);
        if na <= N && all(alpha(n:min(na,N)))
            pkt(na) = true;
            t_stamp(na) = tg;
        end
    end
end
for n = 1:N
    if pkt(n), last_stamp = t_stamp(n); end
    t_stamp(n) = last_stamp;
end
age = t - t_stamp;

A = struct('t',t,'alpha',alpha,'tau_net',tau_net,'pkt',pkt, ...
           't_stamp',t_stamp,'age',age,'windows',W,'draws',draws, ...
           'age_healthy_max',age_healthy,'T_pkt',p.cyber.T_pkt);

A.age_max_observed = max(age);
A.dur_max_observed = 0;
if ~isempty(W), A.dur_max_observed = max(W(:,2) - W(:,1)); end
assert(A.dur_max_observed <= p.cyber.dos_dur_max + 1e-9, ...
    ['dos_generator: attack duration %.3f s exceeds p.cyber.dos_dur_max ' ...
     '= %.3f s. Raise dos_dur_max; the buffer horizon follows from it.'], ...
     A.dur_max_observed, p.cyber.dos_dur_max);
assert(A.age_max_observed <= p.cyber.tau_buffer, ...
    ['dos_generator: observed age of data %.3f s exceeds the delay buffer ' ...
     'horizon %.3f s (attack duration %.3f s + T_pkt %.3f + latency %.3f). ' ...
     'Raise p.cyber.dos_dur_max.'], ...
     A.age_max_observed, p.cyber.tau_buffer, A.dur_max_observed, ...
     p.cyber.T_pkt, p.cyber.tau_nom_hi);
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

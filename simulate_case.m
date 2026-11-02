function out = simulate_case(arm, p, scen, ctx)
%SIMULATE_CASE  Closed-loop simulation of the isolated microgrid over one scenario.

require_params(p);

spec = p.spec;  ix = spec.ix;
Ts   = spec.num.Ts_plant;
Tc   = spec.num.Ts_ctrl;
kC   = round(Tc/Ts);

scen = setdef(scen,'Tend',30);
scen = setdef(scen,'steps',[5 1.2*p.beta*abs(spec.safety.df_min_pu)]);
scen = setdef(scen,'t_ramp',0);
scen = setdef(scen,'res',zeros(0,2));
scen = setdef(scen,'Tout',zeros(0,2));
scen = setdef(scen,'dos',struct('mode','none'));
scen = setdef(scen,'x0',p.x0);

t = 0:Ts:scen.Tend;   N = numel(t);

dPL   = profile1(t, scen.steps, scen.t_ramp);
dPres = profile1(t, scen.res,   scen.t_ramp);
dTout = profile1(t, scen.Tout,  scen.t_ramp);

dosOpt = scen.dos;
mode   = dosOpt.mode;
if isfield(dosOpt,'mode'), dosOpt = rmfield(dosOpt,'mode'); end
A = dos_generator(p, t, mode, dosOpt);

x  = scen.x0(:);
if isfield(scen,'h0') && ~isempty(scen.h0)
    x(ix.df) = spec.safety.df_min_pu + scen.h0;
end

buf   = delay_buffer('init', p, x);
rbuf  = zeros(1, buf.N);   % matching RoCoF history
rhead = 1;

cs  = [];  u = [0;0];  ug = 0;
relay = [];  shed = 0;
rocof_f = 0;
aR = Ts/max(p.etm.T_rocof, Ts);

X   = zeros(spec.nx,N);   U = zeros(2,N);   UG = zeros(1,N);
SIG = zeros(1,N);  H = zeros(1,N);  PSI1 = zeros(1,N);  PSI2 = zeros(1,N);
STG = zeros(1,N);  SAFE = true(1,N);  PHI = zeros(1,N);  DEL = zeros(1,N);
SHED = zeros(1,N); ROC = zeros(1,N); TSOLVE = nan(1,N);
nPktSeen = 0;
cfg0 = struct('a_f',ctx.cfg.a_f,'Delta_cbf',0);

for n = 1:N
    d = [dPres(n); dPL(n)*(1 - shed); dTout(n)];

    F  = cpmg_dynamics(x, ug, d, p, 0);
    rocof_f = rocof_f + aR*(F(ix.df) - rocof_f);

    buf   = delay_buffer('push', buf, t(n), x);
    rhead = mod(rhead, buf.N) + 1;   rbuf(rhead) = rocof_f;

    [x_wan, ~, bi] = delay_buffer('fetch', buf, t(n), A.age(n));
    kback   = mod(rhead - bi.k - 1, buf.N) + 1;
    rocof_w = rbuf(kback);

    if mod(n-1, kC) == 0
        nWin = max(n-kC+1,1):n;
        net  = struct('alpha',A.alpha(n),'pkt',any(A.pkt(nWin)),'age',A.age(n));
        nPktSeen = nPktSeen + any(A.pkt(nWin));
        c   = ctx;  c.rocof_wan = rocof_w;
        [ug, u, cs, ci] = ctrl_bank(arm, t(n), x, x_wan, rocof_f, net, p, c, cs);
        SIG(n) = ci.sigma;
        if ~isempty(ci.qi)
            STG(n)    = ci.qi.stage;
            SAFE(n)   = ci.qi.safe_guaranteed;
            PSI2(n)   = ci.qi.psi2;
            DEL(n)    = ci.qi.delta;
            TSOLVE(n) = ci.qi.solve_time;
        else
            STG(n) = 0;  SAFE(n) = true;
        end
    else
        SIG(n) = SIG(max(n-1,1));  STG(n) = STG(max(n-1,1));
        SAFE(n) = SAFE(max(n-1,1));
    end

    B = barrier_functions(x, ug, d, p, 0, cfg0);
    H(n) = B.f.h;  PSI1(n) = B.f.psi1;  PHI(n) = B.dg.Phi;
    ROC(n) = F(ix.df);

    [shed, relay] = ufls_relay(x(ix.df), t(n), relay, spec);
    SHED(n) = shed;

    X(:,n) = x;  U(:,n) = u;  UG(n) = ug;

    if n < N
        x = cpmg_integrate(x, ug, u, d, p, 0, Ts, 'rk4');
    end
end

out = struct('arm',arm,'t',t,'X',X,'U',U,'ug',UG,'sigma',SIG, ...
    'alpha',A.alpha,'pkt',A.pkt,'age',A.age,'rocof',ROC, ...
    'h',H,'psi1',PSI1,'psi2',PSI2,'stage',STG,'safe_guaranteed',SAFE, ...
    'Phi',PHI,'delta',DEL,'shed',SHED,'solve_time',TSOLVE, ...
    'd',[dPres;dPL;dTout],'scen',scen,'dos',A,'relay',relay, ...
    'buf_stats',delay_buffer('stats',buf), ...
    'pkt_seen',nPktSeen,'pkt_sent',nnz(A.pkt), ...
    'pkt_seen_frac',nPktSeen/max(nnz(A.pkt),1));
end

function require_params(p)
NEED = 'params-2.6';
if ~isfield(p,'params_version') || ~strcmp(p.params_version, NEED)
    got = 'absent';
    if isfield(p,'params_version'), got = p.params_version; end
    error('simulate_case:staleParams', ...
      ['Parameter struct is stale: version "%s", need "%s".\n' ...
       'Rebuild it:  p = params_isolated();   (L and R are unaffected ' ...
       'unless the plant or spec changed.)'], got, NEED);
end
req = {'etm','cyber','uB_rate','uA_rate','SoC_enf_lo','aw_band'};
for k = 1:numel(req)
    if ~isfield(p, req{k})
        error('simulate_case:missingField', ...
          'Parameter struct is missing p.%s. Rebuild:  p = params_isolated();', req{k});
    end
end
reqEtm = {'tau_trigger','n_miss','T_react','T_dwell','n_ok','h_recover', ...
          'T_rocof','T_bump','k_recover'};
for k = 1:numel(reqEtm)
    if ~isfield(p.etm, reqEtm{k})
        error('simulate_case:missingField', ...
          'Parameter struct is missing p.etm.%s. Rebuild:  p = params_isolated();', ...
          reqEtm{k});
    end
end
end

function y = profile1(t, steps, tr)
y = zeros(size(t));
for k = 1:size(steps,1)
    t0 = steps(k,1);  a = steps(k,2);
    if tr <= 0
        y = y + a*(t >= t0);
    else
        r = min(max((t - t0)/tr, 0), 1);
        y = y + a*r;
    end
end
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

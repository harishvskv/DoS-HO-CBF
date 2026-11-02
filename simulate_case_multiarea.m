function out = simulate_case_multiarea(arm, P, scen, ctx)
%SIMULATE_CASE_MULTIAREA  Closed-loop simulation of the coupled N-area system over one scenario.

nA   = numel(P);
spec = P(1).spec;  ix = spec.ix;  nx = spec.nx;  nu = spec.nu;
Ts   = spec.num.Ts_plant;  Tc = spec.num.Ts_ctrl;  kC = round(Tc/Ts);
net_ = P(1).net;

for i = 1:nA, require_params39(P(i)); end

scen = setdef(scen,'Tend',60);
scen = setdef(scen,'steps',zeros(0,3));
scen = setdef(scen,'res',  zeros(0,3));
scen = setdef(scen,'Tout', zeros(0,3));
scen = setdef(scen,'t_ramp',0);
scen = setdef(scen,'dos',struct('mode','none'));

t = 0:Ts:scen.Tend;   N = numel(t);

DPL = zeros(nA,N);  DPR = zeros(nA,N);  DTO = zeros(nA,N);
for i = 1:nA
    DPL(i,:) = prof(t, scen.steps, i, scen.t_ramp);
    DPR(i,:) = prof(t, scen.res,   i, scen.t_ramp);
    DTO(i,:) = prof(t, scen.Tout,  i, scen.t_ramp);
end

dosOpt = scen.dos;  mode = dosOpt.mode;
if isfield(dosOpt,'mode'),  dosOpt = rmfield(dosOpt,'mode');  end
atkAreas = 1:nA;
if isfield(dosOpt,'areas'), atkAreas = dosOpt.areas; dosOpt = rmfield(dosOpt,'areas'); end

for i = nA:-1:1
    if any(atkAreas == i)
        A(i) = dos_generator(P(i), t, mode, dosOpt);
    else
        A(i) = dos_generator(P(i), t, 'none', dosOpt);
    end
end

X0 = zeros(nA*nx,1);
for i = 1:nA, X0((i-1)*nx + (1:nx)) = P(i).x0; end
if isfield(scen,'X0') && ~isempty(scen.X0), X0 = scen.X0(:); end
Xs = X0;

for i = 1:nA
    buf(i)  = delay_buffer('init', P(i), Xs((i-1)*nx + (1:nx)));   % #ok<AGROW>
    rbuf{i} = zeros(1, buf(i).N);   % #ok<AGROW>
    rhead(i) = 1;   % #ok<AGROW>
    cs{i} = [];   % #ok<AGROW>
    relay{i} = [];   % #ok<AGROW>
end
shed  = zeros(nA,1);
nSeen = zeros(nA,1);
rocof_f = zeros(nA,1);  tie_f = zeros(nA,1);
Ug = zeros(nA,1);  U = zeros(nA*nu,1);

XL   = zeros(nA*nx,N);  UL = zeros(nA*nu,N);  UGL = zeros(nA,N);
SIG  = zeros(nA,N);  HL = zeros(nA,N);  STG = zeros(nA,N);
SAFE = true(nA,N);   PHI = zeros(nA,N);  SHED = zeros(nA,N);
TIE  = zeros(nA,N);  ROC = zeros(nA,N);  PSI2 = zeros(nA,N);
aR = Ts/max(P(1).etm.T_rocof, Ts);
cfg0 = struct('a_f',ctx.cfg.a_f,'Delta_cbf',0);

for n = 1:N
    D = [DPR(:,n), DPL(:,n).*(1-shed), DTO(:,n)];

    [F, ~, DG] = cpmg_dynamics_multiarea(Xs, Ug, D, P);
    for i = 1:nA
        r = (i-1)*nx;
        rocof_f(i) = rocof_f(i) + aR*(F(r+ix.df)    - rocof_f(i));
        tie_f(i)   = tie_f(i)   + aR*(F(r+ix.dPtie) - tie_f(i));
        ROC(i,n) = F(r+ix.df);   TIE(i,n) = tie_f(i);

        buf(i)   = delay_buffer('push', buf(i), t(n), Xs(r+(1:nx)));
        rhead(i) = mod(rhead(i), buf(i).N) + 1;
        rbuf{i}(rhead(i)) = rocof_f(i);
    end

    if mod(n-1, kC) == 0
        nWin = max(n-kC+1,1):n;
        for i = 1:nA
            r = (i-1)*nx;  c = (i-1)*nu;
            [x_wan, ~, bi] = delay_buffer('fetch', buf(i), t(n), A(i).age(n));
            kb = mod(rhead(i) - bi.k - 1, buf(i).N) + 1;

            pktOK_i = any(A(i).pkt(nWin));
            nSeen(i) = nSeen(i) + pktOK_i;
            neti = struct('alpha',A(i).alpha(n), ...
                          'pkt',  pktOK_i, ...
                          'age',  A(i).age(n));
            ci_ctx = ctx;
            ci_ctx.tie_in    = tie_f(i);   % measured, survives isolation
            ci_ctx.rocof_wan = rbuf{i}(kb);
            ci_ctx.qp.tie_in = tie_f(i);

            [ug_i, u_i, cs{i}, cinfo] = ctrl_bank(arm, t(n), Xs(r+(1:nx)), ...
                                    x_wan, rocof_f(i), neti, P(i), ci_ctx, cs{i});
            Ug(i) = ug_i;  U(c+(1:nu)) = u_i;
            SIG(i,n) = cinfo.sigma;
            if ~isempty(cinfo.qi)
                STG(i,n)  = cinfo.qi.stage;
                SAFE(i,n) = cinfo.qi.safe_guaranteed;
                PSI2(i,n) = cinfo.qi.psi2;
            end
        end
    else
        SIG(:,n) = SIG(:,max(n-1,1));  STG(:,n) = STG(:,max(n-1,1));
        SAFE(:,n) = SAFE(:,max(n-1,1));
    end

    for i = 1:nA
        r = (i-1)*nx;
        B = barrier_functions(Xs(r+(1:nx)), Ug(i), D(i,:)', P(i), ...
                              DG(i).tie_in, cfg0);
        HL(i,n) = B.f.h;  PHI(i,n) = B.dg.Phi;
        [shed(i), relay{i}] = ufls_relay(Xs(r+ix.df), t(n), relay{i}, spec);
        SHED(i,n) = shed(i);
    end

    XL(:,n) = Xs;  UL(:,n) = U;  UGL(:,n) = Ug;

    if n < N
        Xs = cpmg_integrate_multiarea(Xs, Ug, U, D, P, Ts, 'rk4');
    end
end

out = struct('arm',arm,'t',t,'X',XL,'U',UL,'ug',UGL,'sigma',SIG, ...
    'h',HL,'stage',STG,'safe_guaranteed',SAFE,'Phi',PHI,'shed',SHED, ...
    'tie_rate',TIE,'rocof',ROC,'psi2',PSI2,'nA',nA,'scen',scen, ...
    'dos',A,'relay',{relay},'D',cat(3,DPR,DPL,DTO), ...
    'pkt_seen',nSeen,'pkt_sent',arrayfun(@(q) nnz(q.pkt), A(:)), ...
    'pkt_seen_frac',nSeen./max(arrayfun(@(q) nnz(q.pkt), A(:)),1));
end

function y = prof(t, steps, area, tr)
y = zeros(size(t));
if isempty(steps), return; end
for k = 1:size(steps,1)
    if steps(k,3) ~= area, continue; end
    t0 = steps(k,1);  a = steps(k,2);
    if tr <= 0
        y = y + a*(t >= t0);
    else
        y = y + a*min(max((t - t0)/tr, 0), 1);
    end
end
end

function require_params39(p)
NEED = 'params-2.6';
if ~isfield(p,'params_version') || ~strcmp(p.params_version, NEED)
    got = 'absent';
    if isfield(p,'params_version'), got = p.params_version; end
    error('simulate_case_multiarea:staleParams', ...
      ['Parameter struct is stale: version "%s", need "%s".\n' ...
       'Rebuild it:  P = params_39bus();'], got, NEED);
end
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

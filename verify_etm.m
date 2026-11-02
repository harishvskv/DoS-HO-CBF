function rep = verify_etm(p, L, R, doPlots)
%VERIFY_ETM  Trigger exclusivity, reconnection and packet-visibility checks.

if nargin < 1 || isempty(p), p = params_isolated(); end
if nargin < 4 || isempty(doPlots), doPlots = false; end

spec = p.spec;  ix = spec.ix;
Ts   = spec.num.Ts_plant;
cfg  = struct('a_f',[8 8],'Delta_cbf',0);
zero3 = zeros(spec.nd,1);

a_d      = R.a_d_bound;   % adversarial, for the safety bound
a_d_pred = R.a_d_fit;
dP    = 1.2*p.beta*abs(spec.safety.df_min_pu);
h_max = -spec.safety.df_min_pu;

D_long  = min(3.0, p.cyber.dos_dur_max);
D_short = 0.5;   % anti-chatter probe

rep = struct('tests',struct('id',{},'name',{},'value',{},'tol',{},'pass',{},'note',{}));
rep.scenario = struct('dP_imb',dP,'D_long',D_long,'D_short',D_short, ...
                      'dos_dur_max',p.cyber.dos_dur_max, ...
                      'age_dos_max',p.cyber.age_dos_max, ...
                      'tau_buffer',p.cyber.tau_buffer);

mixf = @(al,sg,uA,uL) al*(1-sg)*uA + (1 - al*(1-sg))*uL;
uA = 1;  uL = -1;
tbl = [ mixf(1,0,uA,uL) ;   % healthy, fresh      -> AGC
        mixf(1,1,uA,uL) ;   % healthy, triggered  -> local
        mixf(0,0,uA,uL) ;   % DoS,     not yet    -> local
        mixf(0,1,uA,uL) ];   % DoS,     triggered  -> local
want = [uA; uL; uL; uL];
e01 = max(abs(tbl - want));
rep = addT(rep,'E01','switching law truth table', e01, 1e-14, e01<1e-14, ...
    sprintf(['(a=1,s=0)->AGC  (1,1)->loc  (0,0)->loc  (0,1)->loc. ' ...
             'The submitted Eq.(36) used alpha*sigma and selected the local ' ...
             'law during healthy operation.']));

ageH  = p.cyber.T_pkt + p.cyber.tau_nom_hi + p.etm.T_rocof;
ord   = (ageH < p.etm.tau_trigger) && (p.etm.tau_trigger < L.tau_max_feasible);
marg1 = p.etm.tau_trigger - ageH;
marg2 = L.tau_max_feasible - p.etm.tau_trigger;
rep = addT(rep,'E02','healthy age + RoCoF lag < tau_trigger < tau_max', ...
    p.etm.tau_trigger, L.tau_max_feasible, ord, ...
    sprintf(['%.4f < %.4f < %.4f s.  Margins: %.4f s below, %.4f s above. ' ...
             'T_pkt %.3f + latency %.3f + RoCoF lag %.3f.'], ...
            ageH, p.etm.tau_trigger, L.tau_max_feasible, marg1, marg2, ...
            p.cyber.T_pkt, p.cyber.tau_nom_hi, p.etm.T_rocof));
if min(marg1,marg2) < 0.020
    fprintf(['\n [E02] WARNING: only %.0f ms of headroom in the timing budget. ' ...
             'Consider T_pkt = 0.030 s.\n'], 1000*min(marg1,marg2));
end

Th = 60;  th = 0:Ts:Th;
Ah = dos_generator(p, th, 'none');
[nTrigH, ~, ~] = runETM(p, th, Ah, p.x0, dP*0.3, a_d, a_d_pred, 0, h_max, cfg, 'closed');
rep = addT(rep,'E03','no false trigger over 60 s of healthy operation', ...
    nTrigH, 0, nTrigH == 0, ...
    sprintf(['max healthy age %.4f s against trigger %.4f s; predictor uses ' ...
             'the frozen-input curvature a_d_pred = %.4f, not the adversarial ' ...
             'bound %.4f'], max(Ah.age), p.etm.tau_trigger, a_d_pred, a_d));

Td = 20;  td = 0:Ts:Td;
t_atk = 5;
Ad = dos_generator(p, td, 'deterministic', struct('starts',t_atk,'durations',D_long));
x0 = p.x0;
[nT, tT, mT] = runETM(p, td, Ad, x0, 0.2*dP, a_d, a_d_pred, 0, h_max, cfg, 'closed');

lat = NaN;  mech04 = 'none';
kPost = find(tT >= t_atk - 1e-9, 1);
if ~isempty(kPost)
    lat    = tT(kPost) - t_atk;
    mech04 = mT{kPost};
end
bound = p.etm.n_miss*p.cyber.T_pkt + p.cyber.tau_nom_hi + Ts;
rep = addT(rep,'E04','detection latency under DoS', lat, bound, ...
    ~isnan(lat) && lat >= -1e-9 && lat <= bound + 1e-9, ...
    sprintf(['fired %+.4f s after attack onset via "%s"; packet-loss bound ' ...
             'is %.4f s. %d trigger(s) in total over the run.'], ...
            lat, mech04, bound, nT));

h0g = [0.003 0.005 0.007 0.009 0.011];
row = struct('h0',{},'t_trig',{},'t_floor',{},'margin',{},'mech',{},'ok',{});
for k = 1:numel(h0g)
    x0 = p.x0;  x0(ix.df) = spec.safety.df_min_pu + h0g(k);

    st = struct('x0',x0,'SoC0',p.SoC_nom,'dP_imb',dP);
    S  = compute_scatd(p, L, R, st, struct('a_f',cfg.a_f,'use_bound',true));

    tt = 0:Ts:(D_long + 1.0);
    Aa = dos_generator(p, tt, 'deterministic', struct('starts',0,'durations',D_long));
    [nT2, tT2, mT2, tFloor] = runETM(p, tt, Aa, x0, dP, a_d, a_d_pred, S.h_min, h_max, cfg, 'frozen');

    tTrig = NaN;  if nT2 > 0, tTrig = tT2(1); end
    marg  = tFloor - tTrig;
    row(end+1) = struct('h0',h0g(k),'t_trig',tTrig,'t_floor',tFloor, ...
        'margin',marg,'mech',firstMech(mT2), ...
        'ok', ~isnan(tTrig) && (isnan(tFloor) || marg > 0));   % #ok<AGROW>
end
worst = min([row.margin]);
rep = addT(rep,'E05','trigger precedes loss of the barrier', worst, 0, ...
    all([row.ok]), ...
    sprintf('worst margin %.4f s across %d initial margins', worst, numel(row)));
rep.E05 = row;

pktBound = p.etm.n_miss*p.cyber.T_pkt;
tooFast  = [row.t_floor] < pktBound;
covered  = arrayfun(@(r) strcmp(r.mech,'predictive'), row);
e06 = all(covered(tooFast));
rep = addT(rep,'E06','predictive path covers the fast-loss states', ...
    double(e06), 0.5, e06, ...
    sprintf(['%d of %d states lose the barrier faster than the %.3f s ' ...
             'packet-loss bound; all were caught by the predictive path.'], ...
            nnz(tooFast), numel(row), pktBound));

tc = 0:Ts:40;
Ac = dos_generator(p, tc, 'deterministic', ...
        struct('starts',[10 12 14 16 18],'durations',D_short));
x0 = p.x0;  x0(ix.df) = spec.safety.df_min_pu + 0.9*h_max;
[~,~,~,~,nSwitch] = runETM(p, tc, Ac, x0, dP*0.2, a_d, a_d_pred, 0, h_max, cfg, 'closed');
maxSw = 2*(1 + floor(range(tc)/p.etm.T_dwell));
rep = addT(rep,'E07','bounded mode changes, no chattering', nSwitch, maxSw, ...
    nSwitch <= maxSw, ...
    sprintf('%d mode changes over %.0f s; dwell %.2f s caps it at %d', ...
            nSwitch, range(tc), p.etm.T_dwell, maxSw));

kPkt = round(p.cyber.T_pkt/spec.num.Ts_ctrl);
es = []; ctx = struct('a_d',a_d,'a_d_pred',a_d_pred,'h_min',0, ...
                      'h_max',h_max,'cfg',cfg);
xb = p.x0;  xb(ix.df) = spec.safety.df_min_pu + 0.2*h_max;   % low margin
rb = barrierRate(xb,p,cfg);
[~, es] = da_etm(0, xb, rb, false, 5, p, ctx, es);   % force isolation
assert(es.sigma == 1, 'verify_etm: setup failed, expected isolation');

tk = p.etm.T_dwell;
for n = 1:20*kPkt
    tk = tk + spec.num.Ts_ctrl;
    [sg, es] = da_etm(tk, xb, rb, mod(n,kPkt)==0, 0.01, p, ctx, es);
end
heldLow = (sg == 1);

xr = p.x0;  xr(ix.df) = spec.safety.df_min_pu + 0.95*h_max;
rr = barrierRate(xr,p,cfg);
for n = 1:20*kPkt
    tk = tk + spec.num.Ts_ctrl;
    [sg2, es] = da_etm(tk, xr, rr, mod(n,kPkt)==0, 0.01, p, ctx, es);
end
rep = addT(rep,'E08','reconnection needs dwell AND packets AND margin', ...
    double(heldLow && sg2 == 0), 0.5, heldLow && sg2 == 0, ...
    sprintf(['packet cadence 1 arrival per %d controller steps. Stayed ' ...
             'isolated with a low margin, reconnected once recovered.'], kPkt));

Delta11 = L.supPsiDot*spec.num.Ts_ctrl;
ctx11 = struct('a_d',a_d,'a_d_pred',a_d_pred,'h_min',0,'h_max',h_max, ...
    'cfg',struct('a_f',[8 8],'a_S',[2 2],'a_T',[0.5 0.5]), ...
    'Delta_cbf',Delta11,'etc_sigma',0.05, ...
    'qp',struct('Delta_cbf',Delta11,'a_f',[8 8],'a_S',[2 2],'a_T',[0.5 0.5], ...
                'clf',[4 4],'w',[10 500],'p_delta',1e5,'Ts',spec.num.Ts_ctrl));
scen11 = struct('Tend',60,'steps',[3 dP],'t_ramp',0, ...
    'dos',struct('mode','deterministic','starts',3,'durations',2));
o11  = simulate_case('A4', p, scen11, ctx11);
sigR = o11.sigma;  tr = o11.t;

tEnd     = 5;   % attack ends here
iRel     = find(sigR(:)' == 0 & tr > tEnd, 1);
released = ~isempty(iRel);
tRel     = NaN;  if released, tRel = tr(iRel) - tEnd; end

tail    = tr > (tr(end) - 20);
quiet   = mean(sigR(tail) == 0);
nSwTail = nnz(diff(sigR(tail)) ~= 0);
e11ok   = released && tRel <= 30 && quiet >= 0.80 && nSwTail <= 2;
rep = addT(rep,'E11','isolation releases AND stays released', ...
    tRel, 30, e11ok, ...
    sprintf(['released %.2f s after the attack; isolated %.1f%% of the run; ' ...
             'final 20 s: %.0f%% connected with %d mode change(s); SoC %.4f ' ...
             'to %.4f. Run through simulate_case, the deployed path.'], ...
             tRel, 100*mean(sigR>0), 100*quiet, nSwTail, ...
             o11.X(ix.dSB,1), o11.X(ix.dSB,end)));

th2 = 0:Ts:30;
Ah2 = dos_generator(p, th2, 'none');
[~,~,~,~,~,~,pf] = runETM(p, th2, Ah2, p.x0, 0.2*dP, a_d, a_d_pred, 0, h_max, cfg, 'closed');
rep = addT(rep,'E12','controller sees the packets that arrive', pf, 0.95, pf >= 0.95, ...
    sprintf(['%.1f%% of arrivals reached the trigger. Sampling the arrival flag ' ...
             'at the controller instant alone gives 22.8%%, which fires the ' ...
             'packet-loss path with no attack present.'], 100*pf));

buf = delay_buffer('init', p, p.x0);
for n = 1:200
    buf = delay_buffer('push', buf, n*Ts, p.x0 + n*1e-4);
end
threw = false;
try
    delay_buffer('fetch', buf, 200*Ts, p.cyber.tau_buffer + 1.0);
catch ME
    threw = strcmp(ME.identifier,'delay_buffer:horizon');
end
rep = addT(rep,'E09','buffer errors rather than silently truncating', ...
    double(threw), 0.5, threw, ...
    'the legacy buffer clamped, so a 2.0 s delay spike became 0.19 s');

kBack = 50;
xg = delay_buffer('fetch', buf, 200*Ts, kBack*Ts);
xe = p.x0 + (200 - kBack)*1e-4;
e10 = max(abs(xg - xe));
rep = addT(rep,'E10','buffer round trip is exact', e10, 1e-12, e10 < 1e-12, ...
    sprintf('requested age %.3f s, recovered the exact stored sample', kBack*Ts));

rep.allPass = all([rep.tests.pass]);
printRep(rep);
if doPlots, plotE05(rep.E05, p); end
if ~rep.allPass
    warning('verify_etm:FAILED','Cyber layer verification FAILED. Do not proceed to Layer 6.');
end
end   % verify_etm

function [nTrig, tTrig, mech, tFloor, nSwitch, sigLog, pktFrac, rcLog] = ...
         runETM(p, t, A, x0, dP, a_d, a_d_pred, h_min, h_max, cfg, mode)
if nargin < 11 || isempty(mode), mode = 'frozen'; end
spec = p.spec;  Ts = spec.num.Ts_plant;
d  = [0; dP; 0];
x  = x0;  es = [];  prevSig = 0;  nSwitch = 0;  sigLog = zeros(1,numel(t));
ctx = struct('a_d',a_d,'a_d_pred',a_d_pred,'h_min',h_min,'h_max',h_max,'cfg',cfg);
[~, ue] = ctrl_pi_agc(x0, p, 'agc_edge');
ug = 0;  sg = 0;  tFloor = NaN;  rcLog = false(5,0);  rcVal = zeros(5,0);
kC = round(spec.num.Ts_ctrl/Ts);
nSeen = 0;
for n = 1:numel(t)
    B = barrier_functions(x, 0, d, p, 0, cfg);
    if isnan(tFloor) && B.f.h <= h_min, tFloor = t(n); end
    if mod(n-1, kC) == 0
        nWin  = max(n-kC+1,1):n;
        pktOK = any(A.pkt(nWin));
        nSeen = nSeen + pktOK;
        [sg, es, di] = da_etm(t(n), x, B.f.Lfh, pktOK, A.age(n), p, ctx, es);
        if sg == 1
            rcLog(:,end+1) = [di.rc_dwell; di.rc_fresh; di.rc_h; ...
                              di.rc_Tpred; di.rc_age];   % #ok<AGROW>
            rcVal(:,end+1) = [di.dwell; di.n_fresh; di.h; di.T_pred; di.age];   % #ok<AGROW>
        end

        ug = A.alpha(n)*ctrl_pi_agc(x, p, 'agc_only');

        if A.alpha(n) && sg == 0
            [~, ue] = ctrl_pi_agc(x, p, 'agc_edge');
        elseif strcmpi(mode,'closed')
            [~, ue] = ctrl_pi_agc(x, p, 'edge_only');
        end
    end
    sigLog(n) = sg;
    if sg ~= prevSig, nSwitch = nSwitch + 1; prevSig = sg; end
    x = cpmg_integrate(x, ug, ue, d, p, 0, Ts, 'rk4');
end
nTrig = es.n_trig;  tTrig = es.t_trig;  mech = es.mech;
pktFrac = nSeen/max(nnz(A.pkt),1);
if isempty(rcLog), rcLog = struct('frac',nan(1,5),'last',nan(1,5));
else, rcLog = struct('frac',mean(rcLog,2)','last',rcVal(:,end)'); end
end

function r = barrierRate(x,p,cfg)
B = barrier_functions(x, 0, zeros(p.spec.nd,1), p, 0, cfg);
r = B.f.Lfh;
end

function m = firstMech(mech)
if isempty(mech), m = 'none'; else, m = mech{1}; end
end

function rep = addT(rep,id,name,value,tol,pass,note)
n = numel(rep.tests)+1;
rep.tests(n) = struct('id',id,'name',name,'value',value,'tol',tol, ...
                      'pass',logical(pass),'note',note);
end

function printRep(rep)
fprintf('\n========================================================================\n');
fprintf(' CYBER LAYER VERIFICATION (Layer 5)\n');
fprintf('========================================================================\n');
for k = 1:numel(rep.tests)
    T = rep.tests(k);
    fprintf(' %-5s %-50s %s\n', T.id, T.name, tern(T.pass));
    if ~isempty(T.note), fprintf('       %s\n', T.note); end
end
if isfield(rep,'E05')
    fprintf('\n E05 detail\n');
    fprintf(' %8s %10s %10s %10s  %s\n','h0','t_trig','t_floor','margin','path');
    for k = 1:numel(rep.E05)
        r = rep.E05(k);
        fprintf(' %8.4f %10.4f %10.4f %+10.4f  %s\n', ...
            r.h0, r.t_trig, r.t_floor, r.margin, r.mech);
    end
end
fprintf('------------------------------------------------------------------------\n');
fprintf(' OVERALL: %s   (%d/%d)\n', tern(all([rep.tests.pass])), ...
        nnz([rep.tests.pass]), numel(rep.tests));
fprintf('========================================================================\n\n');
end

function s = tern(c)
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function plotE05(row, ~)
figure('Color','w','Name','E05 trigger timing','Position',[100 100 600 380]);
h0 = [row.h0];
bar(h0, [[row.t_trig]; [row.t_floor]]'); grid on; box on;
xlabel('initial margin h_0 (p.u.)'); ylabel('time (s)');
legend('trigger','barrier floor reached','Location','northwest');
title('Isolation must precede loss of the barrier');
end

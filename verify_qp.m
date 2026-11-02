function rep = verify_qp(p, L, R, doPlots)
%VERIFY_QP  Feasibility-cascade, slack and recovery checks on the quadratic program.

if nargin < 1 || isempty(p), p = params_isolated(); end
if nargin < 4 || isempty(doPlots), doPlots = false; end

spec = p.spec;  ix = spec.ix;
Ts   = spec.num.Ts_plant;
Tc   = spec.num.Ts_ctrl;
h_max = -spec.safety.df_min_pu;
dP    = 1.2*p.beta*abs(spec.safety.df_min_pu);

Delta = L.supPsiDot*Tc;   % local-loop tightening, Theorem 2
qopt  = struct('Delta_cbf',Delta,'a_f',[8 8],'a_S',[2 2],'a_T',[0.5 0.5], ...
               'clf',[4 4],'w',[10 500],'p_delta',1e5,'Ts',Tc);

rep = struct('tests',struct('id',{},'name',{},'value',{},'tol',{},'pass',{},'note',{}));

zero3 = zeros(spec.nd,1);
nomRate = @(xx) nomRateD(xx, zero3, p, ix);

x0 = p.x0;
[~, ~, qi0] = ctrl_hocbf_qp(x0, nomRate(x0), [0;0], p, [], qopt);
names = qi0.rows;
need  = {'CBF:freq','CBF:SoC_lo','CBF:SoC_hi','CBF:Tin_lo','CBF:Tin_hi'};
have  = all(ismember(need, names));
rep = addT(rep,'Q01','all five barriers are rows on the decision variable', ...
    double(have), 0.5, have, ...
    sprintf(['rows: %s. The submitted QP had only the frequency row; SoC ' ...
             'appeared nowhere and comfort was written as a bound on a STATE.'], ...
            strjoin(names,', ')));

rep = addT(rep,'Q02','CLF slack is zero when the problem is feasible', ...
    qi0.delta, 1e-6, qi0.delta < 1e-6 && qi0.stage == 1, ...
    sprintf('stage %d, delta = %.3e', qi0.stage, qi0.delta));

xf = p.x0;  xf(ix.df) = spec.safety.df_min_pu + 0.98*h_max;
uref = [0.05; -0.02];
csf = struct('u_prev',uref,'n_stage',zeros(1,4),'n_solve',0);
[uf, ~, qif] = ctrl_hocbf_qp(xf, nomRate(xf), uref, p, csf, qopt);
e03 = norm(uf - uref);
rep = addT(rep,'Q03','filter is transparent far from the barriers', ...
    e03, 1e-6, e03 < 1e-6, ...
    sprintf('||u - u_ref|| = %.3e at %.0f%% of the safe set, stage %d', ...
            e03, 98, qif.stage));

[h_crit, hinfo] = compute_hcrit(p, dP, struct('qp',qopt));
rep.h_crit = h_crit;  rep.h_crit_info = hinfo;

stF = struct('x0',p.x0,'SoC0',p.SoC_nom,'dP_imb',dP);
SF  = compute_scatd(p, L, R, stF, struct('a_f',qopt.a_f,'use_bound',true));
rep.h_min_of_SoC = SF.h_min_of_SoC;

h0g = h_crit + (h_max - h_crit)*[0.05 0.20 0.40 0.65 0.90];
row = struct('h0',{},'h_min_traj',{},'floor',{},'st4_frac',{},'minPsi2',{},'ok',{});
Tend = 6.0;
for k = 1:numel(h0g)
    x  = p.x0;  x(ix.df) = spec.safety.df_min_pu + h0g(k);
    st = struct('x0',x,'SoC0',p.SoC_nom,'dP_imb',dP);
    S  = compute_scatd(p, L, R, st, struct('a_f',qopt.a_f,'use_bound',true));

    d  = [0; dP; 0];
    cs = [];  u = [0;0];
    nStep = round(Tend/Ts);  kCtrl = round(Tc/Ts);
    hmin_seen = inf;  psi2min = inf;  n4 = 0;  nC = 0;
    for n = 1:nStep
        if mod(n-1,kCtrl) == 0
            [~, uref] = ctrl_pi_agc(x, p, 'edge_only');
            rc = nomRateD(x,d,p,ix);
            [u, cs, qi] = ctrl_hocbf_qp(x, rc, uref, p, cs, qopt);
            nC = nC + 1;
            if qi.stage == 4, n4 = n4 + 1; end
            if qi.stage <= 3, psi2min = min(psi2min, qi.psi2_margin); end
        end
        Bn = barrier_functions(x, 0, d, p, 0, struct('a_f',qopt.a_f,'Delta_cbf',0));
        hmin_seen = min(hmin_seen, Bn.f.h);
        x = cpmg_integrate(x, 0, u, d, p, 0, Ts, 'rk4');
    end
    row(end+1) = struct('h0',h0g(k),'h_min_traj',hmin_seen,'floor',S.h_min, ...
        'st4_frac',n4/max(nC,1),'minPsi2',psi2min, ...
        'ok', hmin_seen >= S.h_min - 1e-6);   % #ok<AGROW>
end
rep.Q05 = row;

okPsi = min([row.minPsi2]) >= -1e-6;
rep = addT(rep,'Q04','psi2 >= Delta_cbf at the QP solution', ...
    min([row.minPsi2]), 0, okPsi, ...
    'evaluated at every controller step in which the QP reported stage 1 to 3');

okInv  = all([row.ok]);
worstH = min([row.h_min_traj]);
nBreach = nnz(~[row.ok]);
rep = addT(rep,'Q05','closed-loop invariance for h0 >= h_crit', ...
    worstH, 0, okInv, ...
    sprintf(['h_crit = %.5f p.u. (%.1f%% of the safe set); lowest barrier ' ...
             'reached %.5f over %d runs above it, %d breach(es); stage-4 ' ...
             'usage %.1f%% worst case. A pass certifies the TESTED states ' ...
             'under the TESTED disturbance, not a universal guarantee.'], ...
            h_crit, 100*hinfo.frac, worstH, numel(row), nBreach, ...
            100*max([row.st4_frac])));

sub = h_crit*[0.30 0.50 0.70];
rowS = struct('h0',{},'h_min_traj',{},'floor',{},'st4_frac',{}, ...
              'minPsi2',{},'breached',{},'flagged',{},'ok',{});
d2 = [0; dP; 0];  kCtrl = round(Tc/Ts);
for j = 1:numel(sub)
    h_below = max(sub(j), 1e-4);
    xb2 = p.x0;  xb2(ix.df) = spec.safety.df_min_pu + h_below;
    cs2 = [];  u2 = [0;0];  flagged = false;  breached = false;
    hmin_s = inf;  n4 = 0;  nC2 = 0;
    nStep = round(3.0/Ts);
    for n = 1:nStep
        if mod(n-1,kCtrl) == 0
            [~, ur2] = ctrl_pi_agc(xb2, p, 'edge_only');
            [u2, cs2, q2] = ctrl_hocbf_qp(xb2, nomRateD(xb2,d2,p,ix), ur2, p, cs2, qopt);
            nC2 = nC2 + 1;
            if q2.stage == 4, n4 = n4 + 1; end
            if ~q2.safe_guaranteed || q2.stage == 4, flagged = true; end
        end
        Bb = barrier_functions(xb2, 0, d2, p, 0, struct('a_f',qopt.a_f,'Delta_cbf',0));
        hmin_s = min(hmin_s, Bb.f.h);
        if Bb.f.h < 0, breached = true; end
        xb2 = cpmg_integrate(xb2, 0, u2, d2, p, 0, Ts, 'rk4');
    end
    rowS(j) = struct('h0',h_below,'h_min_traj',hmin_s,'floor',0, ...
        'st4_frac',n4/max(nC2,1),'minPsi2',NaN, ...
        'breached',breached,'flagged',flagged, ...
        'ok', ~breached || flagged);   % #ok<AGROW>
end
rep.Q05sub = rowS;
q12ok = all([rowS.ok]);
nB = nnz([rowS.breached]);
rep = addT(rep,'Q12','a breach below h_crit is declared, never silent', ...
    double(q12ok), 0.5, q12ok, ...
    sprintf(['%d sub-critical margins tested (%.0f%%, %.0f%%, %.0f%% of ' ...
             'h_crit = %.5f): %d breached, every breach declared. The ' ...
             'guarantee is conditional and the controller says so.'], ...
            numel(rowS), 100*sub/h_crit, h_crit, nB));

xt = p.x0;  xt(ix.df) = spec.safety.df_min_pu + 0.4*h_max;
taus = [0 0.01 0.05 0.10];
auth = zeros(size(taus));  uBk = zeros(size(taus));
for k = 1:numel(taus)
    q = qopt;  q.Delta_cbf = L.supPsiDot*taus(k);  q.use_rate = false;
    [ut,~,qt] = ctrl_hocbf_qp(xt, nomRateD(xt,[0;dP;0],p,ix), [0;0], p, [], q);
    auth(k) = qt.LgLfh*ut;
    uBk(k)  = ut(1);
end
mono = all(diff(auth) >= -1e-9);
rep = addT(rep,'Q06','larger delay demands more barrier authority', ...
    min(diff(auth)), 0, mono, ...
    sprintf('LgLf h * u across tau = [%s] s: [%s]', ...
            num2str(taus,'%.2f '), num2str(auth,'%.4f ')));

xd = p.x0;
xd(ix.df)  = spec.safety.df_min_pu + 0.10*h_max;   % thin margin
xd(ix.dSB) = p.SoC_min;   % battery exhausted
[~, ~, qd] = ctrl_hocbf_qp(xd, nomRateD(xd,[0;dP;0],p,ix), [0;0], p, [], qopt);
esc = qd.stage >= 2 || ~qd.safe_guaranteed;
rep = addT(rep,'Q07','depleted battery escalates the ladder and is flagged', ...
    qd.stage, 1, esc, ...
    sprintf(['stage %d, Phi = %.4f, safe_guaranteed = %d. Remark 4 of the ' ...
             'submission claimed feasibility is "seamlessly restored"; it is ' ...
             'not, and the ladder records which stage was used.'], ...
            qd.stage, qd.Phi, qd.safe_guaranteed));

xa = p.x0;  xa(ix.df) = spec.safety.df_min_pu + 0.02*h_max;
dBig = 5*dP;
[~, ~, qa] = ctrl_hocbf_qp(xa, nomRateD(xa,[0;dBig;0],p,ix), [0;0], p, [], qopt);
rep = addT(rep,'Q08','imbalance beyond authority is flagged, not smoothed', ...
    double(~qa.safe_guaranteed || qa.stage == 4), 0.5, ...
    (~qa.safe_guaranteed || qa.stage == 4), ...
    sprintf('imbalance %.2f p.u., stage %d, safe_guaranteed = %d', ...
            dBig, qa.stage, qa.safe_guaranteed));

csr = struct('u_prev',[0;0],'n_stage',zeros(1,4),'n_solve',0);
xr  = p.x0;  xr(ix.df) = spec.safety.df_min_pu + 0.05*h_max;
[ur, ~, ~] = ctrl_hocbf_qp(xr, nomRateD(xr,[0;dP;0],p,ix), [1;-1], p, csr, qopt);
dmax = [p.uB_rate; p.uA_rate]*Tc;
viol = max(abs(ur) - dmax);
rep = addT(rep,'Q09','actuator rate rows are respected', viol, 1e-9, viol <= 1e-9, ...
    sprintf('|du| = [%.4f %.4f] against limits [%.4f %.4f]', abs(ur), dmax));

t_reach_B = dP/p.uB_rate;
t_reach_A = dP/p.uA_rate;
hd0   = -dP/p.M;
t_avail = h_max/abs(hd0);   % time from a full margin to zero
reachOK = t_reach_B < 0.5*t_avail;
rep = addT(rep,'Q11','BESS can slew to cover the imbalance in time', ...
    t_reach_B, 0.5*t_avail, reachOK, ...
    sprintf(['imbalance %.3f p.u.: BESS needs %.3f s to slew, AIAC needs ' ...
             '%.1f s, and a FULL margin survives only %.3f s. The AIAC is ' ...
             'therefore irrelevant on this horizon and fast frequency safety ' ...
             'rests entirely on the BESS.'], dP, t_reach_B, t_reach_A, t_avail));

Delta13 = L.supPsiDot*Tc;
ctx13 = struct('a_d',R.a_d_bound,'a_d_pred',R.a_d_fit,'h_min',0,'h_max',h_max, ...
    'cfg',struct('a_f',qopt.a_f,'a_S',qopt.a_S,'a_T',qopt.a_T), ...
    'Delta_cbf',Delta13,'etc_sigma',0.05, ...
    'qp',struct('Delta_cbf',Delta13,'a_f',qopt.a_f,'a_S',qopt.a_S, ...
                'a_T',qopt.a_T,'clf',[4 4],'w',[10 500],'p_delta',1e5,'Ts',Tc));
scen13 = struct('Tend',120,'steps',[3 dP],'t_ramp',0, ...
    'dos',struct('mode','deterministic','starts',3,'durations',2));
o13 = simulate_case('A4', p, scen13, ctx13);

tail13 = o13.t > 0.8*o13.t(end);
PBtail = mean(abs(o13.X(ix.dPB,tail13)));
SoC13  = o13.X(ix.dSB,:);
drift  = (SoC13(1) - SoC13(end))/o13.t(end);
q13ok  = PBtail < 0.02 && abs(drift) < 2e-4;
rep = addT(rep,'Q13','storage returns to rest after the transient', ...
    PBtail, 0.02, q13ok, ...
    sprintf(['mean |dPB| over the last 20%% of a %.0f s run is %.4f p.u.; SoC ' ...
             '%.4f to %.4f, drift %.2e per s. Ki_e = %.2f, K_soc = %.2f, ' ...
             'restoration time constant %.0f s.'], o13.t(end), PBtail, ...
             SoC13(1), SoC13(end), drift, p.agc.Ki_e, p.agc.K_soc, ...
             p.agc.T_soc_restore));

socg = [0.30 0.50 0.70];
hminS = zeros(size(socg));
for j = 1:numel(socg)
    xq2 = p.x0;  xq2(ix.dSB) = socg(j);
    dq2 = [0; dP; 0];  csq2 = [];  uq2 = [0;0];
    nStep = round(6.0/Ts);  kCtrl = round(Tc/Ts);
    hm = inf;
    for n = 1:nStep
        if mod(n-1,kCtrl) == 0
            [~, ur] = ctrl_pi_agc(xq2, p, 'edge_only');
            [uq2, csq2] = ctrl_hocbf_qp(xq2, nomRateD(xq2,dq2,p,ix), ur, p, csq2, qopt);
        end
        Bq = barrier_functions(xq2, 0, dq2, p, 0, struct('a_f',qopt.a_f,'Delta_cbf',0));
        hm = min(hm, Bq.f.h);
        xq2 = cpmg_integrate(xq2, 0, uq2, dq2, p, 0, Ts, 'rk4');
    end
    hminS(j) = hm;
end
spreadS = max(abs(hminS - hminS(2)));
tolS    = 0.02*(-spec.safety.df_min_pu);   % two per cent of the safe set
q14ok   = spreadS <= tolS;
rep = addT(rep,'Q14','barrier is independent of initial state of charge', ...
    spreadS, tolS, q14ok, ...
    sprintf(['min h at SoC0 = 0.30/0.50/0.70 is %.5f / %.5f / %.5f; spread ' ...
             '%.5f against a tolerance of %.5f. Phi is unity throughout this ' ...
             'range, so any slope means charge is coupling into the reference.'], ...
            hminS, spreadS, tolS));

[ua, ~, ~] = ctrl_hocbf_qp(xr, nomRateD(xr,[0;dP;0],p,ix), [0.1;0], p, [], qopt);
[ub2,~, ~] = ctrl_hocbf_qp(xr, nomRateD(xr,[0;dP;0],p,ix), [0.1;0], p, [], qopt);
e10 = norm(ua - ub2);
rep = addT(rep,'Q10','identical inputs give identical outputs', e10, 1e-12, ...
    e10 < 1e-12, '');

rep.allPass = all([rep.tests.pass]);
printRep(rep);
if doPlots, plotQ05(rep, p, h_max); end
if ~rep.allPass
    warning('verify_qp:FAILED','QP verification FAILED. Do not proceed to Layer 7.');
end
end   % verify_qp

function r = nomRateD(x, d, p, ix)
F = cpmg_dynamics(x, 0, d, p, 0);
r = F(ix.df);
end

function rep = addT(rep,id,name,value,tol,pass,note)
n = numel(rep.tests)+1;
rep.tests(n) = struct('id',id,'name',name,'value',value,'tol',tol, ...
                      'pass',logical(pass),'note',note);
end

function printRep(rep)
fprintf('\n========================================================================\n');
fprintf(' QP SAFETY FILTER VERIFICATION (Layer 6)\n');
fprintf('========================================================================\n');
for k = 1:numel(rep.tests)
    T = rep.tests(k);
    fprintf(' %-5s %-50s %s\n', T.id, T.name, tern(T.pass));
    if ~isempty(T.note), fprintf('       %s\n', T.note); end
end
if isfield(rep,'Q05')
    fprintf('\n Q05 detail\n');
    fprintf(' %8s %12s %10s %10s %10s\n','h0','min h(t)','floor','stage4','min psi2');
    for k = 1:numel(rep.Q05)
        r = rep.Q05(k);
        fprintf(' %8.4f %12.5f %10.5f %9.1f%% %10.4f\n', ...
            r.h0, r.h_min_traj, r.floor, 100*r.st4_frac, r.minPsi2);
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

function plotQ05(rep, p, h_max)
row = rep.Q05;
K   = 1e3;
figure('Color','w','Name','Q05 closed-loop invariance','Position',[80 80 1020 400]);

ax1 = subplot(1,2,1);
h0 = [row.h0]*K;  hm = [row.h_min_traj]*K;  fl = [row.floor]*K;
b = bar(h0, hm, 0.55); b.FaceColor = [0.20 0.45 0.75];
hold on; grid on; box on;

ylim([0, max(hm)*1.35]);
xlim([min(h0)-1, max(h0)+1]);

plot(xlim, [mean(fl) mean(fl)], 'r--', 'LineWidth', 1.6);

plot([rep.h_crit*K rep.h_crit*K], ylim, 'm-.', 'LineWidth', 1.6);
text(rep.h_crit*K, max(ylim)*0.97, sprintf(' h_{crit} = %.2f', rep.h_crit*K), ...
     'Color','m','FontSize',9,'VerticalAlignment','top');

xlabel('initial margin h_0  (10^{-3} p.u.)');
ylabel('lowest barrier reached  (10^{-3} p.u.)');
title('Barrier held above the floor for h_0 \geq h_{crit}');

legend([b], {'lowest h reached'}, 'Location','northwest', 'Box','off');

txt = sprintf('feasibility floor = %.3f\\times10^{-3}', mean(fl));
if all(abs(fl) < 1e-12)
    txt = sprintf(['feasibility floor = 0 at SoC = %.2f\n' ...
                   'the QP stays feasible to the barrier itself'], p.SoC_nom);
end
text(0.40, 0.90, txt, 'Units','normalized', 'FontSize',9, ...
     'Color',[0.7 0 0], 'VerticalAlignment','top');

ax2 = subplot(1,2,2);
Sg = linspace(p.SoC_min, p.SoC_max, 400);
fm = arrayfun(rep.h_min_of_SoC, Sg)*K;
plot(Sg, fm, 'r-', 'LineWidth',1.8); grid on; box on; hold on;
plot(xlim, [h_max h_max]*K, 'k:', 'LineWidth',1.2);
text(mean(xlim), h_max*K, ' full safe set', 'VerticalAlignment','bottom', ...
     'FontSize',9);
plot([p.SoC_enf_lo p.SoC_enf_lo], [0 h_max*K*1.1], '--', 'Color',[0.45 0.45 0.45]);
plot([p.SoC_enf_hi p.SoC_enf_hi], [0 h_max*K*1.1], '--', 'Color',[0.45 0.45 0.45]);
xlabel('state of charge');
ylabel('feasibility floor h_{min}  (10^{-3} p.u.)');
title('Floor rises only as the battery depletes');
ylim([0, h_max*K*1.15]);
xlim([p.SoC_min p.SoC_max]);

set([ax1 ax2],'FontSize',9);
end

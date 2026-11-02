function rep = verify_lie(nSamp, doPlots)
%VERIFY_LIE  Relative-degree and Lie-derivative checks on the barrier functions.

if nargin < 1 || isempty(nSamp),   nSamp   = 4000;  end
if nargin < 2 || isempty(doPlots), doPlots = false; end

spec = cpmg_spec();
p    = params_isolated(spec);
ix   = spec.ix;  nx = spec.nx;
rng(20260904,'twister');

cfg = struct('a_f',[8 8],'a_S',[2 2],'a_T',[0.5 0.5],'Delta_cbf',0);

Dbox = zeros(nx,2);
Dbox(ix.df,   :) = [-0.030 , 0.020];
Dbox(ix.dPtie,:) = [-0.300 , 0.300];
Dbox(ix.dPm,  :) = [ p.Pm_min , p.Pm_max ];
Dbox(ix.dPv,  :) = [ p.Pv_min , p.Pv_max ];
Dbox(ix.zACE, :) = [-5.000 , 5.000];
Dbox(ix.dPB,  :) = [ p.uB_min , p.uB_max ];
Dbox(ix.dSB,  :) = [ 0.10 , 1.00 ];
Dbox(ix.dPA,  :) = [ p.uA_min , p.uA_max ];
Dbox(ix.dTin, :) = [-6.0 , 6.0 ];

fields = {'f','Slo','Shi','Tlo','Thi'};
labels = {'freq','SoC_lo','SoC_hi','Tin_lo','Tin_hi'};
nB     = numel(fields);

hh   = 1e-6;
zero3 = zeros(spec.nd,1);

errLfh   = zeros(nB,1);  errLf2h = zeros(nB,1);  errLgLfh = zeros(nB,1);
maxLgh   = 0;            errPsi2 = 0;
nSkip    = 0;  nUsed = 0;
worstPhiGap = 0;

for s = 1:nSamp
    x  = Dbox(:,1) + (Dbox(:,2)-Dbox(:,1)).*rand(nx,1);
    ug = p.agc.u_min + (p.agc.u_max-p.agc.u_min)*rand;
    u  = [ p.uB_min + (p.uB_max-p.uB_min)*rand ;
           p.uA_min + (p.uA_max-p.uA_min)*rand ];

    B = barrier_functions(x, ug, zero3, p, 0, cfg);

    maxLgh = max(maxLgh, max(abs(B.rd.Lgh)));

    if straddlesKink(x, ug, zero3, p, hh)
        nSkip = nSkip + 1;
        continue
    end
    nUsed = nUsed + 1;

    [Fx, Gx] = cpmg_dynamics(x, ug, zero3, p, 0);

    for k = 1:nB
        bk = B.(fields{k});

        gh  = numGrad(@(z) hval(z,fields{k},p,ug,zero3,cfg), x, hh);
        errLfh(k) = max(errLfh(k), relerr(gh*Fx, bk.Lfh));

        gLf = numGrad(@(z) lfhval(z,fields{k},p,ug,zero3,cfg), x, hh);
        errLf2h(k) = max(errLf2h(k), relerr(gLf*Fx, bk.Lf2h));

        errLgLfh(k) = max(errLgLfh(k), relerr(gLf*Gx, bk.LgLfh));
    end

    k = 1;  a = cfg.a_f;
    gp1 = numGrad(@(z) psi1val(z,fields{k},p,ug,zero3,cfg), x, hh);
    lhs = gp1*(Fx + Gx*u) + a(2)*B.(fields{k}).psi1;
    rhs = B.(fields{k}).b + B.(fields{k}).LgLfh*u;
    errPsi2 = max(errPsi2, relerr(lhs, rhs));

    r_f = B.f.LgLfh(1)  /(1/(p.M*p.TB));
    r_S = B.Slo.LgLfh(1)/(-p.eta_B/(p.E_B_s*p.TB));
    worstPhiGap = max(worstPhiGap, abs(r_f - r_S));
end

rep = struct('tests',struct('id',{},'name',{},'value',{},'tol',{},'pass',{},'note',{}));
tol = 1e-5;

rep = add(rep,'L01','Lg h = 0 for all barriers', maxLgh, 1e-14, maxLgh<1e-14, ...
          'relative degree is at least 2 for every safe set');

for k = 1:nB
    rep = add(rep, sprintf('L02.%d',k), sprintf('Lf h  [%s]', labels{k}), ...
              errLfh(k), tol, errLfh(k)<tol, '');
end
for k = 1:nB
    rep = add(rep, sprintf('L03.%d',k), sprintf('Lf^2 h  [%s]', labels{k}), ...
              errLf2h(k), tol, errLf2h(k)<tol, '');
end
for k = 1:nB
    rep = add(rep, sprintf('L04.%d',k), sprintf('LgLf h  [%s]', labels{k}), ...
              errLgLfh(k), tol, errLgLfh(k)<tol, '');
end

x0 = p.x0;
B0 = barrier_functions(x0, 0, zero3, p, 0, cfg);
zS = max(abs(B0.Slo.LgLfh(2)), abs(B0.Shi.LgLfh(2)));   % storage vs u_A
zT = max(abs(B0.Tlo.LgLfh(1)), abs(B0.Thi.LgLfh(1)));   % comfort vs u_B
rep = add(rep,'L05','structural zeros in LgLf h', max(zS,zT), 1e-14, ...
          max(zS,zT)<1e-14, ...
          'storage is reachable from u_B only; comfort from u_A only');

rep = add(rep,'L06','psi2 = d(psi1)/dt + a2*psi1', errPsi2, 1e-4, errPsi2<1e-4, ...
          'HO-CBF recursion is internally consistent');

rep = add(rep,'L07','AIAC channel of LgLf h_f is negative', B0.f.LgLfh(2), 0, ...
          B0.f.LgLfh(2) < 0, ...
          sprintf(['LgLf h_f = [%+.4f , %+.4f]. A negative AIAC entry means ' ...
                   'frequency support requires u_A < 0 (load shedding), per Hui2019.'], ...
                  B0.f.LgLfh));

rep = add(rep,'L08','Phi enters freq and storage BESS channels identically', ...
          worstPhiGap, 1e-12, worstPhiGap<1e-12, ...
          'both barriers lose fast-actuator authority simultaneously as Phi -> 0');

xt = x0; xt(ix.df) = spec.safety.df_min_pu;  Bt = barrier_functions(xt,0,zero3,p,0,cfg);
e1 = abs(Bt.f.h);
xt = x0; xt(ix.dSB) = p.SoC_enf_lo;          Bt = barrier_functions(xt,0,zero3,p,0,cfg);
e2 = abs(Bt.Slo.h);
xt = x0; xt(ix.dTin) = p.Tin_max - p.Tin_ref; Bt = barrier_functions(xt,0,zero3,p,0,cfg);
e3 = abs(Bt.Thi.h);
rep = add(rep,'L09','h = 0 exactly on each set boundary', max([e1 e2 e3]), 1e-14, ...
          max([e1 e2 e3])<1e-14, '');

xt = x0; xt(ix.dSB) = p.SoC_enf_lo;
Bt = barrier_functions(xt,0,zero3,p,0,cfg);
authAtAnchor = abs(Bt.Slo.LgLfh(1))*p.uB_max;
xt = x0; xt(ix.dSB) = p.SoC_min;
Bp = barrier_functions(xt,0,zero3,p,0,cfg);
authAtPhys = abs(Bp.Slo.LgLfh(1))*p.uB_max;
anchorOK = p.Phi_at_anchor >= spec.safety.Phi_floor_at_anchor && ...
           authAtAnchor > 20*max(authAtPhys,eps);
rep = add(rep,'L12','storage barrier anchored where Phi is still usable', ...
          p.Phi_at_anchor, spec.safety.Phi_floor_at_anchor, anchorOK, ...
          sprintf(['Phi = %.3f at the enforced bound SoC = %.2f, giving barrier ' ...
                   'authority %.3e; at the PHYSICAL bound SoC = %.2f authority is ' ...
                   '%.3e. Anchoring at the physical bound would ask the barrier to ' ...
                   'defend a region the actuator cannot reach.'], ...
                  p.Phi_at_anchor, p.SoC_enf_lo, authAtAnchor, p.SoC_min, authAtPhys));

LfT = abs(B0.Thi.Lfh);
dT3 = LfT*3.0 + abs(B0.Thi.LgLfh(2))*p.uA_max*0.5*9.0;
rep = add(rep,'L10','comfort barrier magnitude on a 3 s horizon', dT3, Inf, true, ...
          sprintf(['|Lf h_T| = %.3e degC/s; reachable dTin over 3 s is %.3e degC ' ...
                   'against a %.1f degC margin. The comfort barrier CANNOT bind on ' ...
                   'a short DoS event. State this scope limit in Section IV-B.'], ...
                  LfT, dT3, p.Tin_max - p.Tin_ref));

frac = nSkip/nSamp;
rep = add(rep,'L11','fraction of samples on a non-smooth stencil', frac, 0.35, ...
          frac < 0.35, ...
          sprintf(['%d of %d samples straddled a dead-band, GRC or anti-windup ' ...
                   'kink and were excluded from the differentiability tests. ' ...
                   '%d samples used.'], nSkip, nSamp, nUsed));

rep.allPass = all([rep.tests.pass]);
rep.nUsed   = nUsed;
rep.nSkip   = nSkip;
printRep(rep);

if doPlots, plotLie(p,spec,cfg); end
if ~rep.allPass
    warning('verify_lie:FAILED', ...
        'Lie derivative verification FAILED. Do not proceed to Layer 4.');
end
end   % verify_lie

function v = hval(x, fld, p, ug, d, cfg)
B = barrier_functions(x, ug, d, p, 0, cfg);  v = B.(fld).h;
end

function v = lfhval(x, fld, p, ug, d, cfg)
B = barrier_functions(x, ug, d, p, 0, cfg);  v = B.(fld).Lfh;
end

function v = psi1val(x, fld, p, ug, d, cfg)
B = barrier_functions(x, ug, d, p, 0, cfg);  v = B.(fld).psi1;
end

function g = numGrad(f, x, h)
n = numel(x);  g = zeros(1,n);
for i = 1:n
    xp = x; xm = x; xp(i) = xp(i)+h; xm(i) = xm(i)-h;
    g(i) = (f(xp) - f(xm))/(2*h);
end
end

function e = relerr(a, b)
e = norm(a(:)-b(:))/max(1, norm(b(:)));
end

function tf = straddlesKink(x, ug, d, p, h)
n = numel(x);
[~,~,dg0] = cpmg_dynamics(x, ug, d, p, 0);
ref = [dg0.db_active, dg0.grc_active, dg0.lim_active];
tf  = false;
for i = 1:n
    for sgn = [-1 1]
        xs = x; xs(i) = xs(i) + sgn*h;
        [~,~,dgs] = cpmg_dynamics(xs, ug, d, p, 0);
        if ~isequal([dgs.db_active, dgs.grc_active, dgs.lim_active], ref)
            tf = true; return
        end
    end
end
if any(dg0.aw_sigma > 0 & dg0.aw_sigma < 1), tf = true; end
end

function rep = add(rep,id,name,value,tol,pass,note)
n = numel(rep.tests)+1;
rep.tests(n) = struct('id',id,'name',name,'value',value,'tol',tol, ...
                      'pass',logical(pass),'note',note);
end

function printRep(rep)
fprintf('\n');
fprintf('========================================================================\n');
fprintf(' LIE DERIVATIVE VERIFICATION REPORT\n');
fprintf('========================================================================\n');
fprintf(' %-7s %-44s %12s %-6s\n','ID','QUANTITY','REL ERR','RESULT');
fprintf('------------------------------------------------------------------------\n');
for k = 1:numel(rep.tests)
    T = rep.tests(k);
    if isfinite(T.tol) && T.tol > 0
        fprintf(' %-7s %-44s %12.3e %-6s\n', T.id, T.name, T.value, tern(T.pass));
    else
        fprintf(' %-7s %-44s %12s %-6s\n', T.id, T.name, '-', tern(T.pass));
    end
    if ~isempty(T.note), fprintf('         %s\n', T.note); end
end
fprintf('------------------------------------------------------------------------\n');
fprintf(' OVERALL: %s   (%d/%d)   samples used %d, skipped %d\n', ...
        tern(all([rep.tests.pass])), nnz([rep.tests.pass]), numel(rep.tests), ...
        rep.nUsed, rep.nSkip);
fprintf('========================================================================\n\n');
end

function s = tern(c)
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function plotLie(p,spec,cfg)
figure('Color','w','Name','Layer 3 diagnostics','Position',[100 100 950 400]);
S = linspace(0.15,0.95,400);
Af = zeros(size(S)); AS = zeros(size(S));
for k = 1:numel(S)
    x = p.x0; x(spec.ix.dSB) = S(k);
    B = barrier_functions(x,0,zeros(spec.nd,1),p,0,cfg);
    Af(k) = abs(B.f.LgLfh(1))*p.uB_max;
    AS(k) = abs(B.Slo.LgLfh(1))*p.uB_max;
end
subplot(1,2,1);
plot(S,Af,'LineWidth',1.8); grid on; box on; hold on;
yline(abs(cfg.a_f(1))*0,'k:');
xline(p.SoC_min,'r--'); xline(p.SoC_max,'r--');
xlabel('State of charge'); ylabel('BESS barrier authority (p.u./s^2)');
title('Frequency barrier: fast-actuator authority');
subplot(1,2,2);
plot(S,AS,'LineWidth',1.8); grid on; box on;
xline(p.SoC_min,'r--'); xline(p.SoC_max,'r--');
xlabel('State of charge'); ylabel('SoC barrier authority (1/s^2)');
title('Storage barrier: same \Phi, same collapse');
end

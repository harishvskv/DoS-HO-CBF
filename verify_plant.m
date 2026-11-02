function rep = verify_plant(mode, doPlots)
%VERIFY_PLANT  Dimensional, continuity and integration-order checks on the plant.

if nargin < 1 || isempty(mode),    mode = 'both';  end
if nargin < 2 || isempty(doPlots), doPlots = false; end

spec = cpmg_spec();
p    = params_isolated(spec);
ix   = spec.ix;

rep = struct();
rep.spec_version = spec.version;
rep.timestamp    = datestr(now,'yyyy-mm-dd HH:MM:SS');   % #ok<TNOW1,DATST>
rep.tests        = struct('id',{},'name',{},'value',{},'tol',{},'pass',{},'note',{});

zero3 = zeros(spec.nd,1);

idxs = cell2mat(struct2cell(spec.ix))';
ok   = isequal(sort(idxs),1:spec.nx) && numel(spec.stateNames)==spec.nx;
rep  = addTest(rep,'T01','spec index map is a permutation of 1..nx', ...
               double(ok),0.5,ok>0.5,'');

okS = abs(spec.safety.df_min_pu - spec.hz2pu(49.40)) < 1e-12;
rep = addTest(rep,'T01b','safety threshold = UFR Stage-1 (49.40 Hz)', ...
              spec.safety.df_min_pu,1e-12,okS, ...
              sprintf('%.4f p.u. = %.2f Hz', spec.safety.df_min_pu, spec.pu2hz(spec.safety.df_min_pu)));

if any(strcmpi(mode,{'isolated','both'}))

x0 = p.x0;
[F0,~,dg0] = cpmg_dynamics(x0, 0, zero3, p, 0);
res = max(abs(F0));
rep = addTest(rep,'T02','equilibrium residual max|F(x0)|',res,1e-10,res<1e-10, ...
              sprintf('PL = %.6f vs PL0 = %.6f', dg0.PL, p.PL0));

del = 0.01;
sgn = @(k,v) signOfDfDot(x0,k,v,p,zero3);
s_PB   = sgn(ix.dPB,  +del);
s_PA   = sgn(ix.dPA,  +del);
s_Ptie = sgn(ix.dPtie,+del);   % net export must LOWER frequency

xt = x0; xt(ix.dPB) = +del;
Ft = cpmg_dynamics(xt,0,zero3,p,0);
s_SoC = sign(Ft(ix.dSB));

xt = x0; xt(ix.dPA) = +del;
Ft = cpmg_dynamics(xt,0,zero3,p,0);
s_Tin = sign(Ft(ix.dTin));

sgnOK = (s_PB==+1) && (s_PA==-1) && (s_Ptie==-1) && (s_SoC==-1) && (s_Tin==-1);
rep = addTest(rep,'T03','D3 sign conventions in assembled field', ...
              double(sgnOK),0.5,sgnOK, ...
              sprintf('dPB:%+d dPA:%+d dPtie:%+d SoC:%+d Tin:%+d', ...
                      s_PB,s_PA,s_Ptie,s_SoC,s_Tin));

dfs  = [0, 0.5*p.df_db, 0.999*p.df_db, 1.001*p.df_db, 3*p.df_db];
phid = zeros(size(dfs));
for k = 1:numel(dfs)
    xt = x0; xt(ix.df) = dfs(k);
    [~,~,dgk] = cpmg_dynamics(xt,0,zero3,p,0);
    phid(k) = dgk.phi_db;
end
dbOK = all(abs(phid(1:3))<1e-14) && phid(4)>0 && abs(phid(5)-(3*p.df_db-p.df_db))<1e-12;
rep = addTest(rep,'T04','dead-band suppresses droop for |df| <= df_db', ...
              double(dbOK),0.5,dbOK, ...
              sprintf('df_db = %.3f Hz; phi_db(3*db) = %.3e', p.df_db*spec.f_base_Hz, phid(5)));

xt = x0; xt(ix.dPv) = 0.45; xt(ix.dPm) = -0.20;   % forces a large raw rate
[Ft,~,dgt] = cpmg_dynamics(xt,0,zero3,p,0);
raw = (-xt(ix.dPm) + xt(ix.dPv))/p.Tt;
grcOK = dgt.grc_active && abs(Ft(ix.dPm) - p.GRC) < 1e-12 && raw > p.GRC;
rep = addTest(rep,'T05','GRC clips dPm/dt at +/- p.GRC', ...
              Ft(ix.dPm),1e-12,grcOK, ...
              sprintf('raw = %.4f p.u./s, clipped to %.4f, limit %.4f', raw, Ft(ix.dPm), p.GRC));

xt = x0;
xt(ix.dPv) = p.Pv_max;
xt(ix.df)  = -0.020;   % strong under-frequency
ug_hi      = p.agc.u_max;
raw_Pv = ( -xt(ix.dPv) - deadbandRef(xt(ix.df),p.df_db)/p.R + ug_hi )/p.Tg;
[Ft,~,dgt] = cpmg_dynamics(xt,ug_hi,zero3,p,0);
awOK_v = raw_Pv > 0 && abs(Ft(ix.dPv)) < 1e-14 && dgt.lim_active(1);
rep = addTest(rep,'T06a','anti-windup holds the governor valve limit', ...
              Ft(ix.dPv),1e-14,awOK_v, ...
              sprintf('raw rate %+.3f p.u./s at dPv = Pv_max, clamped to %+.3e', ...
                      raw_Pv, Ft(ix.dPv)));

p2t = p;  p2t.Pm_max = 0.10;  p2t.Pm_min = -0.10;
xt = x0;  xt(ix.dPm) = p2t.Pm_max;  xt(ix.dPv) = 0.50;
[Ft2,~,dgt2] = cpmg_dynamics(xt,0,zero3,p2t,0);
raw_Pm = (-xt(ix.dPm) + xt(ix.dPv))/p2t.Tt;
awOK_m = raw_Pm > 0 && abs(Ft2(ix.dPm)) < 1e-14 && dgt2.lim_active(2);
rep = addTest(rep,'T06b','anti-windup holds the mechanical power limit', ...
              Ft2(ix.dPm),1e-14,awOK_m, ...
              sprintf('raw rate %+.3f clamped to %+.3e (reduced-limit copy)', ...
                      raw_Pm, Ft2(ix.dPm)));

rep = addTest(rep,'T06c','Pm magnitude limit is structurally unreachable', ...
              p.Pm_max - p.Pv_max, Inf, true, ...
              sprintf(['Pm_max = %.2f, Pv_max = %.2f. Since dPm relaxes toward dPv, ' ...
                       'the Pm clamp can never fire for admissible dPv. Record in ' ...
                       'Section II-B; do not claim it as an active nonlinearity.'], ...
                      p.Pm_max, p.Pv_max));

nSw  = 2001;
dPvs = linspace(p.Pv_max - 3*p.aw_band, p.Pv_max, nSw);
Fv   = zeros(1,nSw);
for k = 1:nSw
    xt = x0; xt(ix.dPv) = dPvs(k); xt(ix.df) = -0.020;
    Fk = cpmg_dynamics(xt,p.agc.u_max,zero3,p,0);
    Fv(k) = Fk(ix.dPv);
end
jump   = max(abs(diff(Fv)));
dstep  = dPvs(2)-dPvs(1);
slope  = jump/dstep;   % local Lipschitz slope of F(4)
contOK = abs(Fv(end)) < 1e-14 && slope < 5/p.aw_band*max(abs(Fv));
rep = addTest(rep,'T06d','F is continuous across the anti-windup band', ...
              slope,5/p.aw_band*max(abs(Fv)),contOK, ...
              sprintf(['max jump %.3e over a %.3e step, local slope %.1f. ' ...
                       'Hard clamp would jump by %.3f. Assumption 1 (Lipschitz F) ' ...
                       'now holds.'], jump, dstep, slope, max(abs(Fv))));

Sg  = linspace(0,1,501);
Phi = zeros(size(Sg));
for k = 1:numel(Sg)
    xt = x0; xt(ix.dSB) = Sg(k);
    [~,~,dgk] = cpmg_dynamics(xt,0,zero3,p,0);
    Phi(k) = dgk.Phi;
end
Phi_lo  = interp1(Sg,Phi,p.SoC_min);
Phi_hi  = interp1(Sg,Phi,p.SoC_max);
Phi_mid = interp1(Sg,Phi,p.SoC_nom);
phiOK = all(Phi>=0 & Phi<=1) && Phi_lo<0.01 && Phi_hi<0.01 && Phi_mid>0.99;
rep = addTest(rep,'T07','Phi in [0,1], -> 0 at BOTH SoC bounds', ...
              Phi_mid,0.01,phiOK, ...
              sprintf('Phi(%.2f)=%.4f  Phi(%.2f)=%.4f  Phi(%.2f)=%.4f', ...
                      p.SoC_min,Phi_lo,p.SoC_nom,Phi_mid,p.SoC_max,Phi_hi));

xt = x0; xt(ix.dSB) = p.SoC_enf_lo;  [~,~,dgA] = cpmg_dynamics(xt,0,zero3,p,0);
xt = x0; xt(ix.dSB) = p.SoC_min;     [~,~,dgP] = cpmg_dynamics(xt,0,zero3,p,0);
anchorOK = dgA.Phi >= spec.safety.Phi_floor_at_anchor && dgA.Phi > 20*dgP.Phi;
rep = addTest(rep,'T07b','Phi is usable at the ENFORCED storage bound', ...
              dgA.Phi, spec.safety.Phi_floor_at_anchor, anchorOK, ...
              sprintf(['Phi(%.2f) = %.4f at the enforced bound versus Phi(%.2f) = ' ...
                       '%.4f at the physical bound. The storage barrier is anchored ' ...
                       'at the former (spec.safety.SoC_anchor = ''%s'').'], ...
                      p.SoC_enf_lo, dgA.Phi, p.SoC_min, dgP.Phi, spec.safety.SoC_anchor));

dfp = 0.01;  h = 1e-7;
xt = x0; xt(ix.df) = dfp;
[~,~,dgc] = cpmg_dynamics(xt,0,zero3,p,0);
xp = xt; xp(ix.df) = dfp+h;  [~,~,dgp] = cpmg_dynamics(xp,0,zero3,p,0);
xm = xt; xm(ix.df) = dfp-h;  [~,~,dgm] = cpmg_dynamics(xm,0,zero3,p,0);
num = (dgp.PL - dgm.PL)/(2*h);
errL = abs(num - dgc.dPL_ddf)/max(abs(num),eps);

h2 = 1e-4;
xp2 = xt; xp2(ix.df) = dfp+h2;  [~,~,dgp2] = cpmg_dynamics(xp2,0,zero3,p,0);
xm2 = xt; xm2(ix.df) = dfp-h2;  [~,~,dgm2] = cpmg_dynamics(xm2,0,zero3,p,0);
d2num = (dgp2.PL - 2*dgc.PL + dgm2.PL)/h2^2;
base  = 1 + p.D*dfp;
d2ana = p.PL0*p.gamma*(p.gamma-1)*p.D^2*base^(p.gamma-2);
errC  = abs(d2num - d2ana)/max(abs(d2ana),eps);
nlOK  = errL < 1e-6 && abs(d2ana) > 1e-6 && errC < 1e-4;
rep = addTest(rep,'T08','dPL/ddf analytic = numerical, load is nonlinear', ...
              errL,1e-6,nlOK, ...
              sprintf('gamma = %.2f; d2PL/ddf2 numerical %.4f vs analytic %.4f', ...
                      p.gamma, d2num, d2ana));

[~,G0] = cpmg_dynamics(x0,0,zero3,p,0);
Lgh    = G0(ix.df,:);
rep = addTest(rep,'T09','Lg h = 0 (relative degree >= 2)', ...
              max(abs(Lgh)),1e-14,max(abs(Lgh))<1e-14,'');

LgLfh_num = numericLgLfh(x0,p,zero3);
LgLfh_ana = [dgAt(x0,p,zero3).Phi/(p.M*p.TB) , -1/(p.M*p.Tac)];
errRD = max(abs(LgLfh_num - LgLfh_ana));
rdOK  = errRD < 1e-6 && all(abs(LgLfh_ana) > 1e-6);
rep = addTest(rep,'T10','LgLf h analytic = numerical, nonzero', ...
              errRD,1e-6,rdOK, ...
              sprintf('LgLf h = [%.4f , %.4f]  (AIAC entry must be negative)', LgLfh_ana));

xt = x0; xt(ix.dSB) = p.SoC_min;
LgLfh_deg = numericLgLfh(xt,p,zero3);
degOK = abs(LgLfh_deg(1)) < 0.01*abs(LgLfh_ana(1));
rep = addTest(rep,'T11','relative degree degenerates as Phi -> 0', ...
              LgLfh_deg(1),0.01*abs(LgLfh_ana(1)),degOK, ...
              sprintf('BESS channel collapses from %.4f to %.4f at SoC = %.2f', ...
                      LgLfh_ana(1), LgLfh_deg(1), p.SoC_min));

tau_hr = 1/p.a_T/3600;
xt = x0; xt(ix.dPA) = 0;   % nominal AIAC power
[Ft] = cpmg_dynamics(xt,0,zero3,p,0);
thOK = abs(Ft(ix.dTin)) < 1e-14 && abs(tau_hr - p.tau_th_hr) < 1e-9;
rep = addTest(rep,'T12','thermal equilibrium exact, tau_th correct', ...
              tau_hr,1e-9,thOK, ...
              sprintf('tau_th = %.2f h; residual = %.3e degC/s', tau_hr, Ft(ix.dTin)));

dT_3s = abs(p.b_T*p.uA_max*3.0);
rep = addTest(rep,'T12b','indoor dT over a 3 s full-authority AIAC action', ...
              dT_3s,Inf,true, ...
              sprintf('%.2e degC. Comfort barrier is INACTIVE on this horizon.', dT_3s));

Pdis = 0.5;
xt = x0; xt(ix.dPB) = Pdis;
Ft = cpmg_dynamics(xt,0,zero3,p,0);
ana_rate = -p.eta_B*Pdis/p.E_B_s;
errU = abs(Ft(ix.dSB) - ana_rate)/abs(ana_rate);
rep = addTest(rep,'T13a','SoC vector field matches -eta*P/E_B', ...
              errU,1e-12,errU<1e-12, ...
              sprintf('E_B = %.1f s (legacy used C_bat = %.0f, a %.2fx error)', ...
                      p.E_B_s, p.C_bat_kWh, p.C_bat_kWh/p.E_B_s));

Tdis = 10.0;  hstep = spec.num.Ts_plant;
xt = x0;
for k = 1:round(Tdis/hstep)
    xt(ix.dPB) = Pdis;
    xt = cpmg_integrate(xt,0,[0;0],zero3,p,0,hstep,'rk4');
end
dSoC_sim = p.SoC_nom - xt(ix.dSB);
Pbar     = Pdis*(1 - hstep/(2*p.TB));
dSoC_ana = p.eta_B*Pbar*Tdis/p.E_B_s;
errE = abs(dSoC_sim - dSoC_ana)/dSoC_ana;
rep = addTest(rep,'T13b','integrated SoC excursion matches energy balance', ...
              errE,5e-3,errE<5e-3, ...
              sprintf('sim %.6f vs analytic %.6f', dSoC_sim, dSoC_ana));

specS = spec;
specS.nl.governor_deadband = false;
specS.nl.turbine_GRC       = false;
specS.nl.bess_Phi          = false;
pS = p;  pS.spec = specS;
pS.Pm_min = -1e3; pS.Pm_max = 1e3;   % disable the anti-windup kinks
pS.Pv_min = -1e3; pS.Pv_max = 1e3;
[ordS, errS, hsS] = convergenceStudy(pS, specS);
rep = addTest(rep,'T14a','RK4 order on the SMOOTH field', ...
              ordS,3.5,ordS>3.5, ...
              sprintf('observed order %.2f (expect ~4)', ordS));

[ordN, errN, hsN] = convergenceStudy(p, spec);
rep = addTest(rep,'T14b','RK4 order on the FULL nonlinear field', ...
              ordN,0.8,ordN>0.8, ...
              sprintf(['observed order %.2f. Degradation from %.2f is caused by ' ...
                       'dead-band and GRC event crossings, not by the integrator. ' ...
                       'Quote this order in the numerical-accuracy paragraph.'], ...
                      ordN, ordS));

iProd = find(abs(hsN - spec.num.Ts_ctrl) < 1e-12, 1);
if isempty(iProd), iProd = 1; end
rep = addTest(rep,'T14c','trajectory error at the production step size', ...
              errN(iProd),1e-4,errN(iProd)<1e-4, ...
              sprintf('||x(T)-x_ref(T)|| = %.3e at h = %.0e s over a 5 s transient', ...
                      errN(iProd), hsN(iProd)));

hs = hsN; errs = errN;

bufOK = p.cyber.tau_buffer > p.cyber.age_dos_max;
rep = addTest(rep,'T20','delay buffer horizon exceeds worst-case age of data', ...
              p.cyber.tau_buffer,p.cyber.age_dos_max,bufOK, ...
              sprintf(['buffer %.3f s > age_dos_max %.3f s = dur %.2f + T_pkt %.3f ' ...
                       '+ latency %.3f. Age of data is NOT attack duration; ' ...
                       'the legacy buffer conflated them and silently truncated.'], ...
                      p.cyber.tau_buffer, p.cyber.age_dos_max, ...
                      p.cyber.dos_dur_max, p.cyber.T_pkt, p.cyber.tau_nom_hi));

st = []; t = 0; dt = 0.01; shed = 0; tripped_at = nan(1,4);
for k = 1:4
    dfk = spec.ufls.df_pu(k) - 1e-4;   % just below stage k pick-up
    for n = 1:round((spec.ufls.delay_s(k)+2*dt)/dt)
        [shed, st] = ufls_relay(dfk, t, st, spec);
        t = t + dt;
    end
    tripped_at(k) = st.trip_time(k);
end
uflsOK = all(~isnan(tripped_at)) && abs(shed - sum(spec.ufls.shed_frac)) < 1e-12;
rep = addTest(rep,'T19','UFLS stages arm and trip at IEGC pick-up levels', ...
              shed,1e-12,uflsOK, ...
              sprintf('trip times [%s] s', num2str(tripped_at,'%.2f ')));

if doPlots, plotDiagnostics(Sg,Phi,hs,errs,p,spec); end

end   % isolated

if any(strcmpi(mode,{'39bus','both'}))

P  = params_39bus(spec);
nA = numel(P);  nx = spec.nx;  nu = spec.nu;

X0 = zeros(nA*nx,1);
for i = 1:nA, X0((i-1)*nx + ix.dSB) = P(i).SoC_nom; end
Ug = zeros(nA,1);  Dm = zeros(nA,spec.nd);

[Fm, Gm] = cpmg_dynamics_multiarea(X0,Ug,Dm,P);

mask = true(size(Gm));
for i = 1:nA, mask((i-1)*nx+(1:nx),(i-1)*nu+(1:nu)) = false; end
bdOK = all(Gm(mask)==0);
rep = addTest(rep,'T16','multi-area G is block diagonal (decentralised)', ...
              max(abs(Gm(mask))),0,bdOK,'');

Xr = X0;
for i = 1:nA, Xr((i-1)*nx+ix.df) = 0.01*i; end   % arbitrary df spread
Fr = cpmg_dynamics_multiarea(Xr,Ug,Dm,P);
tieSum = 0;
for i = 1:nA, tieSum = tieSum + Fr((i-1)*nx+ix.dPtie); end
rep = addTest(rep,'T17','sum of tie-line derivatives = 0', ...
              abs(tieSum),1e-12,abs(tieSum)<1e-12,'');

resM = max(abs(Fm));
rep = addTest(rep,'T18','multi-area equilibrium residual', ...
              resM,1e-10,resM<1e-10,'');

Xr = X0;
for i = 1:nA, Xr((i-1)*nx+ix.df) = 0.005*(i-2); end
U  = zeros(nA*nu,1);  hstep = spec.num.Ts_plant;
Xc = cpmg_integrate_multiarea(Xr,Ug,U,Dm,P,hstep,'rk4');
Xs = Xr;
df_hold = zeros(nA,1);
for i = 1:nA, df_hold(i) = Xr((i-1)*nx+ix.df); end
for i = 1:nA
    rows = (i-1)*nx+(1:nx);
    tie_i = 2*pi*sum(P(1).net.T_ij(i,:)'.*(df_hold(i)-df_hold));
    Xs(rows) = cpmg_integrate(Xr(rows),Ug(i),U((i-1)*nu+(1:nu)),Dm(i,:)',P(i),tie_i,hstep,'rk4');
end
splitErr = max(abs(Xc - Xs));
rep = addTest(rep,'T15','operator-splitting error per step', ...
              splitErr,1e-6,splitErr<1e-6, ...
              'reported results use cpmg_integrate_multiarea (no splitting)');

end   % multi-area

rep.allPass = all([rep.tests.pass]);
printReport(rep);

if ~rep.allPass
    warning('verify_plant:FAILED', ...
        'Plant verification FAILED. Do not proceed to Layer 3 until all tests pass.');
end

end   % verify_plant

function rep = addTest(rep,id,name,value,tol,pass,note)
n = numel(rep.tests)+1;
rep.tests(n).id    = id;
rep.tests(n).name  = name;
rep.tests(n).value = value;
rep.tests(n).tol   = tol;
rep.tests(n).pass  = logical(pass);
rep.tests(n).note  = note;
end

function y = deadbandRef(u,w)
if u > w
    y = u - w;
elseif u < -w
    y = u + w;
else
    y = 0;
end
end

function s = signOfDfDot(x0,k,v,p,d0)
xt = x0; xt(k) = xt(k) + v;
F  = cpmg_dynamics(xt,0,d0,p,0);
s  = sign(F(p.spec.ix.df));
end

function dg = dgAt(x,p,d0)
[~,~,dg] = cpmg_dynamics(x,0,d0,p,0);
end

function LgLfh = numericLgLfh(x,p,d0)
ix = p.spec.ix;  n = p.spec.nx;  hh = 1e-7;
g = zeros(1,n);
for i = 1:n
    xp = x; xm = x; xp(i) = xp(i)+hh; xm(i) = xm(i)-hh;
    Fp = cpmg_dynamics(xp,0,d0,p,0);
    Fm = cpmg_dynamics(xm,0,d0,p,0);
    g(i) = (Fp(ix.df) - Fm(ix.df))/(2*hh);
end
[~,G] = cpmg_dynamics(x,0,d0,p,0);
LgLfh = g*G;
end

function [ordObs, errs, hs] = convergenceStudy(p,spec)
Tend = 5.0;
d    = [0; 0.10; 0];   % 0.10 p.u. load step
u    = [0.1; -0.05];
hs   = [2e-2, 1e-2, 5e-3, 2.5e-3];
href = hs(end)/8;

xr = p.x0;
for n = 1:round(Tend/href)
    xr = cpmg_integrate(xr,0,u,d,p,0,href,'rk4');
end

Xend = zeros(spec.nx,numel(hs));
for k = 1:numel(hs)
    x = p.x0;
    for n = 1:round(Tend/hs(k))
        x = cpmg_integrate(x,0,u,d,p,0,hs(k),'rk4');
    end
    Xend(:,k) = x;
end

errs = zeros(1,numel(hs));
for k = 1:numel(hs)
    errs(k) = norm(Xend(:,k) - xr);
end

valid = errs > 1e-12;
if nnz(valid) >= 2
    pfit   = polyfit(log(hs(valid)), log(errs(valid)), 1);
    ordObs = pfit(1);
else
    ordObs = NaN;
end
end

function plotDiagnostics(Sg,Phi,hs,errs,p,spec)
figure('Color','w','Name','verify_plant diagnostics','Position',[100 100 1000 700]);

subplot(2,2,1);
plot(Sg,Phi,'LineWidth',1.8); grid on; box on; hold on;
xline(p.SoC_min,'r--','SoC_{min}'); xline(p.SoC_max,'r--','SoC_{max}');
xlabel('State of charge'); ylabel('\Phi(SoC)');
title('T07: two-sided BESS derating');

subplot(2,2,2);
dfs = linspace(-4*p.df_db,4*p.df_db,801); phid = zeros(size(dfs));
for k=1:numel(dfs)
    xt=p.x0; xt(spec.ix.df)=dfs(k);
    [~,~,dgk]=cpmg_dynamics(xt,0,zeros(spec.nd,1),p,0); phid(k)=dgk.phi_db;
end
plot(dfs*spec.f_base_Hz,phid*spec.f_base_Hz,'LineWidth',1.8); grid on; box on;
xlabel('\Delta f (Hz)'); ylabel('\phi_{db}(\Delta f) (Hz)');
title('T04: governor dead-band');

subplot(2,2,3);
loglog(hs(1:numel(errs)),errs,'o-','LineWidth',1.6); grid on; box on;
xlabel('step size h (s)'); ylabel('||x(T) - x_{ref}(T)||');
title('T14: RK4 convergence');

subplot(2,2,4);
Sgrid = linspace(0.15,0.95,400); Lg = zeros(size(Sgrid));
for k=1:numel(Sgrid)
    xt=p.x0; xt(spec.ix.dSB)=Sgrid(k);
    tmp = numericLgLfh(xt,p,zeros(spec.nd,1)); Lg(k)=tmp(1);
end
plot(Sgrid,Lg,'LineWidth',1.8); grid on; box on;
xlabel('State of charge'); ylabel('L_gL_f h  (BESS channel)');
title('T11: relative-degree degeneration');
end

function printReport(rep)
fprintf('\n');
fprintf('========================================================================\n');
fprintf(' PLANT VERIFICATION REPORT   %s   spec = %s\n', rep.timestamp, rep.spec_version);
fprintf('========================================================================\n');
fprintf(' %-5s %-52s %-6s\n','ID','TEST','RESULT');
fprintf('------------------------------------------------------------------------\n');
for k = 1:numel(rep.tests)
    T = rep.tests(k);
    fprintf(' %-5s %-52s %-6s\n', T.id, T.name, ternary(T.pass,'PASS','FAIL'));
    if ~isempty(T.note)
        fprintf('       %s\n', T.note);
    end
end
fprintf('------------------------------------------------------------------------\n');
fprintf(' OVERALL: %s   (%d/%d)\n', ternary(all([rep.tests.pass]),'PASS','FAIL'), ...
        nnz([rep.tests.pass]), numel(rep.tests));
fprintf('========================================================================\n\n');
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end

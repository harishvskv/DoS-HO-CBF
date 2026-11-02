function [h_crit, info] = compute_hcrit(p, dP, opt)
%COMPUTE_HCRIT  Critical initial barrier margin below which no admissible input prevents a breach.

if nargin < 3, opt = struct(); end
opt = setdef(opt,'tol',1e-4);
opt = setdef(opt,'Tend',6.0);
opt = setdef(opt,'verbose',false);
opt = setdef(opt,'qp',struct());

spec = p.spec;  ix = spec.ix;
Ts   = spec.num.Ts_plant;
Tc   = spec.num.Ts_ctrl;
h_max = -spec.safety.df_min_pu;

info = struct();
info.authority_exceeded = dP >= p.uB_max;

if info.authority_exceeded
    info.t_star = Inf;
    info.est_lag = Inf; info.est_slew = Inf; info.est_zoh = Inf;
    info.est_total = Inf;
    h_crit = Inf;
    info.frac = Inf;
    if opt.verbose, printH(h_crit, info, p, dP, h_max); end
    return
end

t_star        = -p.TB*log(1 - dP/p.uB_max);
info.t_star   = t_star;
info.est_lag  = ( (dP - p.uB_max)*t_star + dP*p.TB )/p.M;
info.est_slew = (dP/(2*p.M))*(dP/p.uB_rate);
info.est_zoh  = (dP/p.M)*Tc;
info.est_total = info.est_lag + info.est_slew + info.est_zoh;

holds = @(h0) barrierHeld(h0, dP, p, opt, Ts, Tc, ix, spec);

lo = 0;  hi = h_max;
if ~holds(hi)
    h_crit = Inf;   % even a full margin is insufficient
else
    while (hi - lo) > opt.tol
        mid = 0.5*(lo + hi);
        if holds(mid), hi = mid; else, lo = mid; end
    end
    h_crit = hi;
end
info.frac = h_crit/h_max;

if opt.verbose, printH(h_crit, info, p, dP, h_max); end
end   % compute_hcrit

function ok = barrierHeld(h0, dP, p, opt, Ts, Tc, ix, spec)
x = p.x0;  x(ix.df) = spec.safety.df_min_pu + h0;
d = [0; dP; 0];
cs = [];  u = [0;0];
cfg0 = struct('a_f',[8 8],'Delta_cbf',0);
nStep = round(opt.Tend/Ts);  kCtrl = round(Tc/Ts);
ok = true;
for n = 1:nStep
    if mod(n-1,kCtrl) == 0
        [~, uref] = ctrl_pi_agc(x, p, 'edge_only');
        F  = cpmg_dynamics(x, 0, d, p, 0);
        [u, cs] = ctrl_hocbf_qp(x, F(ix.df), uref, p, cs, opt.qp);
    end
    B = barrier_functions(x, 0, d, p, 0, cfg0);
    if B.f.h < 0, ok = false; return; end
    x = cpmg_integrate(x, 0, u, d, p, 0, Ts, 'rk4');
end
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

function printH(h_crit, info, p, dP, h_max)
fprintf('\n--- critical initial margin, dP = %.4f p.u. ---\n', dP);
if info.authority_exceeded
    fprintf('  imbalance exceeds u_B_max = %.2f p.u.: no initial margin suffices\n\n', p.uB_max);
    return
end
fprintf('  t* (BESS covers deficit) = %.4f s   [T_B = %.3f s]\n', info.t_star, p.TB);
fprintf('  loss  lag   %.5f\n', info.est_lag);
fprintf('  loss  slew  %.5f   [u_B_rate = %.1f p.u./s]\n', info.est_slew, p.uB_rate);
fprintf('  loss  ZOH   %.5f\n', info.est_zoh);
fprintf('  estimate    %.5f\n', info.est_total);
fprintf('  MEASURED h_crit = %.5f p.u.  (%.1f%% of the safe set %.4f)\n\n', ...
        h_crit, 100*info.frac, h_max);
end

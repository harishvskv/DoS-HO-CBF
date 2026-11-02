function [ug, u_edge, cs, ci] = ctrl_bank(arm, t, x, x_wan, rocof, net, p, ctx, cs)
%CTRL_BANK  Controller arms A0 to A6 for controlled comparison under identical actuator authority.

if isempty(cs)
    cs = struct('qp',[],'etm',[],'u_prev',[0;0], ...
                'x_last_tx',x,'n_tx',0,'sigma',0, ...
                'lam_prev',[],'u_ref_prev',[0;0],'bump',[0;0],'t_bump',-inf, ...
                'ug_prev',0);
end
ci = struct('arm',arm,'sigma',0,'qi',[],'mech','none','T_pred',NaN,'tau_pred',NaN);

switch upper(arm)

case 'A0'   % legacy baseline, no edge
    [ug, ~] = ctrl_pi_agc(x_wan, p, 'agc_only');
    ug     = net.alpha*ug;
    u_edge = [0;0];

case 'A1'   % baseline WITH equal authority
    ug     = net.alpha*ctrl_pi_agc(x_wan, p, 'agc_only');
    [~,ue] = ctrl_pi_agc(x, p, 'edge_only');
    u_edge = ue;

case 'A5'   % event-triggered PI-AGC
    sig_e = getfielddef(ctx,'etc_sigma',0.05);
    abs_e = getfielddef(ctx,'etc_abs',  1e-3);
    Tmax_e= getfielddef(ctx,'etc_Tmax', 1.0);
    if ~isfield(cs,'t_last_tx') || isempty(cs.t_last_tx), cs.t_last_tx = -inf; end

    io   = [p.spec.ix.df, p.spec.ix.dPtie];
    dev  = norm(x(io) - cs.x_last_tx(io));
    thr  = sig_e*max(norm(cs.x_last_tx(io)), 1e-9) + abs_e;
    fire = (dev > thr) || ((t - cs.t_last_tx) >= Tmax_e);
    if fire && net.alpha
        cs.x_last_tx = x;
        cs.t_last_tx = t;
        cs.n_tx      = cs.n_tx + 1;
    end
    [ug, ue] = ctrl_pi_agc(cs.x_last_tx, p, 'agc_edge');
    ug     = net.alpha*ug;
    u_edge = ue;
    ci.n_tx = cs.n_tx;

case 'A6'   % delay-compensated predictive AGC
    tau_p = min(net.age, p.cyber.tau_pred_max);
    xp    = x_wan;
    if tau_p > 0
        nP = max(1, ceil(tau_p/p.spec.num.Ts_plant));
        hP = tau_p/nP;
        dP_hat = [0; 0; 0];   % nominal: no disturbance forecast
        for q = 1:nP
            xp = cpmg_integrate(xp, cs.ug_prev, cs.u_prev, dP_hat, p, 0, hP, 'rk4');
        end
    end
    ug     = net.alpha*ctrl_pi_agc(xp, p, 'agc_only');   % central path compensated
    [~,ue] = ctrl_pi_agc(x,  p, 'edge_only');   % edge on LOCAL state, as A1
    u_edge = ue;
    cs.ug_prev = ug;
    ci.tau_pred = tau_p;

case {'A2','A3','A4'}

    if strcmpi(arm,'A2')
        sigma = 0;  ci.mech = 'disabled';
    else
        ectx = struct('a_d',ctx.a_d,'a_d_pred',ctx.a_d_pred, ...
                      'h_min',ctx.h_min,'h_max',ctx.h_max,'cfg',ctx.cfg, ...
                      'tie_in',getfielddef(ctx,'tie_in',0));
        [sigma, cs.etm, edi] = da_etm(t, x, rocof, net.pkt, net.age, p, ectx, cs.etm);
        ci.mech   = edi.mech;
        ci.T_pred = edi.T_pred;
    end
    cs.sigma = sigma;  ci.sigma = sigma;

    gs = struct('sigma', sigma);
    [ug_raw, u_agc] = ctrl_pi_agc(x_wan, p, 'agc_edge', gs);   % delayed, WAN side
    [~,      u_loc] = ctrl_pi_agc(x,     p, 'edge_only', gs);   % local, fresh

    if strcmpi(arm,'A2')
        lam = 1;
    else
        lam = net.alpha*(1 - sigma);
    end
    ug = net.alpha*ug_raw;

    u_raw = lam*u_agc + (1 - lam)*u_loc;
    if isempty(cs.lam_prev) || abs(lam - cs.lam_prev) > 1e-9
        cs.bump   = cs.u_ref_prev - u_raw;
        cs.t_bump = t;
    end
    dec   = exp(-(t - cs.t_bump)/p.etm.T_bump);
    u_ref = u_raw + cs.bump*dec;

    cs.lam_prev   = lam;
    cs.u_ref_prev = u_ref;

    switch upper(arm)
        case 'A3', Delta = 0;   % ablate Theorem 2
        otherwise, Delta = ctx.Delta_cbf;
    end

    qopt = ctx.qp;  qopt.Delta_cbf = Delta;
    qopt.tie_in   = getfielddef(ctx,'tie_in',0);
    [u_edge, cs.qp, qi] = ctrl_hocbf_qp(x, rocof, u_ref, p, cs.qp, qopt);
    ci.qi    = qi;
    ci.u_ref = u_ref;
    ci.lam   = lam;

otherwise
    error('ctrl_bank:arm','Unknown ablation arm "%s"', arm);
end

u_edge = [ min(max(u_edge(1), p.uB_min), p.uB_max) ;
           min(max(u_edge(2), p.uA_min), p.uA_max) ];
cs.u_prev = u_edge;
ci.u_edge = u_edge;
ci.ug     = ug;
end

function v = getfielddef(s, f, dflt)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

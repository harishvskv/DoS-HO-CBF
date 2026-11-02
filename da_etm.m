function [sigma, es, di] = da_etm(t, x, rocof, pkt_arrived, age, p, ctx, es)
%DA_ETM  Delay-aware event-triggering mechanism with predictive, packet-loss and age paths.

if isempty(es)
    es = struct('sigma',0,'t_iso',-inf,'t_last_pkt',t,'n_fresh',0, ...
                'n_trig',0,'mech',{{}},'t_trig',[],'mech_last','none');
end

ctx = setdef(ctx,'a_d',   1.7960);
ctx = setdef(ctx,'a_d_pred', ctx.a_d);
ctx = setdef(ctx,'tie_in', 0);
ctx = setdef(ctx,'h_min', 0);
ctx = setdef(ctx,'h_max', -p.spec.safety.df_min_pu);
ctx = setdef(ctx,'cfg',   struct('a_f',[8 8],'Delta_cbf',0));

d_hat = estimate_disturbance(x, rocof, p);
B     = barrier_functions(x, 0, d_hat, p, ctx.tie_in, ctx.cfg);
h     = B.f.h;
hdot  = B.f.Lfh;

T_pred = time_to_floor(h, hdot, ctx.a_d_pred, ctx.h_min);

if pkt_arrived
    es.t_last_pkt = t;
    es.n_fresh    = es.n_fresh + 1;
elseif (t - es.t_last_pkt) > 1.5*p.cyber.T_pkt
    es.n_fresh    = 0;
end
n_missed = floor((t - es.t_last_pkt)/p.cyber.T_pkt);

m1 = T_pred   <  p.etm.T_react;   % predictive, safety-critical
m2 = n_missed >= p.etm.n_miss;   % packet loss
m3 = age      >  p.etm.tau_trigger;   % age-of-data backstop

mech = 'none';
if es.sigma == 0
    if m1
        mech = 'predictive';
    elseif m2
        mech = 'packet_loss';
    elseif m3
        mech = 'age_of_data';
    end
    if ~strcmp(mech,'none')
        es.sigma     = 1;
        es.t_iso     = t;
        es.n_trig    = es.n_trig + 1;
        es.mech{end+1} = mech;
        es.t_trig(end+1) = t;
        es.mech_last = mech;
    end
else
    dwell_ok  = (t - es.t_iso) >= p.etm.T_dwell;
    fresh_ok  = es.n_fresh     >= p.etm.n_ok;
    hi_ok     = h      >= p.etm.h_recover*ctx.h_max;
    Tp_ok     = T_pred >= p.etm.k_recover*p.etm.T_react;
    age_ok    = age    <= p.etm.tau_trigger;
    margin_ok = hi_ok && Tp_ok && age_ok;
    if dwell_ok && fresh_ok && margin_ok
        es.sigma = 0;
        mech = 'reconnect';
        es.mech_last = mech;
    end
end

sigma = es.sigma;

di = struct('mech',mech,'T_pred',T_pred,'h',h,'hdot',hdot,'w_hat',d_hat(1), ...
            'n_missed',n_missed,'n_fresh',es.n_fresh, ...
            'age',age,'m1',m1,'m2',m2,'m3',m3, ...
            'dwell',t - es.t_iso, ...
            'rc_dwell', (t - es.t_iso) >= p.etm.T_dwell, ...
            'rc_fresh', es.n_fresh >= p.etm.n_ok, ...
            'rc_h',     h      >= p.etm.h_recover*ctx.h_max, ...
            'rc_Tpred', T_pred >= p.etm.k_recover*p.etm.T_react, ...
            'rc_age',   age    <= p.etm.tau_trigger);
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

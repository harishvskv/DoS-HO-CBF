function m = metrics(out, p)
%METRICS  Frequency, barrier and energy metrics for one closed-loop run.

spec = p.spec;  ix = spec.ix;
t  = out.t(:);   dt = t(2) - t(1);
df = out.X(ix.df,:)';
ptie = out.X(ix.dPtie,:)';
ACE  = p.beta*df + ptie;

m.freq.IAE  = trapz(t, abs(df));
m.freq.ISE  = trapz(t, df.^2);
m.freq.ITAE = trapz(t, t.*abs(df));
m.freq.nadir_pu = min(df);
m.freq.nadir_Hz = spec.pu2hz(min(df));
m.freq.peak_pu  = max(df);

band = spec.band.df_lo_pu;   % IEGC normal band, performance only
inb  = abs(df) <= abs(band);
k    = find(~inb, 1, 'last');
if isempty(k), m.freq.t_settle = 0; else, m.freq.t_settle = t(min(k+1,numel(t))); end

m.ace.IAE  = trapz(t, abs(ACE));
m.ace.ISE  = trapz(t, ACE.^2);
m.ace.ITAE = trapz(t, t.*abs(ACE));

h = out.h(:);
m.safety.min_h          = min(h);
m.safety.violation      = min(h) < 0;
m.safety.guarantee_void = any(~out.safe_guaranteed);
m.safety.guarded        = any(out.stage > 0);
m.safety.silent_failure = m.safety.guarded && m.safety.violation && ...
                          ~m.safety.guarantee_void;
m.safety.unguarded_viol = ~m.safety.guarded && m.safety.violation;
m.safety.t_unsafe       = dt*nnz(h < 0);
m.safety.ufls_trip      = any(~isnan(out.relay.trip_time));
m.safety.shed_final     = out.shed(end);
m.safety.min_h_Hz       = spec.pu2hz(min(h) + spec.safety.df_min_pu);

st = out.stage(:);
m.safety.stage_frac = arrayfun(@(k) nnz(st==k)/max(nnz(st>0),1), 1:4);

PB = out.X(ix.dPB,:)';  PA = out.X(ix.dPA,:)';
m.effort.E_bess  = trapz(t, abs(PB));
m.effort.E_aiac  = trapz(t, abs(PA));
m.effort.peak_uB = max(abs(out.U(1,:)));
m.effort.peak_uA = max(abs(out.U(2,:)));
m.effort.sat_frac = nnz( abs(out.U(1,:)) >= p.uB_max - 1e-9 | ...
                         abs(out.U(2,:)) >= p.uA_max - 1e-9 )/numel(t);

S = out.X(ix.dSB,:)';
m.soc.min = min(S);  m.soc.max = max(S);  m.soc.final = S(end);
m.soc.excursion = max(S) - min(S);
m.soc.hit_enforced = min(S) <= p.SoC_enf_lo || max(S) >= p.SoC_enf_hi;
m.soc.min_Phi = min(out.Phi);

Tin = p.Tin_ref + out.X(ix.dTin,:)';
dev = max( max(Tin - p.Tin_max, 0), max(p.Tin_min - Tin, 0) );
m.comfort.J_comf = trapz(t, dev);
m.comfort.max_dev_C = max(dev);
m.comfort.band_violation = any(dev > 0);
m.comfort.Tin_range = [min(Tin) max(Tin)];

m.cyber.iso_frac = mean(out.sigma > 0);
m.cyber.max_age  = max(out.age);
m.cyber.mean_age = mean(out.age);
m.cyber.dos_frac = mean(~out.alpha);
if isfield(out,'pkt_seen_frac'), m.cyber.pkt_seen_frac = out.pkt_seen_frac; end

ts = out.solve_time(~isnan(out.solve_time));
if isempty(ts)
    m.qp = struct('mean_solve_ms',NaN,'max_solve_ms',NaN,'n_solve',0);
else
    m.qp = struct('mean_solve_ms',1e3*mean(ts),'max_solve_ms',1e3*max(ts), ...
                  'n_solve',numel(ts));
end

m.arm = out.arm;
end

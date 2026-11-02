function [d_hat, info] = estimate_disturbance(x, rocof, p)
%ESTIMATE_DISTURBANCE  Bounded estimate of the unmeasured power and temperature disturbance.

spec = p.spec;  ix = spec.ix;
eps0 = spec.num.eps_guard;

df    = x(ix.df);
dPtie = x(ix.dPtie);
dPm   = x(ix.dPm);
dPB   = x(ix.dPB);
dPA   = x(ix.dPA);

if spec.nl.load_frequency
    base    = max(1 + p.D*df, eps0);
    PL_loc  = p.PL0*base^p.gamma;
else
    PL_loc  = p.PL0 + p.D*df;
end

measured = dPm + dPB - dPA - dPtie - (PL_loc - p.PL0);

w = p.M*rocof - measured;

d_hat = zeros(spec.nd,1);
d_hat(spec.id.dPres) = w;

info = struct('w',w,'measured_terms',measured,'PL_local',PL_loc,'rocof',rocof);
end

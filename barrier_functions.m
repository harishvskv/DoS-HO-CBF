function B = barrier_functions(x, ug, d, p, tie_in, cfg)
%BARRIER_FUNCTIONS  Barrier functions and Lie derivatives for the frequency, storage and thermal sets.

if nargin < 5 || isempty(tie_in), tie_in = 0; end
if nargin < 6, cfg = struct(); end
cfg = setdef(cfg,'a_f',[8 8]);
cfg = setdef(cfg,'a_S',[2 2]);
cfg = setdef(cfg,'a_T',[0.5 0.5]);
cfg = setdef(cfg,'Delta_cbf',0);
cfg = setdef(cfg,'enable',true(1,5));

spec = p.spec;  ix = spec.ix;

[F, G, dg] = cpmg_dynamics(x, ug, d, p, tie_in);

dPtie = x(ix.dPtie);
dPm   = x(ix.dPm);
dPB   = x(ix.dPB);
SoC   = x(ix.dSB);
dPA   = x(ix.dPA);
dTin  = x(ix.dTin);

Phi     = dg.Phi;
dPm_dot = dg.dPm_dot;   % GRC- and clamp-limited, exact
bT_eff  = p.b_T*(dg.COP_eff/p.COP);

Bf.h    = x(ix.df) - spec.safety.df_min_pu;
Bf.Lfh  = F(ix.df);
Bf.Lf2h = ( dPm_dot ...
          + (-dPB/p.TB) ...
          - (-dPA/p.Tac) ...
          - tie_in ...
          - dg.dPL_ddf*F(ix.df) )/p.M;
Bf.LgLfh = [ Phi/(p.M*p.TB) , -1/(p.M*p.Tac) ];

cS = p.eta_B/p.E_B_s;

switch lower(spec.safety.SoC_anchor)
    case 'effective'
        S_lo_enf = p.SoC_enf_lo;  S_hi_enf = p.SoC_enf_hi;
    case 'physical'
        S_lo_enf = p.SoC_min;     S_hi_enf = p.SoC_max;
    otherwise
        error('barrier_functions:socAnchor','Unknown SoC_anchor "%s"', ...
              spec.safety.SoC_anchor);
end

BSlo.h     = SoC - S_lo_enf;
BSlo.Lfh   = -cS*dPB;
BSlo.Lf2h  = +cS*dPB/p.TB;
BSlo.LgLfh = [ -cS*Phi/p.TB , 0 ];

BShi.h     = S_hi_enf - SoC;
BShi.Lfh   = +cS*dPB;
BShi.Lf2h  = -cS*dPB/p.TB;
BShi.LgLfh = [ +cS*Phi/p.TB , 0 ];

BSlo.h_phys = SoC - p.SoC_min;
BShi.h_phys = p.SoC_max - SoC;
BSlo.anchor = S_lo_enf;
BShi.anchor = S_hi_enf;

Tin_abs = p.Tin_ref + dTin;

BTlo.h     = Tin_abs - p.Tin_min;
BTlo.Lfh   = F(ix.dTin);
BTlo.Lf2h  = -p.a_T*F(ix.dTin) + bT_eff*dPA/p.Tac;
BTlo.LgLfh = [ 0 , -bT_eff/p.Tac ];

BThi.h     = p.Tin_max - Tin_abs;
BThi.Lfh   = -F(ix.dTin);
BThi.Lf2h  = +p.a_T*F(ix.dTin) - bT_eff*dPA/p.Tac;
BThi.LgLfh = [ 0 , +bT_eff/p.Tac ];

Bf   = hocbf(Bf  , cfg.a_f);
BSlo = hocbf(BSlo, cfg.a_S);
BShi = hocbf(BShi, cfg.a_S);
BTlo = hocbf(BTlo, cfg.a_T);
BThi = hocbf(BThi, cfg.a_T);

B.f = Bf;  B.Slo = BSlo;  B.Shi = BShi;  B.Tlo = BTlo;  B.Thi = BThi;

all_b  = {Bf, BSlo, BShi, BTlo, BThi};
names  = {'CBF:freq','CBF:SoC_lo','CBF:SoC_hi','CBF:Tin_lo','CBF:Tin_hi'};
tighten = [cfg.Delta_cbf, 0, 0, 0, 0];

A = []; bq = []; nm = {};
for k = 1:5
    if ~cfg.enable(k), continue; end
    A  = [A ; -all_b{k}.LgLfh];   % #ok<AGROW>
    bq = [bq;  all_b{k}.b - tighten(k)];   % #ok<AGROW>
    nm = [nm ; names(k)];   % #ok<AGROW>
end
B.A = A;  B.bq = bq;  B.names = nm;
B.dg = dg;
B.cfg = cfg;

B.rd.Lgh       = G(ix.df,:);   % must be exactly [0 0]
B.rd.LgLfh     = [Bf.LgLfh ; BSlo.LgLfh ; BShi.LgLfh ; BTlo.LgLfh ; BThi.LgLfh];
B.rd.names     = names;
B.rd.Phi       = Phi;
B.rd.Phi_degen = dg.Phi_degen;
B.rd.bess_authority = abs(Bf.LgLfh(1))*max(abs([p.uB_min p.uB_max]));
B.rd.aiac_authority = abs(Bf.LgLfh(2))*max(abs([p.uA_min p.uA_max]));
B.rd.soc_anchor     = [BSlo.anchor, BShi.anchor];
B.rd.soc_physical   = [p.SoC_min,   p.SoC_max];
B.rd.h_soc_phys     = [BSlo.h_phys, BShi.h_phys];

end   % barrier_functions

function s = hocbf(s, a)
s.a    = a;
s.psi1 = s.Lfh + a(1)*s.h;
s.b    = s.Lf2h + (a(1)+a(2))*s.Lfh + a(1)*a(2)*s.h;
end

function c = setdef(c, f, v)
if ~isfield(c,f) || isempty(c.(f)), c.(f) = v; end
end

function P = params_39bus(spec)
%PARAMS_39BUS  Parameter set and tie-line matrix for the three-area system.

if nargin < 1 || isempty(spec), spec = cpmg_spec(); end

nA = 3;

M_a     = [10.0 ,  8.0 ,  6.0 ];
D_a     = [ 1.0 ,  1.2 ,  0.9 ];
Tt_a    = [ 0.30,  0.33,  0.35];   % turbine time constant       [s]
Tg_a    = [ 0.08,  0.07,  0.07];   % governor time constant      [s]
R_a     = [ 0.05,  0.05,  0.05];
gam_a   = [ 1.5 ,  1.5 ,  1.5 ];   % load exponent                [-]
PL0_a   = [ 0.90,  0.72,  1.08];   % nominal load                [p.u.]
Sbase_a = [ 100 ,  100 ,  100 ];   % area power base             [MW]
GRCm_a  = [ 10.0,  10.0,  10.0];   % GRC                    [p.u./min]
Cbat_a  = [5000 , 5000 , 5000 ];   % BESS energy capacity        [kWh]

for i = nA:-1:1
    p = struct();
    p.name       = sprintf('IEEE 39-bus Area %d', i);
    p.params_version = 'params-2.6';
    p.spec       = spec;
    p.areaIndex  = i;
    p.nAreas     = nA;
    p.isolated   = false;

    p.S_base_MW = Sbase_a(i);
    p.S_base_kW = p.S_base_MW*1e3;
    p.f_base_Hz = spec.f_base_Hz;

    p.M     = M_a(i);
    p.D     = D_a(i);
    p.gamma = gam_a(i);
    p.PL0   = PL0_a(i);
    p.Pref  = PL0_a(i);

    p.Tg = Tg_a(i);
    p.Tt = Tt_a(i);
    p.R  = R_a(i);

    p.df_db        = 0.03/p.f_base_Hz;   % +/- 0.03 Hz dead-band
    p.db_smooth_k  = 2e3;
    p.GRC_pu_per_min = GRCm_a(i);
    p.GRC          = p.GRC_pu_per_min/60;
    p.grc_smooth_k = 1.0;
    p.Pv_min = -0.50;  p.Pv_max = 0.50;
    p.aw_band = 0.02;

    p.Pm_min = -0.50;  p.Pm_max = 0.50;

    p.TB        = 0.10;
    p.C_bat_kWh = Cbat_a(i);
    p.eta_B     = 0.95;
    p.E_B_s     = p.C_bat_kWh*3600/p.S_base_kW;   % [s] -- see (E1) note

    p.SoC_nom  = spec.safety.SoC_nom;
    p.SoC_min  = spec.safety.SoC_min;
    p.SoC_max  = spec.safety.SoC_max;
    p.SoC_marg = 0.05;
    p.SoC_lo   = p.SoC_min + p.SoC_marg;
    p.SoC_hi   = p.SoC_max - p.SoC_marg;
    p.ks       = 60.0;
    p.Phi_min  = 1e-3;
    p.SoC_enf_lo = p.SoC_lo;
    p.SoC_enf_hi = p.SoC_hi;
    phiAt = @(S) 0.25*(1+tanh(p.ks*(S-p.SoC_lo)))*(1-tanh(p.ks*(S-p.SoC_hi)));
    p.Phi_at_anchor = min(phiAt(p.SoC_enf_lo), phiAt(p.SoC_enf_hi));
    assert(p.Phi_at_anchor >= spec.safety.Phi_floor_at_anchor, ...
        ['params: Phi at the enforced storage bound is %.3f, below the required ' ...
         'floor %.3f. Widen SoC_marg or reduce ks.'], ...
         p.Phi_at_anchor, spec.safety.Phi_floor_at_anchor);
    p.uB_min   = -0.60;  p.uB_max = 0.60;
    p.uB_rate  = 12.0;
    p.uB_full_traverse_s = (p.uB_max - p.uB_min)/p.uB_rate;

    p.Tac = 0.20;
    p.COP = 2.50;
    p.tau_th_hr = 4.0;
    p.a_T       = 1/(p.tau_th_hr*3600);
    p.Tout_nom  = 32.0;
    p.Tin_ref   = spec.safety.Tin_ref_C;
    p.PA_nom    = 0.20;
    p.k_th_pu = (p.Tout_nom - p.Tin_ref)/p.PA_nom;
    p.b_T     = p.k_th_pu*p.a_T;
    p.Rth_eq_CperkW  = p.k_th_pu/(p.COP*p.S_base_kW);
    p.Cth_eq_kWhperC = p.tau_th_hr/p.Rth_eq_CperkW;

    p.COP_dTnom = 10.0;  p.COP_slope = 0.02;
    p.COP_min   = 1.5;   p.COP_max   = 4.0;

    p.Tin_min  = spec.safety.Tin_min_C;
    p.Tin_max  = spec.safety.Tin_max_C;
    p.uA_min   = -0.30;  p.uA_max = 0.30;
    p.uA_rate  = 0.10;
    p.uA_full_traverse_s = (p.uA_max - p.uA_min)/p.uA_rate;

    p.beta = 1/p.R + p.D;   % natural frequency bias [p.u./p.u.]

    p.agc.Kp = 0.40;
    p.agc.Ki = 0.50;
    p.agc.u_min = -0.50;  p.agc.u_max = 0.50;

    p.agc.Kp_e  = 0.50;
    p.agc.Ki_e  = 0.00;
    p.agc.K_soc = 0.50;
    p.agc.T_soc_restore = p.E_B_s/(p.eta_B*max(p.agc.K_soc,eps));

    p.cyber.T_pkt      = 0.050;
    p.cyber.tau_pred_max = 1.000;
    p.cyber.tau_nom_lo = 0.020;   % nominal WAN latency, lower  [s]
    p.cyber.tau_nom_hi = 0.040;   % nominal WAN latency, upper  [s]

    p.cyber.age_healthy_max = p.cyber.T_pkt + p.cyber.tau_nom_hi;

    p.cyber.dos_dur_max = 3.500;

    p.cyber.age_dos_max = p.cyber.dos_dur_max + p.cyber.T_pkt + p.cyber.tau_nom_hi;
    p.cyber.tau_dos_hi  = p.cyber.age_dos_max;   % retained name, same quantity

    p.cyber.tau_buffer = 1.5*p.cyber.age_dos_max;
    assert(p.cyber.tau_buffer > p.cyber.age_dos_max, ...
        'params: delay buffer horizon must exceed the worst-case age of data');
    assert(p.cyber.age_healthy_max < p.cyber.tau_nom_hi + p.cyber.T_pkt + eps, ...
        'params: healthy age of data is inconsistent with T_pkt and tau_nom_hi');

    p.etm.tau_trigger = 0.120;
    p.etm.n_miss      = 2;
    p.etm.T_react     = 0.150;   % predictive-trigger horizon [s]
    p.etm.T_dwell     = 1.000;
    p.etm.n_ok        = 3;
    p.etm.h_recover   = 0.70;
    p.etm.k_recover   = 3.0;
    p.etm.T_rocof     = 0.020;   % RoCoF measurement filter lag [s]
    p.etm.T_bump      = 2.000;
    assert(p.cyber.age_healthy_max < p.etm.tau_trigger, ...
        'params: trigger would fire on healthy jitter');

    x0 = zeros(spec.nx,1);
    x0(spec.ix.dSB) = p.SoC_nom;
    p.x0 = x0;

    P(i) = p;   % #ok<AGROW>
end

net = struct();
net.nAreas = nA;

T = [ 0.00 , 0.20 , 0.25 ;
      0.20 , 0.00 , 0.15 ;
      0.25 , 0.15 , 0.00 ];
assert(isequal(T,T'), 'params_39bus: T_ij must be symmetric');
assert(all(diag(T)==0), 'params_39bus: T_ii must be zero');
net.T_ij = T;

net.dist_share = [1.00 , 0.80 , 1.20];

net.gidx = @(area,local) (area-1)*spec.nx + local;

P(1).net = net;
for i = 2:nA, P(i).net = net; end

end

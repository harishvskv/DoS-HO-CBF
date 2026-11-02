function p = params_isolated(spec)
%PARAMS_ISOLATED  Parameter set for the isolated single-area microgrid.

if nargin < 1 || isempty(spec), spec = cpmg_spec(); end
p = struct();
p.name   = 'Isolated nonlinear CPMG';

p.params_version = 'params-2.6';
p.spec   = spec;
p.nAreas = 1;
p.isolated = true;

p.S_base_MW  = 10.0;
p.S_base_kW  = p.S_base_MW*1e3;
p.f_base_Hz  = spec.f_base_Hz;

p.M      = 8.0;
p.D      = 1.2;
p.gamma  = 1.5;
p.PL0    = 0.90;
p.Pref   = 0.85;

p.Tg     = 0.20;
p.Tt     = 0.30;
p.R      = 0.04;

p.df_db      = 0.03/p.f_base_Hz;   % = 6.0e-4 p.u.  (+/- 0.03 Hz)
p.db_smooth_k = 2e3;   % steepness of the smooth surrogate

p.GRC_pu_per_min = 20.0;
p.GRC = p.GRC_pu_per_min/60;   % [p.u./s]
p.grc_smooth_k = 1.0;

p.Pv_min = -0.50;   p.Pv_max = 0.50;   % valve travel        [p.u.]
p.aw_band = 0.02;

p.Pm_min = -0.50;   p.Pm_max = 0.50;   % mechanical power    [p.u.]

p.TB       = 0.10;
p.C_bat_kWh = 500.0;
p.eta_B    = 0.95;

p.E_B_s = p.C_bat_kWh*3600/p.S_base_kW;   % = 180 s

p.SoC_nom   = spec.safety.SoC_nom;
p.SoC_min   = spec.safety.SoC_min;
p.SoC_max   = spec.safety.SoC_max;
p.SoC_marg  = 0.05;   % knee margin
p.SoC_lo    = p.SoC_min + p.SoC_marg;
p.SoC_hi    = p.SoC_max - p.SoC_marg;
p.ks        = 60.0;   % Phi steepness
p.Phi_min   = 1e-3;   % degeneracy detection floor

p.SoC_enf_lo = p.SoC_lo;
p.SoC_enf_hi = p.SoC_hi;
phiAt = @(S) 0.25*(1+tanh(p.ks*(S-p.SoC_lo)))*(1-tanh(p.ks*(S-p.SoC_hi)));
p.Phi_at_anchor = min(phiAt(p.SoC_enf_lo), phiAt(p.SoC_enf_hi));
assert(p.Phi_at_anchor >= spec.safety.Phi_floor_at_anchor, ...
    ['params: Phi at the enforced storage bound is %.3f, below the required ' ...
     'floor %.3f. Widen SoC_marg or reduce ks.'], ...
     p.Phi_at_anchor, spec.safety.Phi_floor_at_anchor);

p.uB_min = -0.60;  p.uB_max = 0.60;   % [p.u.]
p.uB_rate  = 12.0;
p.uB_full_traverse_s = (p.uB_max - p.uB_min)/p.uB_rate;

p.Tac    = 0.20;
p.COP    = 2.50;

p.tau_th_hr = 4.0;
p.a_T = 1/(p.tau_th_hr*3600);   % [1/s]

p.Tout_nom = 32.0;   % nominal outdoor temperature [degC]
p.Tin_ref  = spec.safety.Tin_ref_C;   % comfort set point           [degC]
p.PA_nom   = 0.20;

p.k_th_pu = (p.Tout_nom - p.Tin_ref)/p.PA_nom;   % [degC per p.u.]
p.b_T     = p.k_th_pu*p.a_T;   % [degC/(s*p.u.)]

p.Rth_eq_CperkW  = p.k_th_pu/(p.COP*p.S_base_kW);   % [degC/kW]
p.Cth_eq_kWhperC = p.tau_th_hr/p.Rth_eq_CperkW;   % [kWh/degC]

p.COP_dTnom = 10.0;
p.COP_slope = 0.02;
p.COP_min   = 1.5;  p.COP_max = 4.0;

p.Tin_min   = spec.safety.Tin_min_C;
p.Tin_max   = spec.safety.Tin_max_C;

p.uA_min = -0.30;  p.uA_max = 0.30;
p.uA_rate  = 0.10;
p.uA_full_traverse_s = (p.uA_max - p.uA_min)/p.uA_rate;

p.beta = 1/p.R + p.D;
p.T_ij = 0;   % no synchronising coefficients

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
p.etm.T_rocof     = 0.020;
p.etm.T_bump      = 2.000;
assert(p.cyber.age_healthy_max < p.etm.tau_trigger, ...
    'params: trigger would fire on healthy jitter');

ix = spec.ix;
x0 = zeros(spec.nx,1);
x0(ix.df)    = 0.0;
x0(ix.dPtie) = 0.0;
x0(ix.dPm)   = 0.0;
x0(ix.dPv)   = 0.0;
x0(ix.zACE)  = 0.0;
x0(ix.dPB)   = 0.0;
x0(ix.dSB)   = p.SoC_nom;   % ABSOLUTE SoC, not a deviation
x0(ix.dPA)   = 0.0;
x0(ix.dTin)  = 0.0;   % deviation from Tin_ref
p.x0 = x0;

end

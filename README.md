# DSC-SEC-HO-CBF

**Decentralized Safety-Critical Edge-Control for Cyber-Physical Microgrids Under Denial-of-Service Attacks: A Time-Delay High-Order Control Barrier Function Approach**

MATLAB implementation accompanying manuscript **SEGAN-D-26-02975**, *Sustainable Energy, Grids and Networks*.

This repository contains the complete code for every result reported in the manuscript: the nonlinear plant, the barrier-certificate machinery, the delay-aware event-triggering mechanism, the controller bank, the verification suite, and the study drivers that produce each table and figure.

No proprietary toolbox is required. There is no YALMIP, no MOSEK and no CVX anywhere in the code. The quadratic programs are solved with `quadprog` (Optimization Toolbox), with a dense active-set fallback for installations without it.

---

## 1. Scope of the model

The plant is **nonlinear** and is never linearised. No small-signal surrogate is substituted for it in any study. The nonlinearities retained in the vector field are

- governor dead-band (GDB) on the speed-governor input,
- generation rate constraint (GRC) on the turbine power derivative,
- composite frequency-dependent load with non-unity exponent `gamma`,
- state-of-charge dependent derating `Phi(S)` of the storage inverter,
- actuator saturation and rate limits on both edge channels.

Each area carries the nine-state vector

```
x = [ df  dPtie  dPm  dPv  z  dPB  S  dPA  dTin ]'   in R^9
```

Under-frequency load shedding is modelled as a staged relay **external** to the vector field (`ufls_relay.m`), so that the drift remains continuous where the barrier analysis requires it.

Two configurations are provided: the isolated single-area microgrid, and the three-area interconnected system with tie-line coupling used for the multi-area validation.

---

## 2. Requirements

| Item | Notes |
|---|---|
| MATLAB | R2023b or later |
| Optimization Toolbox | required, for `quadprog` |
| Statistics and Machine Learning Toolbox | optional, exact binomial confidence intervals |

Developed and run on Windows 11, Intel Core i7-12700 (2.10 GHz), 15.6 GB usable RAM, matching the environment stated in Section V of the manuscript.

---

## 3. Reproducing the results

```matlab
addpath(genpath(pwd));
S = run_all();          % verification gates, then every study
make_figures(S);        % writes figures/, one plot per file
```

Partial execution:

```matlab
S = run_all('verify');                              % Layers 0-6, the gates only
S = run_all('studies', struct('S', S, 'n_mc', 50)); % Layer 7, reusing S
S = run_all('all', struct('useCache', true));       % reuse cached constants
```

`run_all` halts at the **first** failed gate and returns `S.gate_failed` naming it. Each layer consumes the outputs of the previous one, so running a study against a stale Lipschitz or contraction-rate struct produces numbers that look plausible and are wrong. The driver makes the dependency order explicit; do not call the study scripts out of order.

Indicative wall-clock times are printed in the header of `run_all.m`. A complete cold run, including the Lipschitz sweep, is of the order of fifty minutes.

---

## 4. Scope of the reported evidence

Forward invariance is **conditional**. It holds while the quadratic program remains feasible, while the age of data stays within the certified bound, and while the initial margin exceeds the critical value. Below that margin a breach is unavoidable for any admissible input, and the controller declares it rather than preventing it.

The stochastic study comprises **fifty** independent trials. The proposed scheme records no violation in those fifty trials, which gives an exact ninety-five per cent Clopper-Pearson interval of [0.00, 7.11] per cent on the underlying violation rate. An observed zero over fifty trials is not evidence of a zero underlying failure probability, and no claim in this repository or the manuscript should be read that way.

The three-area study is additional validation of the formulation in a multi-area setting with tie-line coupling. It is not a demonstration of scalability to large networks.

---

## 5. Citation

Please cite the published version of the manuscript. A BibTeX entry will be added on acceptance.

## 6. Contact

V.S.K.V. Harish, Department of Electrical Engineering, Netaji Subhas University of Technology, Dwarka, New Delhi 110078, India. `vskvharish@ieee.org`

## 7. License

BSD 3-Clause. Copyright (c) 2026, V.S.K.V. Harish. See `LICENSE`.

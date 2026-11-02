function [u, cs, qi] = ctrl_hocbf_qp(x, rocof, u_ref, p, cs, opt)
%CTRL_HOCBF_QP  High-order barrier quadratic program with a Lyapunov objective and slack.

spec = p.spec;  ix = spec.ix;

if nargin < 6, opt = struct(); end
opt = setdef(opt,'Delta_cbf',0);
opt = setdef(opt,'a_f',[8 8]);
opt = setdef(opt,'a_S',[2 2]);
opt = setdef(opt,'a_T',[0.5 0.5]);
opt = setdef(opt,'clf',[4 4]);
opt = setdef(opt,'w',[10 500]);
opt = setdef(opt,'p_delta',1e5);
opt = setdef(opt,'enable',true(1,5));
opt = setdef(opt,'use_rate',true);
opt = setdef(opt,'Ts',spec.num.Ts_ctrl);
opt = setdef(opt,'soft_penalty',1e4);
opt = setdef(opt,'tie_in',0);

if isempty(cs)
    cs = struct('u_prev',[0;0],'n_stage',zeros(1,4),'n_solve',0);
end

tSolve = tic;

d_hat = estimate_disturbance(x, rocof, p);
cfg   = struct('a_f',opt.a_f,'a_S',opt.a_S,'a_T',opt.a_T, ...
               'Delta_cbf',opt.Delta_cbf,'enable',opt.enable);
B = barrier_functions(x, 0, d_hat, p, opt.tie_in, cfg);

nRows = size(B.A,1);
Acbf  = [B.A , zeros(nRows,1)];   % barriers do not touch the slack
bcbf  = B.bq;

c1 = opt.clf(1);  c2 = opt.clf(2);
z1 = x(ix.df);
z2 = B.f.Lfh + c1*z1;
V  = 0.5*z1^2 + 0.5*z2^2;

Aclf = [ z2*B.f.LgLfh , -1 ];
bclf = -c2*V - z1*(z2 - c1*z1) - z2*( B.f.Lf2h + c1*(z2 - c1*z1) );

W = diag(opt.w);
H = blkdiag(W, opt.p_delta);
f = [ -W*u_ref(:) ; 0 ];
H = (H + H')/2;   % symmetry guard for quadprog

lb = [ p.uB_min ; p.uA_min ; 0   ];
ub = [ p.uB_max ; p.uA_max ; Inf ];

if opt.use_rate
    dmax = [p.uB_rate; p.uA_rate]*opt.Ts;
    lb(1:2) = max(lb(1:2), cs.u_prev - dmax);
    ub(1:2) = min(ub(1:2), cs.u_prev + dmax);
    lb(1:2) = min(lb(1:2), ub(1:2));   % guard against an empty box
end

qopt = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');

stage = 0;  delta = 0;  soft = zeros(1,4);  exitflag = -99;  z = [];

A1 = [Acbf; Aclf];  b1 = [bcbf; bclf];
[z, ~, exitflag] = quadprog(H, f, A1, b1, [], [], lb, ub, [], qopt);
if exitflag == 1, stage = 1; end

if stage == 0
    [z, ~, exitflag] = quadprog(H, f, Acbf, bcbf, [], [], lb, ub, [], qopt);
    if exitflag == 1, stage = 2; end
end

if stage == 0
    keep = strcmp(B.names,'CBF:freq');
    nS   = nnz(~keep);
    if nS > 0
        Hs = blkdiag(H, opt.soft_penalty*eye(nS));
        fs = [f; zeros(nS,1)];
        Ah = [B.A(keep,:) , zeros(nnz(keep),1) , zeros(nnz(keep),nS)];
        As = [B.A(~keep,:), zeros(nS,1)        , -eye(nS)];
        As2 = [Ah; As];
        bs2 = [B.bq(keep); B.bq(~keep)];
        lbs = [lb; zeros(nS,1)];
        ubs = [ub; inf(nS,1)];
        [zs, ~, exitflag] = quadprog(Hs, fs, As2, bs2, [], [], lbs, ubs, [], qopt);
        if exitflag == 1
            stage = 3;
            z = zs(1:3);
            soft(1:min(nS,4)) = zs(4:min(3+nS,7))';
        end
    end
end

if stage == 0
    stage = 4;
    g = B.f.LgLfh(:);
    ubest = zeros(2,1);
    for j = 1:2
        if g(j) >= 0, ubest(j) = ub(j); else, ubest(j) = lb(j); end
    end
    z = [ubest; 0];
    exitflag = 0;
end

u     = z(1:2);
delta = z(3);

psi2 = B.f.b + B.f.LgLfh*u;   % must be >= Delta_cbf
sat  = [ abs(u(1) - lb(1)) < 1e-9 || abs(u(1) - ub(1)) < 1e-9 , ...
         abs(u(2) - lb(2)) < 1e-9 || abs(u(2) - ub(2)) < 1e-9 ];

qi = struct( ...
    'stage',            stage, ...
    'exitflag',         exitflag, ...
    'delta',            delta, ...
    'soft',             soft, ...
    'h',                B.f.h, ...
    'psi1',             B.f.psi1, ...
    'psi2',             psi2, ...
    'psi2_margin',      psi2 - opt.Delta_cbf, ...
    'V',                V, ...
    'LgLfh',            B.f.LgLfh, ...
    'Phi',              B.dg.Phi, ...
    'Phi_degen',        B.dg.Phi_degen, ...
    'sat',              sat, ...
    'w_hat',            d_hat(spec.id.dPres), ...
    'Delta_cbf',        opt.Delta_cbf, ...
    'safe_guaranteed',  stage <= 3 && psi2 >= opt.Delta_cbf - 1e-9, ...
    'rows',             {B.names}, ...
    'solve_time',       toc(tSolve));

cs.u_prev  = u;
cs.n_stage(stage) = cs.n_stage(stage) + 1;
cs.n_solve = cs.n_solve + 1;
end

function s = setdef(s,f,v)
if ~isfield(s,f) || isempty(s.(f)), s.(f) = v; end
end

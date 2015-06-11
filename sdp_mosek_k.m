function [u, U, obj] = sdp_mosek_k(Hi, A, k, V) 

[B, ~] = size(A);
 
% building the linear cost matrix
C = [[-A, Hi/2]; [Hi'/2, 0]];
n = norm(C , 'fro');
C = C / n;

% SDP cost
prob.bardim = B + 1;
[subk, subl, val] = find(tril(C));
prob.barc.subj = ones(1, length(subk));
prob.barc.subk = subk';
prob.barc.subl = subl';
prob.barc.val  = val';

% lower bounds for the constraints
prob.blc = [ones(1, B+1), 2*(2*k-V)];
prob.buc = [ones(1, B+1), 2*(2*k-V)];

% no linear constraints on conic variables
prob.a = sparse(B+2, 0);

% adding the diagonal constraint
dsubi = 1:(B+1);
dsubj = ones(1, B+1);
dsubk = 1:(B+1);
dsubl = 1:(B+1);
dval  = ones(1, B+1);

% adding the linear constraint
lsubi = (B+2) * ones(1, B); 
lsubj = ones(1, B);
lsubl = 1:B;
lsubk = (B+1) * ones(1, B);
lval  = ones(1, B);

% adding all constraints together
prob.bara.subi = [dsubi, lsubi];
prob.bara.subj = [dsubj, lsubj];
prob.bara.subk = [dsubk, lsubk];
prob.bara.subl = [dsubl, lsubl];
prob.bara.val  = [dval,  lval ];

params = [];
params.MSK_IPAR_LOG = 0;
if isdeployed()
    params.MSK_IPAR_NUM_THREADS = 1;
else
    params.MSK_IPAR_NUM_THREADS = 4;
end

% solving the problem
[~, res] = mosekopt('maximize info echo(0)', prob, params);

% getting the SDP result
[ii, jj, ~] = find(tril(ones(B+1)));
idx = (ii + (B+1)*(jj-1));
X = zeros(B+1);
X(idx) = res.sol.itr.barx;
X = X + tril(X, -1)';

U = X(1:B, 1:B);
u = X(B+1, 1:B);

obj = res.sol.itr.pobjval * n;

end

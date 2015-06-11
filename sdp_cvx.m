function [u, U, obj] = sdp_cvx(H, A)

[B, ~] = size(A);

% building the linear cost matrix
C = [[-A, H/2]; [H'/2, 0]];
n = norm(C , 'fro');
C = C / n;

cvx_begin quiet
    variable M(B+1, B+1);
    maximize(trace(C * M));
    subject to
        M == semidefinite(B+1);
        diag(M) == ones(B+1, 1);
cvx_end

% getting back the u and U from M
U = M(1:B, 1:B);
u = M(B+1, 1:B);

% res-scaling the objective
obj = cvx_optval * n;

end
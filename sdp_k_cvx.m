function [u, U, obj] = sdp_k_cvx(Hi, A, k, V)

[B, ~] = size(A);

% building the linear cost matrix
C = [[-A, Hi/2]; [Hi'/2, 0]];
n = norm(C , 'fro');
C = C / n;

% indicator variables to get sum of u_y=k
I1 = zeros(B+1, 1);
I2 = zeros(B+1, 1);
I1(1:B) = 1;
I2(end) = 1;

cvx_begin quiet
    variable M(B+1, B+1);
    maximize(trace(C * M));
    subject to
        M == semidefinite(B+1);
        diag(M) == ones(B+1, 1);
        I1' * M * I2 == 2 * k - V;
cvx_end

% getting back the u and U from M
U = M(1:B, 1:B);
u = M(B+1, 1:B);

% res-scaling the objective
obj = cvx_optval * n;

end
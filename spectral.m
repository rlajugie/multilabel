function [u, U, obj] = spectral(H, A)

B   = size(A,1);
B1  = (1 / B);

% norm of A
nrm_A = norm(A, 'fro');

% normalizing A and H
A = A ./ sqrt(nrm_A);
H = H ./ sqrt(nrm_A);

C = [A, -eye(B); - B1 * 0.25 * (H * H'), A];
lambdas = eig(C, 'nobalance');

[lambda, idx]  = min(real(lambdas));

u = (A - lambda * eye(B)) \ (H/2);
u = u';

u(isnan(u)) = 0;

obj = u * H - u * A * u';
obj = obj * sqrt(nrm_A);

U = u' * u;

end
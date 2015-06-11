function [ u, U, obj ] = sdp_lowrank( H, A )
%SDP_LOWRANK Summary of this function goes here
%   Detailed explanation goes here

[B, ~] = size(A);

% building the linear cost matrix
C = [[-A, H/2]; [H'/2, 0]];
n = norm(C , 'fro');
C = C / n;

N = size(C, 1);
fun_set=@functions_elliptope;
fun_obj=@functions_max_cut;
param{1} = -C;
y0 = feval(fun_set,'retraction', randn(N,2), zeros(N,2),param);
[y,obj] = low_rank_optim(fun_set,fun_obj,param,y0);

obj = - obj * n;

M = y * y';

U = M(1:B, 1:B);
u = M(B+1, 1:B);

end
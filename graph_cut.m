function [ u, U, obj ] = graph_cut( b, A )
%GRAPH_CUT Summary of this function goes here
%   Detailed explanation goes here

C = - triu(A, 1);

c = 0.5 * b;

V = length(c);

[i,j,s] = find(C);
z = zeros(size(s,1),1);
D = [i,j,z,s,s,z];  % [i,j,e00,e01,e10,e11]

h = BK_Create();
BK_AddVars(h, V);


BK_SetUnary(h, [c'; zeros(1, V)]);
BK_SetPairwise(h, D);
e = BK_Minimize(h);

res = BK_GetLabeling(h);

BK_Delete(h);

u = 2*(double(res')-1)-1;

U = u' * u;

obj = u * b - u * A * u';

end
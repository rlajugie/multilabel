function [ us, Us, objective ] = optimize_f1( xi, yi, W, A, params )
%OPTIMIZE_OVER_K Summary of this function goes here
%   Detailed explanation goes here

V = length(yi);
u = cell(V+1, 1);
U = cell(V+1, 1);
obj = zeros(V+1, 1);

k = 0:V;

for i = 1:V+1
    Hi = W' * xi' - yi' / ( V + (2*k(i)) + sum(yi));
    
    if numel(A) > 1
        switch params.relaxation
            case 'sdp'
                switch params.solver
                    case 'cvx'
                        [u{i}, U{i}, obj(i)] = sdp_k_cvx(Hi, A, k(i), V);
                    case 'mosek'
                        [u{i}, U{i}, obj(i)] = sdp_k_mosek(Hi, A, k(i), V);
                end
            case {'spectral', 'graph-cut'}
                [ u{i}, U{i}, obj(i) ] = no_need_for_lagrange(Hi, A, k(i), V, params);
        end
    else
        u{i}   = sign(Hi');
        U{i}   = u' * u;
        obj(i) = compute_hamming(u, yi) + u * W' * xi'; 
    end
    
    obj(i) = obj(i) + V / ( V + (2 * k(i)) + sum(yi));
end

[objective, id] = max(obj);

us = u{id};
Us = U{id};

end

function [ u, U, obj ] = optimize_hamming( xi, yi, W, A, params )

V = length(yi);
H = W' * xi' - yi' / (2 * V);

if numel(A) > 1
    switch params.relaxation
        case 'graph-cut'
            [u, U, obj] = graph_cut(H, A);
        case 'sdp'
            switch params.solver
                case 'cvx'
                    [u, U, obj] = sdp_cvx(H, A);
                case 'mosek'
                    [u, U, obj] = sdp_mosek(H, A);
                case 'low-rank'
                    [u, U, obj] = sdp_lowrank(H, A);
            end
        case 'spectral'
            [u, U, obj] = spectral(H, A);
    end
else
    u   = sign(H');
    U   = u' * u;
    obj = compute_hamming(u, yi) + u * W' * xi'; 
end

obj = obj + 0.5;


end
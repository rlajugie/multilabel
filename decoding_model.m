function y_pred = decoding_model(W, A, x, params)

% getting the size of the problem
[N, d]  = size(x);
[~, V]  = size(W);

y_pred  = zeros(N, V);

fprintf('performing the decoding and rounding...\n');

for i = 1:N
    
    if mod(i, 100)==0
        fprintf('computing the decoding for example %8d/%8d\n', i, N);
    end
    
    % selecting an example
    xi  = x(i, :);
    H   = W' * xi';
    
    % solving the decoding with sdp because it is tighter than spectral
    switch params.relaxation
        case 'graph-cut'
            [u, U, obj] = graph_cut(H, A);
        otherwise
            switch params.solver
                case 'cvx'
                    [u, U, obj] = sdp_cvx(H, A);
                case 'mosek'
                    [u, U, obj] = sdp_mosek(H, A);
                case 'low-rank'
                    [u, U, obj] = sdp_lowrank(H, A);
            end
    end
    
    % performing the rounding
    y_pred(i, :) = rounding_sdp(U, u, A, W, xi, V, params);
    
end

end
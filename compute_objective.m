function [obj, o_spec, o_mi] = compute_objective(W, A, X, Y, ...
    lambda_w, lambda_a, params)

% This function computes the objective function of the paper

[N, ~] = size(X);
[~, V] = size(Y);

o_reg = lambda_w / 2 * trace(W'*W) + lambda_a / 2 * trace(A' * A);

o_mi    = zeros(1, N);
o_spec  = zeros(1, N);


for i = 1:N
    
    if mod(i, 100)==0
        fprintf('computing objective - sample %8d/%8d\n', i, N);
    end
    
    % picking example i
    xi = X(i, :);
    yi = Y(i, :);
    
    switch params.loss
        case 'hamming'
            [~, ~, o_spec(i)] = optimize_hamming(xi, yi, W, A, params);
        case 'f1'
            [~, ~, o_spec(i)] = optimize_f1(xi, yi, W, A, params);
    end
    
    o_mi(i) = yi * W' * xi' - yi * A * yi';
    
    obj = o_reg + sum(o_spec) / N - sum(o_mi) / N;
    
end

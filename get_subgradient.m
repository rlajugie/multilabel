function [dW, dA] = get_subgradient(xi, yi, W, A, lambda_w, lambda_a, params)

switch params.loss
    case 'hamming'
        [u, U] = optimize_hamming(xi, yi, W, A, params);
    case 'f1'
        [u, U] = optimize_f1(xi, yi, W, A, params);
end

dW = lambda_w * W + xi' * (u - yi);

if numel(A) > 1
    switch params.relaxation
        case 'sdp'
            dA = lambda_a * A - U + yi' * yi;
        case {'spectral', 'graph-cut'}
            dA = lambda_a * A - u' * u + yi' * yi;
    end
else
    dA = 0;
end

end
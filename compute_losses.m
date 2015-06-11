function [f1, r, p, h, y_pred] = compute_losses(W, A, X, Y, params)

if numel(A) > 1
    y_pred  = decoding_model(W, A, X, params);
else     
    y_pred  = sign(X*W);
end

[f1, r, p]  = compute_f1(y_pred, Y);
h           = compute_hamming(y_pred, Y);
end
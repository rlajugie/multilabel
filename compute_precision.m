function prec = compute_precision(y, y_gt)

prec = diag(y * y_gt') ./ diag(y * y');

end
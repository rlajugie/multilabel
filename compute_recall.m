function rec = compute_recall(y, y_gt)

rec = diag(y * y_gt') ./ diag(y_gt * y_gt');

end
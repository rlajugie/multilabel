function h = compute_hamming(y, y_gt)

V = size(y, 2);
h = 1/ (2 * V) * (V - diag(y * y_gt'));
h = mean(h);

end
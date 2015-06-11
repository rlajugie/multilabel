function [ f1, r, p ] = compute_f1( y, y_gt )
%COMPUTE_F1 Summary of this function goes here
%   Detailed explanation goes here

% Input labels are in {-1, 1}, changing to {0,1}
y       = (y + 1) / 2;
y_gt    = (y_gt + 1) / 2;

r = compute_recall(y, y_gt);
p = compute_precision(y, y_gt);
r(isnan(r)) = 1;
p(isnan(p)) = 1;

% computing the f1 score
f1_score = 1 - 2 * (r .* p) ./ (r + p);
f1_score(isnan(f1_score)) = 1;

% getting the mean
f1  = mean(f1_score);
r   = mean(r);
p   = mean(p);

end


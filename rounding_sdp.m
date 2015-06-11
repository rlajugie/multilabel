function y_pred = rounding_sdp(U, u, A, W, x, V, params)

% initialize m
m = -Inf;

% getting the error covariance
U1      = U - u' * u;
R       = mySqrtm(U1);

for i = 1:params.max_trials
    
    % getting u with normal distribution
    u_rounded = sign(u + randn(size(u))*R);
    
    % computing the corresponding objective
    y = u_rounded(1:V);
    score = y * W' * x' - u_rounded * A * u_rounded';
    
    % storing the max
    if score > m
        y_pred = y;
        m = score;
    end
    
end

end
function [gamma_a, gamma_w] = get_stepsize(lambda_w, lambda_a, gamma0_w, gamma0_a, t, params)

switch params.botou
    case 0
        gamma_w = 1 / ( 1 + lambda_w * (t));
        gamma_a = 1 / ( 1 + lambda_a * (t));
    case 1
        gamma_w = gamma0_w / ( gamma0_w + lambda_w * (t));
        gamma_a = gamma0_a / ( gamma0_a + lambda_a * (t));
    case 2
        gamma_w = gamma0_w / ...
            ( gamma0_w + lambda_w * (t)) .^ 0.75;
        gamma_a = gamma0_a / ...
            ( gamma0_a + lambda_a * (t)) .^ 0.75;
end
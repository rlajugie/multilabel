function [gamma_a, gamma_w] = get_stepsize(lambda_w, lambda_a, gamma0_w, gamma0_a, t, params)

gamma_w = gamma0_w * ((gamma0_w + lambda_w * t) ^ params.be);
gamma_a = gamma0_a * ((gamma0_a + lambda_a * t) ^ params.be);

end
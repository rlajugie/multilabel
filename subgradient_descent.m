function [outputs] = subgradient_descent( lambda_w, lambda_a, g_w, g_a, params )
%Implementation of the plain stochastic subgradient descent version for the code.
% T is the number of outer loop iteration in CCCP algorithm
% J is the number of subgradient step
% set Hamming to 1 to use the Hamming loss instead of the F1 score

% if the inner loop is graph-cut only negative A is possible
if strcmp(params.relaxation, 'graph-cut')
    params.proj_A = 'negative';
end

% if we use a low-rank sdp optimization use the hamming loss
if strcmp(params.solver, 'low-rank')
    params.loss = 'hamming';
end

% getting the gamma0 parameter
gamma0_w = g_w * lambda_w;
gamma0_a = g_a * lambda_a;

% setting the mosek license file
if strcmp(params.solver,'mosek')
    setenv('MOSEKLM_LICENSE_FILE', params.mosek_license);
end

% setting the random seed
rng(params.seed);

% loading the data
[x_train, x_val, x_test, y_train, y_val, y_test] = get_data(params);

% getting the size of the problem
[N, ~] = size(x_train);

% storage of the objective, train f1 and hamming
objective = inf(1, round(params.T / params.time));

% Intialization
[W, A] = initialize(x_train, y_train, x_val, y_val, lambda_w, params);

% Initialization of the counter
k           = 1;
mean_time   = 0;
params.broken = false;





fprintf('Optimization in W and A started...\n');

for t = 1:params.T
    
    % report average iteration duration every 1000 iterations
    if mod(t, 1000)==0
        fprintf('%8d mean-time=%8.4f\n', t, mean_time);
    end
    
    % compute the losses over whole dataset only every params.time iterations
    if mod(t, params.time)==1
        
        objective(k) = compute_objective(W, A, x_train, y_train, ...
            lambda_w, lambda_a, params);
        
        [f1, r, p, h, ~]        = compute_losses(W, A, x_train, y_train, params);
        [f1v, rv, pv, hv, ~]    = compute_losses(W, A, x_val, y_val, params);
        
        fprintf('iter=%05d obj=%5.3e \n', t, objective(k));
        fprintf('TRAIN : f1=%4.2f r=%4.2f p=%4.2f acc=%4.2f\n', f1, r, p, h);
        fprintf('VAL   : f1=%4.2f r=%4.2f p=%4.2f acc=%4.2f\n', f1v, rv, pv, hv);
        
        k = k + 1;
    end
    
    tic;
    
    % select the random example
    i   = randi(N);
    xi  = x_train(i, :);
    yi  = y_train(i, :);
    
    % get subgradient
    [dW, dA] = get_subgradient(xi, yi, W, A, lambda_w, lambda_a, params);
    
    % get the stepsize
    [gamma_a, gamma_w] = get_stepsize(lambda_w, lambda_a, gamma0_w, gamma0_a, t, params);
    
    % update W
    W = W - gamma_w * dW;
    A = A - gamma_a * dA;
    
    % project the matrix A
    A = project_A(A, params);
    
    % if A or W blow up, stop the opitmization
    if norm(W)>1e15 || norm(A)>1e15
        break;
    end
    
    % compute average iteration time
    mean_time = ((t-1) * mean_time + toc) / t;
    
end




% computing the val and test loss
[f1_train, ~, ~, h_train, y_pred_train] = compute_losses(W, A, x_train, y_train, params);
[f1_val, ~, ~, h_val, y_pred_val]       = compute_losses(W, A, x_val, y_val, params);
[f1_test, ~, ~, h_test, y_pred_test]    = compute_losses(W, A, x_test, y_test, params);

% putting the outputs in a structure
outputs.obj         = objective;
outputs.W           = W;
outputs.A           = A;

outputs.f1_train    = f1_train;
outputs.h_train     = h_train;
outputs.y_pred_train = y_pred_train;

outputs.f1_val      = f1_val;
outputs.h_val       = h_val;
outputs.y_pred_val  = y_pred_val;

outputs.f1_test     = f1_test;
outputs.h_test      = h_test;
outputs.y_pred_test = y_pred_test;

end
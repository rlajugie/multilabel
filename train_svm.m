function [models, f1, f1_all, h_all] = train_svm(x_train, y_train, x_val, y_val)

% getting the size of the problem
[N, d]  = size(x_train);
[~, K]  = size(y_train);

% getting the regularization parameters
lambdas = logspace(-7, 2, 10);

% storing the f1 loss and the corresponding model
f1 = zeros(K, length(lambdas));
f1_all = zeros(1, length(lambdas));
h_all = zeros(1, length(lambdas));
models = cell(K, length(lambdas));


for j = 1:length(lambdas)
    fprintf('lambda(%3d) = %5.3e\n', j, lambdas(j));
    
    y_p = zeros(size(y_val, 1), K);
    
    for i = 1:K
        % setting the parameters
        C = 2 / (lambdas(j) * N);
        options = sprintf('-B 0 -e 10e-7 -q -c %2.7f', C);
        
        % training the model
        model = train(y_train(:, i), x_train, options);
        
        % predicting
        y_p(:, i) = predict(y_val(:, i), x_val, model);
        
        % evaluating
        f1(i, j) = compute_f1(y_p(:, i), y_val(:, i));
        models{i, j} = model;
    end
    
    f1_all(j) = compute_f1(y_p, y_val);
    h_all(j) = compute_hamming(y_p, y_val);
end

end
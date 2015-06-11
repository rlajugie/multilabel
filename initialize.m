function [W, A] = initialize(x_train, y_train, x_val, y_val, lambda, params)

% getting the size of the problem
[~, d] = size(x_train);
[~, V] = size(y_train);

W = zeros(d, V);

if strcmp(params.init_w, 'svm')
    [root, fname, ~] = fileparts(params.data_path);
    init_fname = fullfile(root, sprintf('%s_init_w.mat', fname));
    
    if exist(init_fname, 'file')
        temp = load(init_fname);
        W = temp.W;
    else
        fprintf('initializing with svms...\n');
        [models, f1_loss] = train_svm(sparse(x_train), y_train, sparse(x_val), y_val);
        
        [~, idx] = min(f1_loss, [], 2);
        for i = 1:V
            if models{i, idx(i)}.Label(1) > 0
                W(:, i) = models{i, idx(i)}.w(1:end-1)';
            else
                W(:, i) = - models{i, idx(i)}.w(1:end-1)';
            end
        end
        save(init_fname, 'W');
    end
else
    W =  ((x_train' * x_train + lambda * eye(d)) \ (x_train')) * y_train;
end

if params.quadratic
    switch params.init_a
        case 'rand'
            A     = -0.000001*(2*rand(V)-1);
            A     = A + A';
        case 'eye'
            A     = eye(V);
    end
else
    A = 0;
end

A = project_A(A, params);

end
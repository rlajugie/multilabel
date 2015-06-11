function [x_train, x_val, x_test, y_train, y_val, y_test] = get_data(params)

% loading the data
load(params.data_path);

% making the data matrix full (in case it is sparse)
x_train = full(x_train);
x_test  = full(x_test);

% adding a constant term to the features to get a bias
x_train = cat(2, x_train,   10 * ones(size(x_train, 1), 1));
x_test  = cat(2, x_test,    10 * ones(size(x_test, 1), 1));

% forcing the labels to be in [-1, 1]
if sum(y_train(:)==0)>0
    y_train = 2*full(y_train)-1;
    y_test  = 2*full(y_test)-1;
end

% splitting the validation set
[Ntrain, ~] = size(x_train);
perm     = randperm(Ntrain);
breaker = round(Ntrain * 0.8);
x_val   = x_train(perm((breaker+1):end), :);
y_val   = y_train(perm((breaker+1):end), :);
x_train = x_train(perm(1:breaker), :);
y_train = y_train(perm(1:breaker), :);

end
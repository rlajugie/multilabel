
mosek_path      = '/sequoia/data1/bojanows/local/mosek/7/toolbox/r2009b';
cvx_path        = '/scratch/bojanows/local/cvx';
gc_path         = '~/Documents/thesis/NIPS2014/submod/Code/graph-cuts/matlab_wrapper';
lr_sdp_path     = '~/Documents/thesis/NIPS2014/submod/Code/low-rank-sdp';
liblinear_path  = '/sequoia/data1/bojanows/local/liblinear-1.96/matlab';

% setting up the path
path_setup;

%%

% method hyperparameters
lambda_w    = 10e-4;
lambda_a    = 10e-2;
g_w         = 1;
g_a         = 1;

% parameter structure
params                  = [];

params.seed             = 1;    % random seed
params.max_trials       = 200;  % number of samples for sdp rounding

params.loss             = 'hamming';    % loss on labelings [f1, hamming]
params.relaxation       = 'graph-cut';   % relaxation type [graph-cut, sdp, spectral]
params.solver           = 'low-rank';   % sdp solver [cvx, mosek, low-rank]

params.data_path        = 'datasets/yeast_dataset.mat'; % path to data

params.T                = 1000000;  % number of subgradient steps
params.time             = 50000;    % compute loss every [time] steps
params.be               = -0.5;     % step-size exponent see function get_stepsize for details [-1, -0.75, -0.5]

params.quadratic        = true;         % use quadratic penalty (A) or not
params.proj_A           = 'negative';   % type of projection for A [positive, negative, none]
params.init_a           = 'rand';       % initialization for A [rand, eye]
params.init_w           = 'rand';       % initialization for W [rand, svm]

params.mosek_license    = '/sequoia/data1/bojanows/local/mosek/7/licenses'; % mosek license file

% launching the subgradient descent
outputs = subgradient_descent(lambda_w, lambda_a, g_w, g_a, params);

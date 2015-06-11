
mosek_path      = '/sequoia/data1/bojanows/local/mosek/7/toolbox/r2009b';
cvx_path        = '';
gc_path         = '~/Documents/thesis/NIPS2014/submod/Code/graph-cuts/matlab_wrapper/';
lr_sdp_path     = '~/Documents/thesis/NIPS2014/submod/Code/low-rank-sdp/';

% setting up the path
path_setup;

% method hyperparameters
lambda_w    = 10e-4;
lambda_a    = 10e-2;
g_w         = 1;
g_a         = 1;
proj_A      = 'positive';

% parameter structure
params                  = [];
params.max_trials       = 200;
params.solver           = 'low-rank';
params.data_path        = '/sequoia/data1/bojanows/NIPS2014/mulan/bibtex_dataset.mat';
params.T                = 10000;
params.time             = 10000;
params.botou            = 1;
params.quadratic        = true;
params.seed             = 1;
params.init_a           = 'rand';
params.init_w           = 'svm';
params.loss             = 'hamming';
params.relaxation       = 'sdp';
params.res_dir          = '/sequoia/data1/bojanows/NIPS2014';
params.experiment       = 'experiment_07_06_2014_16_08_spectral';
params.mosek_license    = '/sequoia/data1/bojanows/local/mosek/7/licenses';
params.save_string      = 'sprintf(''lw%4.2e_la%4.2e_gw%4.2e_ga%4.2e_%s_%s.mat'', lambda_w, lambda_a, g_w, g_a, dataset, proj_A)';

% launching the subgradient descent
subgradient_descent(lambda_w, lambda_a, g_w, g_a, proj_A, params);

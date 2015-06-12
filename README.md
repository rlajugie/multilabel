# Semidefinite and Spectral Relaxations for Multi-Label Classification 

## Description

## Dependencies

The code requires some external libraries depending on the needs.

#### MOSEK

https://www.mosek.com/

#### CVX

http://cvxr.com/cvx/download/

#### Max-flow / min-cut

Matlab Wrapper by Andrew Delong.
http://vision.csd.uwo.ca/wiki/vision/upload/d/d7/Bk_matlab.zip

#### LIBLINEAR

http://www.csie.ntu.edu.tw/~cjlin/liblinear/

#### Low-Rank Optimization on the Cone of Positive Semidefinite Matrices

Implementation by Journee.
http://www.montefiore.ulg.ac.be/~journee/Low_rank_optimization.zip

## Running the code

An example of running is provided in the file main.m.

## Details and parameters

Our method has 4 hyperparameters that have to be (cross-)validated. Additional parameters can also been used.

### Hyperparameters

lambda_w is the regularization parameter for the classifiers w.
lambda_a is the regularization parameter fo


### Parameters

params.seed is the random seed
params.max_trials is the number of samples for sdp rounding (when A can be any matrix, this is important)
params.loss is the loss on labelings and is among [f1, hamming]
params.relaxation is the type of relaxation wanted, it can be[graph-cut, sdp, spectral]
params.solver is the external solver called to solve the SDP when chosen, it can be [cvx, mosek, low-rank]
params.data_path is the path to data
params.T is the number of subgradient steps
params.time compute loss every [time] steps
params.be is the step-size exponent see function get_stepsize for details [-1, -0.75, -0.5]
params.quadratic set the use of quadratic penalty (A) or not (otherwise we are just doing one versus rest SVM using Hamming or F1 loss)
params.proj_A type of affinities in A, can be [positive, negative, none], recall that negative call graph-cut solvers
params.init_a is the initialization for A and can be [rand, eye]
params.init_w is the initialization for W [rand, svm]

